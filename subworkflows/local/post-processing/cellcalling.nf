//
// Subworkflows shared by 'mapping_workflow.nf' and 'post_mapping_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SECONDDERIV_CELLCALLING        } from '../../../modules/local/custom/dashboard/2nd_deriv/cellcalling/main'
include { FILTER_MATRICES                } from '../../../modules/local/custom/dashboard/2nd_deriv/filter_mtx/main'
include { SECONDDERIV_CELLCALLING_ALEVIN } from '../../../modules/local/custom/dashboard/2nd_deriv/cellcalling_alevin/main'
include { FILTER_MATRICES_ALEVIN         } from '../../../modules/local/custom/dashboard/2nd_deriv/filter_mtx_alevin/main'
include { SUBSET_VELOCYTO_MATRICES       } from '../../../modules/local/custom/manipulate/subset_velocyto/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO CALL CELLS ON STARSOLO'S RAW MATRIX
        "second_derivative" and "manual_cutoff" both re-call cells on a UMI threshold
        and filter the raw matrix on it; they differ only in where that threshold comes
        from. Every other cellfilter_method leaves the mapper's own cell call standing.

        Whichever call wins, the velocity matrices are then subset to its barcodes, so
        every matrix published for a sample describes the same cells.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow cellcalling_starsolo_workflow {
    take:
        raw_matrix          // channel: [meta, GeneFull_Ex50pAS/raw]
        umi_per_cell        // channel: [meta, UMIperCellSorted.txt]
        cellreads_stats     // channel: [meta, CellReads.stats]
        star_filtered       // channel: [meta, GeneFull_Ex50pAS/filtered], STARsolo's own cell call
        velocyto_raw        // channel: [meta, Velocyto/raw], empty without params.perform_velocity

    main:
        // Initialize reporting channels
        def ch_filtered_mtx       = Channel.empty()
        def ch_secondderiv_knee   = Channel.empty()
        def ch_secondderiv_stats  = Channel.empty()
        def ch_secondderiv_cutoff = Channel.empty()

        if (params.cellfilter_method in ["second_derivative", "manual_cutoff"]) {
            SECONDDERIV_CELLCALLING(umi_per_cell)

            // Join on meta so every sample is filtered on its own cutoff
            def ch_filter_input = raw_matrix
                .join(SECONDDERIV_CELLCALLING.out.cutoff)
                .join(cellreads_stats)

            FILTER_MATRICES(ch_filter_input)

            ch_filtered_mtx           = FILTER_MATRICES.out.filtered_matrix
            ch_secondderiv_stats      = FILTER_MATRICES.out.filter_stats
            ch_secondderiv_knee       = SECONDDERIV_CELLCALLING.out.json_data
            ch_secondderiv_cutoff     = SECONDDERIV_CELLCALLING.out.cutoff
        } else {
            ch_filtered_mtx = star_filtered
        }

        // The velocity matrices take their cell set from whichever GeneFull_Ex50pAS cell call is
        // in force, rather than a cutoff of their own, so every published matrix for a sample
        // describes the same cells. This runs under every cellfilter_method.
        def ch_velocyto_input = velocyto_raw
            .join(ch_filtered_mtx)
            .map { meta, raw_dir, filtered_dir -> [meta, raw_dir, file("${filtered_dir}/barcodes.tsv")] }

        SUBSET_VELOCYTO_MATRICES(ch_velocyto_input)

    emit:
        filtered_matrix    = ch_filtered_mtx
        velocyto_filtered  = SUBSET_VELOCYTO_MATRICES.out.filtered_matrix
        secondderiv_knee   = ch_secondderiv_knee
        secondderiv_stats  = ch_secondderiv_stats
        secondderiv_cutoff = ch_secondderiv_cutoff
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO CALL CELLS ON ALEVIN-FRY'S QUANTIFIED MATRIX
        alevin-fry has no UMI-per-cell curve of its own, so the module derives one from
        the matrix before applying the cutoff. The filter runs after alevin-fry's knee
        rather than replacing it.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow cellcalling_alevin_workflow {
    take:
        af_mtx              // channel: [meta, <sample>_counts/alevin]

    main:
        // Initialize channels
        def ch_af_mtx = af_mtx
        def ch_filtered_mtx       = Channel.empty()
        def ch_secondderiv_umis   = Channel.empty()
        def ch_secondderiv_knee   = Channel.empty()
        def ch_secondderiv_stats  = Channel.empty()
        def ch_secondderiv_cutoff = Channel.empty()

        if (params.cellfilter_method in ["second_derivative", "manual_cutoff"]) {
            SECONDDERIV_CELLCALLING_ALEVIN(ch_af_mtx)

            // Join on meta so every sample is filtered on its own cutoff
            def ch_filter_input = ch_af_mtx
                .join(SECONDDERIV_CELLCALLING_ALEVIN.out.cutoff)

            FILTER_MATRICES_ALEVIN(ch_filter_input)

            ch_filtered_mtx       = FILTER_MATRICES_ALEVIN.out.filtered_matrix
            ch_secondderiv_stats  = FILTER_MATRICES_ALEVIN.out.filter_stats
            ch_secondderiv_umis   = SECONDDERIV_CELLCALLING_ALEVIN.out.umi_per_cell
            ch_secondderiv_knee   = SECONDDERIV_CELLCALLING_ALEVIN.out.json_data
            ch_secondderiv_cutoff = SECONDDERIV_CELLCALLING_ALEVIN.out.cutoff
        }

    emit:
        filtered_matrix    = ch_filtered_mtx
        umi_per_cell       = ch_secondderiv_umis
        secondderiv_knee   = ch_secondderiv_knee
        secondderiv_stats  = ch_secondderiv_stats
        secondderiv_cutoff = ch_secondderiv_cutoff
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
