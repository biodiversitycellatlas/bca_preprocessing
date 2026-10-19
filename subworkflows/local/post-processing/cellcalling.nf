//
// Subworkflow shared by 'mapping_workflow.nf' and 'post_mapping_workflow.nf', for both mappers
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SECONDDERIV_CELLCALLING        } from '../../../modules/local/custom/dashboard/2nd_deriv/cellcalling/main'
include { FILTER_MATRICES                } from '../../../modules/local/custom/dashboard/2nd_deriv/filter_mtx/main'
include { SUBSET_VELOCYTO_MATRICES       } from '../../../modules/local/custom/manipulate/subset_velocyto/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO CALL CELLS ON A MAPPER'S MATRIX, FOR STARSOLO AND ALEVIN-FRY ALIKE
        "second_derivative" and "manual_cutoff" both re-call cells on a UMI threshold
        and filter the matrix on it; they differ only in where that threshold comes
        from. Every other cellfilter_method leaves the mapper's own cell call standing.

        The two mappers differ only in what they hand over:
          - STARsolo: its raw matrix, its own UMI-per-cell curve, CellReads.stats for the
            read-level statistics, its own filtered matrix, and the Velocyto matrices.
          - alevin-fry: its quantified matrix only. It has no UMI-per-cell curve of its own,
            so the module derives one from the matrix. alevin-fry keeps every whitelisted
            barcode (an unfiltered permit list), so the cutoff is found on the whole curve,
            as for STARsolo; without a usable whitelist it runs after alevin-fry's knee.
            Its velocity object is built from the USA matrix downstream, so there is no
            Velocyto subset, and no cell call of its own to fall back on here.

        Whichever call wins, the velocity matrices are then subset to its barcodes, so
        every matrix published for a sample describes the same cells.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow cellcalling_workflow {
    take:
        mapper              // value: 'starsolo' or 'alevin-fry'
        matrix              // channel: [meta, GeneFull_Ex50pAS/raw | <sample>_counts/alevin], the matrix the cutoff filters
        umi_curve           // channel: [meta, UMIperCellSorted.txt | <sample>_counts/alevin], what the curve is read or derived from
        cellreads_stats     // channel: [meta, CellReads.stats], empty for alevin-fry
        mapper_filtered     // channel: [meta, GeneFull_Ex50pAS/filtered], the mapper's own cell call, empty for alevin-fry
        velocyto_raw        // channel: [meta, Velocyto/raw], empty for alevin-fry and without params.perform_velocity

    main:
        // Initialize reporting channels
        def ch_filtered_mtx       = channel.empty()
        def ch_umi_per_cell       = channel.empty()
        def ch_secondderiv_knee   = channel.empty()
        def ch_secondderiv_stats  = channel.empty()
        def ch_secondderiv_cutoff = channel.empty()
        def ch_versions           = channel.empty()

        if (params.cellfilter_method in ["second_derivative", "manual_cutoff"]) {
            SECONDDERIV_CELLCALLING(umi_curve.map { meta, curve -> [meta, mapper, curve] })

            // Join on meta so every sample is filtered on its own cutoff. alevin-fry has no
            // CellReads.stats, so it is given an empty slot rather than joined on an empty channel.
            def ch_with_cutoff = matrix.join(SECONDDERIV_CELLCALLING.out.cutoff)
            def ch_filter_input = mapper == 'starsolo'
                ? ch_with_cutoff.join(cellreads_stats)
                : ch_with_cutoff.map { meta, mtx_dir, cutoff -> [meta, mtx_dir, cutoff, []] }

            FILTER_MATRICES(ch_filter_input.map { meta, mtx_dir, cutoff, cellreads -> [meta, mapper, mtx_dir, cutoff, cellreads] })

            ch_filtered_mtx       = FILTER_MATRICES.out.filtered_matrix
            ch_secondderiv_stats  = FILTER_MATRICES.out.filter_stats
            ch_umi_per_cell       = SECONDDERIV_CELLCALLING.out.umi_per_cell
            ch_secondderiv_knee   = SECONDDERIV_CELLCALLING.out.json_data
            ch_secondderiv_cutoff = SECONDDERIV_CELLCALLING.out.cutoff
            ch_versions           = ch_versions.mix(SECONDDERIV_CELLCALLING.out.versions, FILTER_MATRICES.out.versions)
        } else {
            ch_filtered_mtx = mapper_filtered
        }

        // The velocity matrices take their cell set from whichever GeneFull_Ex50pAS cell call is
        // in force, rather than a cutoff of their own, so every published matrix for a sample
        // describes the same cells. This runs under every cellfilter_method.
        def ch_velocyto_input = velocyto_raw
            .join(ch_filtered_mtx)
            .map { meta, raw_dir, filtered_dir -> [meta, raw_dir, file("${filtered_dir}/barcodes.tsv")] }

        SUBSET_VELOCYTO_MATRICES(ch_velocyto_input)
        ch_versions = ch_versions.mix(SUBSET_VELOCYTO_MATRICES.out.versions)

    emit:
        filtered_matrix    = ch_filtered_mtx
        velocyto_filtered  = SUBSET_VELOCYTO_MATRICES.out.filtered_matrix
        umi_per_cell       = ch_umi_per_cell        // alevin-fry's derived curve; STARsolo's comes from the mapper
        secondderiv_knee   = ch_secondderiv_knee
        secondderiv_stats  = ch_secondderiv_stats
        secondderiv_cutoff = ch_secondderiv_cutoff
        versions           = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
