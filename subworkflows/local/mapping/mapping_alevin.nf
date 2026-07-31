//
// Subworkflow with functionality specific to the workflow 'mapping_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SALMON_SPLICI                 } from '../../../modules/local/tools/salmon_alevin/salmon_splici/main'
include { SALMON_INDEX                  } from '../../../modules/local/tools/salmon_alevin/salmon_index/main'
include { ALEVIN_FRY                    } from '../../../modules/local/tools/salmon_alevin/alevin-fry/main'
include { ALEVIN_QC                     } from '../../../modules/local/tools/alevinQC/main'
include { SECONDDERIV_CELLCALLING_ALEVIN } from '../../../modules/local/custom/dashboard/2nd_deriv/cellcalling_alevin/main'
include { FILTER_MATRICES_ALEVIN        } from '../../../modules/local/custom/dashboard/2nd_deriv/filter_mtx_alevin/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO RUN ALEVIN-FRY MAPPING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow mapping_alevin_workflow {
    take:
        data_output
        bc_whitelist

    main:
        // Initialize reporting channels
        def ch_filtered_mtx       = Channel.empty()
        def ch_secondderiv_umis   = Channel.empty()
        def ch_secondderiv_knee   = Channel.empty()
        def ch_secondderiv_stats  = Channel.empty()
        def ch_secondderiv_cutoff = Channel.empty()

        // Build one Salmon splici and index for all samples
        SALMON_SPLICI(data_output.first())
        SALMON_INDEX(SALMON_SPLICI.out.splici_index_reference)

        // Run Alevin-fry for each sample
        ALEVIN_FRY(
            data_output,
            bc_whitelist,
            SALMON_SPLICI.out.splici_index_reference,
            SALMON_INDEX.out.salmon_index
        )

        // Quality control report for Alevin-fry outputs
        ALEVIN_QC(ALEVIN_FRY.out.mapping_files)

        // Re-call cells on the second-derivative cutoff, if selected.
        if (params.cellfilter_method == "second_derivative") {
            SECONDDERIV_CELLCALLING_ALEVIN(ALEVIN_FRY.out.af_mtx)

            // Join on meta so every sample is filtered on its own cutoff
            def ch_filter_input = ALEVIN_FRY.out.af_mtx
                .join(SECONDDERIV_CELLCALLING_ALEVIN.out.cutoff)

            FILTER_MATRICES_ALEVIN(ch_filter_input)

            ch_filtered_mtx       = FILTER_MATRICES_ALEVIN.out.filtered_matrix
            ch_secondderiv_stats  = FILTER_MATRICES_ALEVIN.out.filter_stats
            ch_secondderiv_umis   = SECONDDERIV_CELLCALLING_ALEVIN.out.umi_per_cell
            ch_secondderiv_knee   = SECONDDERIV_CELLCALLING_ALEVIN.out.json_data
            ch_secondderiv_cutoff = SECONDDERIV_CELLCALLING_ALEVIN.out.cutoff
        }

    emit:
        mapping_files      = ALEVIN_FRY.out.mapping_files
        af_meta_info       = ALEVIN_FRY.out.af_meta_info
        af_quant_json      = ALEVIN_FRY.out.af_quant_json
        af_cell_meta       = ALEVIN_FRY.out.af_cell_meta
        af_mtx             = ALEVIN_FRY.out.af_mtx
        af_filtered_mtx    = ch_filtered_mtx
        af_umipercell      = ch_secondderiv_umis
        secondderiv_knee   = ch_secondderiv_knee
        secondderiv_stats  = ch_secondderiv_stats
        secondderiv_cutoff = ch_secondderiv_cutoff
        qc_reports         = ALEVIN_QC.out.alevinQC_report
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
