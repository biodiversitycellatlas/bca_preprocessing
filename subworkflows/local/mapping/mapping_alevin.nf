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
include { cellcalling_workflow          } from '../post-processing/cellcalling'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO RUN ALEVIN-FRY MAPPING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow mapping_alevin_workflow {
    take:
        data_output
        bc_whitelist
        ref_gtf
        ref_fasta
        remap_geneext

    main:
        // The GeneExt re-mapping builds its own splici and index, kept apart in genome/
        def index_suffix = remap_geneext == 'true' ? '_geneext' : ''

        // Build one Salmon splici and index for all samples, on the same reference STARsolo uses
        SALMON_SPLICI(data_output.first(), ref_gtf, ref_fasta, index_suffix)
        SALMON_INDEX(SALMON_SPLICI.out.splici_index_reference, index_suffix)

        // Run Alevin-fry for each sample
        ALEVIN_FRY(
            data_output,
            bc_whitelist,
            SALMON_SPLICI.out.splici_index_reference,
            SALMON_INDEX.out.salmon_index,
            MappingRateCheck.flagDir(workflow).toString()
        )

        // Quality control report for Alevin-fry outputs
        ALEVIN_QC(ALEVIN_FRY.out.mapping_files)

        // Re-call cells on a UMI cutoff, if selected. The matrix is both what the cutoff filters
        // and what the UMI curve is derived from; alevin-fry has no CellReads.stats, no cell call
        // of its own to fall back on, and no Velocyto matrices.
        cellcalling_workflow(
            'alevin-fry',
            ALEVIN_FRY.out.af_mtx,
            ALEVIN_FRY.out.af_mtx,
            channel.empty(),
            channel.empty(),
            channel.empty()
        )

        def ch_versions = SALMON_SPLICI.out.versions.mix(
            SALMON_INDEX.out.versions,
            ALEVIN_FRY.out.versions,
            ALEVIN_QC.out.versions,
            cellcalling_workflow.out.versions
        )

    emit:
        mapping_files      = ALEVIN_FRY.out.mapping_files
        af_meta_info       = ALEVIN_FRY.out.af_meta_info
        af_quant_json      = ALEVIN_FRY.out.af_quant_json
        af_cell_meta       = ALEVIN_FRY.out.af_cell_meta
        af_mtx             = ALEVIN_FRY.out.af_mtx
        af_filtered_mtx    = cellcalling_workflow.out.filtered_matrix
        af_umipercell      = cellcalling_workflow.out.umi_per_cell
        secondderiv_knee   = cellcalling_workflow.out.secondderiv_knee
        secondderiv_stats  = cellcalling_workflow.out.secondderiv_stats
        secondderiv_cutoff = cellcalling_workflow.out.secondderiv_cutoff
        qc_reports         = ALEVIN_QC.out.alevinQC_report
        versions           = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
