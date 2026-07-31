#!/usr/bin/env nextflow

/*
==============================================================================
BCA Pre-processing Pipeline
==============================================================================
This pipeline handles the analysis of single-cell RNA sequencing data, including
quality control, demultiplexing, mapping, and reporting.

Pre-requisites:
- Created a samplesheet in CSV format (see conf/example_samplesheet.csv)
- Configured the custom config file (config/custom.config)
- Conda & Nextflow available in base environment

Run:
nextflow run -profile <institution_config>,conda -c /path/to/custom_config
------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { preprocessing_workflow    } from './workflows/preprocessing_workflow'
include { QC_mapping_workflow       } from './workflows/mapping_workflow'
include { filtering_workflow        } from './workflows/filtering_workflow'
include { reporting_workflow        } from './workflows/reporting_workflow'

include { PIPELINE_INITIALISATION   } from './subworkflows/local/utils_nfcore_bca_pipeline'
include { PIPELINE_COMPLETION       } from './subworkflows/local/utils_nfcore_bca_pipeline'

include { SAVE_RUN_CONFIG           } from './modules/local/custom/save_configs/main'
include { MAPPING_STATS             } from './modules/local/custom/dashboard/mapping_stats/main'
include { MULTIQC                   } from './modules/local/tools/multiqc/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
        Selects pre-processing workflow depending on the sequencing technique
        and returns the pre-processed FASTQ files, and possibly results from
        the equivalent commercial pipeline (depending on if the path to the
        local installation is given). The pre-processed files are then used
        for mapping and quality control, and once all outputs are finished,
        the pipeline triggers MultiQC and the filtering workflow.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow BCA_PREPROCESSING {

    take:
        samplesheet     // channel: samplesheet read in from --input

    main:
        // Initialize reporting channels
        def multiqc_report_ch   = Channel.empty()

        // Save run configurations
        SAVE_RUN_CONFIG(samplesheet.first())

        // Pre-processing workflow
        preprocessing_workflow(samplesheet)

        // Runs mapping for "standard", "geneext_only"
        if ( params.run_method != "external_pipeline_only" ) {
            // Mapping using STARsolo, Alevin, and/or comparison to commercial pipelines
            QC_mapping_workflow(preprocessing_workflow.out.data_output, preprocessing_workflow.out.bc_whitelist)
        }

        // Continue with filtering and MultiQC only with "standard" run_method
        if (params.run_method == "standard") {

            // Filtering raw matrices of ambient RNA
            filtering_workflow(QC_mapping_workflow.out.starsolo_genefull50_raw, QC_mapping_workflow.out.starsolo_genefull50_filtered, QC_mapping_workflow.out.af_mtx)

            reporting_workflow(
                QC_mapping_workflow.out.mapped_samplesheet,
                SAVE_RUN_CONFIG.out.samplesheet,
                QC_mapping_workflow.out.ref_gtf,
                SAVE_RUN_CONFIG.out.run_config,
                QC_mapping_workflow.out.star_final_log,
                QC_mapping_workflow.out.star_summaries,
                QC_mapping_workflow.out.star_log,
                QC_mapping_workflow.out.starsolo_bam,
                QC_mapping_workflow.out.star_solodir,
                QC_mapping_workflow.out.saturation_logs,
                QC_mapping_workflow.out.star_cellreads,
                QC_mapping_workflow.out.af_meta_info,
                QC_mapping_workflow.out.af_quant_json,
                QC_mapping_workflow.out.af_cell_meta,
                QC_mapping_workflow.out.af_mtx,
                QC_mapping_workflow.out.pavian_sankey,
                QC_mapping_workflow.out.saturation_imgs,
                QC_mapping_workflow.out.saturation_residual_imgs,
                QC_mapping_workflow.out.star_umipercell,
                QC_mapping_workflow.out.featurecount_txt,
                QC_mapping_workflow.out.secondderiv_knee,
                QC_mapping_workflow.out.secondderiv_stats,
                QC_mapping_workflow.out.secondderiv_cutoff,
                filtering_workflow.out.cs_ambient_hist_plot,
                filtering_workflow.out.cs_umap_comparison_plot,
                filtering_workflow.out.cs_top_genes
            )

            multiqc_report_ch = reporting_workflow.out.multiqc_report
        }

    emit:
        preprocs_output         = preprocessing_workflow.out.data_output
        multiqc_report          = multiqc_report_ch
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    BCA_PREPROCESSING (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        BCA_PREPROCESSING.out.multiqc_report
    )
}



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOWS TO DISPLAY RUNTIME INFORMATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow.onComplete {
    summary = """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
    println summary
}

workflow.onError {
    println "Error: Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
