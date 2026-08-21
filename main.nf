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
include { post_mapping_workflow     } from './workflows/post_mapping_workflow'
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

        With run_method = "post_mapping" the first two stages are skipped
        entirely: a previous run's published mapping results are read back in
        and everything after mapping is redone on them, which is how the cells
        can be re-called without paying for mapping again.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow BCA_PREPROCESSING {

    take:
        samplesheet     // channel: samplesheet read in from --input

    main:
        // Initialize reporting channels
        def multiqc_report_ch   = Channel.empty()
        def preprocs_output_ch  = Channel.empty()

        // The mapping results the downstream workflows run on, either freshly mapped or
        // read back from a previous run; both branches emit the same channel names
        def mapping_out         = null
        def run_downstream      = false

        // Save run configurations
        SAVE_RUN_CONFIG(samplesheet.first())

        // "post_mapping" resumes a previous run from the point mapping ended, so neither
        // pre-processing nor mapping is repeated
        if (params.run_method == "post_mapping") {

            post_mapping_workflow(samplesheet)
            mapping_out    = post_mapping_workflow.out
            run_downstream = true

        } else {

            // Pre-processing workflow
            preprocessing_workflow(samplesheet)
            preprocs_output_ch = preprocessing_workflow.out.data_output

            // Runs mapping for "standard", "geneext_only"
            if ( params.run_method != "external_pipeline_only" ) {
                // Mapping using STARsolo, Alevin, and/or comparison to commercial pipelines
                QC_mapping_workflow(preprocessing_workflow.out.data_output, preprocessing_workflow.out.bc_whitelist)
                mapping_out = QC_mapping_workflow.out

                // Continue with filtering and MultiQC only with "standard" run_method
                run_downstream = (params.run_method == "standard")
            }
        }

        if (run_downstream) {

            // Filtering raw matrices of ambient RNA
            filtering_workflow(mapping_out.starsolo_genefull50_raw, mapping_out.starsolo_genefull50_filtered, mapping_out.af_mtx, mapping_out.af_filtered_mtx)

            reporting_workflow(
                mapping_out.mapped_samplesheet,
                SAVE_RUN_CONFIG.out.samplesheet,
                mapping_out.ref_gtf,
                SAVE_RUN_CONFIG.out.run_config,
                mapping_out.star_final_log,
                mapping_out.star_summaries,
                mapping_out.star_log,
                mapping_out.starsolo_bam,
                mapping_out.star_solodir,
                mapping_out.starsolo_genefull50_filtered,
                mapping_out.saturation_logs,
                mapping_out.star_cellreads,
                mapping_out.af_meta_info,
                mapping_out.af_quant_json,
                mapping_out.af_cell_meta,
                mapping_out.af_mtx,
                mapping_out.af_umipercell,
                mapping_out.pavian_sankey,
                mapping_out.saturation_imgs,
                mapping_out.saturation_residual_imgs,
                mapping_out.star_umipercell,
                mapping_out.featurecount_txt,
                mapping_out.secondderiv_knee,
                mapping_out.secondderiv_stats,
                mapping_out.secondderiv_cutoff,
                filtering_workflow.out.cs_ambient_hist_plot,
                filtering_workflow.out.cs_umap_comparison_plot,
                filtering_workflow.out.cs_top_genes,
                mapping_out.geneext_report,
                mapping_out.geneext_log
            )

            multiqc_report_ch = reporting_workflow.out.multiqc_report
        }

    emit:
        preprocs_output         = preprocs_output_ch
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
