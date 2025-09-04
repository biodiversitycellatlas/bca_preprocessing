#!/usr/bin/env nextflow

/*
==============================================================================
BCA Pre-processing Pipeline
==============================================================================
This pipeline handles the analysis of single-cell RNA sequencing data, including
quality control, demultiplexing, mapping, and filtering.

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
    samplesheet // channel: samplesheet read in from --input

    main:
    // Save run configurations
    SAVE_RUN_CONFIG(samplesheet)

    // Pre-processing workflow
    preprocessing_workflow(samplesheet)

    // Mapping using STARsolo, Alevin, and/or comparison to commercial pipelines
    QC_mapping_workflow(
        preprocessing_workflow.out.data_output,
        preprocessing_workflow.out.bc_whitelist
    )

    // Filtering raw matrices of ambient RNA and detecting doublets
    filter_out = filtering_workflow(QC_mapping_workflow.out.mapping_files)

    // Collect all outputs into a single channel and create trigger
    all_outputs = preprocessing_workflow.out.data_output.mix(QC_mapping_workflow.out.all_outputs)
    
    // Define a trigger that waits for both mapping_files and all_outputs
    mapping_stats_trigger = all_outputs.mix(QC_mapping_workflow.out.mapping_files).collect().map { it -> true }
    
    // MultiQC and mapping statistics, only triggered after all outputs are finished
    MAPPING_STATS(mapping_stats_trigger)
    ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    MULTIQC(mapping_stats_trigger, ch_multiqc_config)

    emit:
    multiqc_report = MULTIQC.out
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
