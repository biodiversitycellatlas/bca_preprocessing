#!/usr/bin/env nextflow

/*
==============================================================================
BCA Pre-processing Pipeline
==============================================================================
This pipeline handles the analysis of single-cell RNA sequencing data, including
quality control, demultiplexing, mapping, and filtering.

Pre-requisites:
- Configured the custom config file (config/custom.config)
- Added custom config as profile in the main config file (config/main.config)
- Added profile to the command line option in the submit_nextflow.sh script

Run:
sbatch submit_nextflow.sh main.nf
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

include { SAVE_RUN_CONFIG           } from './modules/custom/save_configs/main'
include { MAPPING_STATS             } from './modules/custom/dashboard/mapping_stats/main'
include { MULTIQC                   } from './modules/tools/multiqc/main'

     
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SETUP CHANNEL 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Set up the sampleID channel from the samplesheet
Channel
  .fromPath(params.input)
  .splitCsv(header: true, sep: ',')
  .map { row ->
    // Concat metadata into a Map object
    def meta = [
      id                : row.sample,
      expected_cells    : row.expected_cells ? row.expected_cells.toInteger() : null,
      p5                : row.p5 ? row.p5 : null,
      p7                : row.p7 ? row.p7 : null,
      rt                : row.rt ? row.rt : null,
    ]
    
    // Assign each FASTQ to its own variable (or null if missing)
    def fastq_cDNA  = row.fastq_cDNA  ? file(row.fastq_cDNA)  : null
    def fastq_BC_UMI = row.fastq_BC_UMI ? file(row.fastq_BC_UMI) : null

    // Return a single map with named entries
    [
        meta         : meta,
        fastq_cDNA   : fastq_cDNA,
        fastq_BC_UMI : fastq_BC_UMI
    ]
  }
  .set { ch_samplesheet }



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
workflow {
    // Save run configurations
    SAVE_RUN_CONFIG()

    // Pre-processing workflow
    preprocessing_workflow(ch_samplesheet)

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
    MULTIQC(mapping_stats_trigger)
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
