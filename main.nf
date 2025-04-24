#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

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
include { preprocessing_workflow    } from './workflows/preprocessing_workflow.nf'
include { QC_mapping_workflow       } from './workflows/mapping_workflow.nf'
include { filtering_workflow        } from './workflows/filtering_workflow.nf'

include { MULTIQC                   } from './modules/multiqc'
include { MAPPING_STATS             } from './modules/mapping_statistics'

     
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SETUP CHANNEL 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Set up the sampleID channel from the fastq_dir
// sample_ids = Channel.fromPath("${params.fastq_dir}/*_R1_001.fastq.gz")
//     .map { file -> 
//         file.name.replaceAll(/_R.*_001\.fastq\.gz$/, '')
//     }
//     .distinct()

// Set up the sampleID channel from the samplesheet
Channel
  .fromPath(params.input)
  .splitCsv(header: true, sep: ',')
  .map { row ->
    // Concat metadata into a Map object
    def meta = [
      id            : row.sample,
      expected_cells: row.expected_cells ? row.expected_cells.toInteger() : null,
      i5_index      : row.i5_index ? row.i5_index : null,
    ]
    // List of fastq paths
    def fastqs = [ file(row.fastq_cDNA), file(row.fastq_BC_UMI) ].findAll { it } 

    // Emit as a 2-element List
    [ meta, fastqs ]
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
    // Pre-processing workflow
    preprocessing_workflow(ch_samplesheet)
    
    // Mapping using STARsolo, Alevin, and/or comparison to commercial pipelines
    QC_mapping_workflow(preprocessing_workflow.out)

    // Filtering raw matrices of ambient RNA and detecting doublets
    filter_out = filtering_workflow(QC_mapping_workflow.out.mapping_files)

    // Collect all outputs into a single channel and create trigger
    all_outputs = preprocessing_workflow.out.mix(QC_mapping_workflow.out.all_outputs)
    mapping_stats_trigger = all_outputs.collect().map { it -> true }
    
    // MultiQC and mapping statistics, only triggered after all outputs are finished
    MULTIQC(mapping_stats_trigger)
    MAPPING_STATS(mapping_stats_trigger, ch_samplesheet.collect()) 
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
