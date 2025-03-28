#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
==============================================================================
BCA Pre-processing Pipeline
==============================================================================
This pipeline handles the analysis of single-cell RNA sequencing data, including
data downloading, quality control, genome indexing, demultiplexing, mapping, and
SATURATION analysis.

Pre-requisites:
- /fastq/ directory
- Annotation files (FASTA & GTF/GFF) in the /genome/ directory
- Configured the config.yaml file
- Configured the submit_nextflow.sh file

Optional:
- seqspec (.yaml) file

Run:
sbatch submit_nextflow.sh main.nf


Help:
To print the output of processes:
println DEMULTIPLEX.out.splitted_files

------------------------------------------------------------------------------
*/


// Set up the sampleID channel
sample_ids = Channel.fromPath("${params.resDir}/fastq/*_R1_001.fastq.gz")
    .map { file -> 
        file.name.replaceAll(/_R.*_001\.fastq\.gz$/, '')
    }
    .distinct()


// Import sub-workflows
include { parse_workflow } from './workflows/parse_workflow'
include { bd_rhapsody_workflow } from './workflows/bd_rhapsody_workflow'
include { QC_mapping_workflow } from './workflows/QC_mapping_workflow'
include { oak_seq_workflow } from './workflows/oak_seq_workflow'
include { filtering_workflow } from './workflows/filtering_workflow'

// Import processes
include { MULTIQC } from './modules/multiqc'
include { MAPPING_STATS } from './modules/mapping_statistics'


/* 
 * MAIN WORKFLOW
 * 
 * Selects pre-processing workflow depending on the    
 * sequencing technique and returns the pre-processed
 * FASTQ files, and possibly results from the equivalent 
 * commercial pipeline (depending on if the path to the 
 * local installation is given). The pre-processed files
 * are then used for mapping and quality control, and once
 * all outputs are finished, the pipeline triggers MultiQC
 * and the filtering workflow.  
*/                                     
workflow {
    // Sequencing-specific analysis of data:
    //  - Parse Bioscience: Demultiplexing using groups of wells and mapping using split-pipe
    //  - BD Rhapsody: Removing variable bases and mapping using BD rhapsody pipeline
    //  - OAK seq: Mapping using CellRanger
    if (params.seqTech == 'parse_biosciences') {     
        data_output = parse_workflow(sample_ids)
    } else if (params.seqTech == 'bd_rhapsody') {
        data_output = bd_rhapsody_workflow(sample_ids)
    } else if (params.seqTech == 'oak_seq') {
        data_output = oak_seq_workflow(sample_ids)
    } else {
        error "Invalid sequencing technology specified. Use 'parse_biosciences' or 'bd_rhapsody'."
    }
    
    // Mapping using STARsolo, Alevin, and/or comparison to commercial pipelines
    qc_output = QC_mapping_workflow(data_output)

    // Filtering raw matrices of ambient RNA and detecting doublets
    filter_out = filtering_workflow(qc_output.mapping_files)

    // Collect all outputs into a single channel and create trigger
    all_outputs = data_output.mix(qc_output.all_outputs)
    mapping_stats_trigger = all_outputs.collect().map { it -> true }
    
    // MultiQC and mapping statistics, only triggered after all outputs are finished
    MULTIQC(mapping_stats_trigger)
    MAPPING_STATS(mapping_stats_trigger, sample_ids.collect()) 
}


// Runtime Information
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

