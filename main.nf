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
- Accession list file with sample IDs (e.g., accession_lists/species_accessions.txt)
- Annotation files (FASTA & GTF/GFF) in the genome directory

Optional:
- seqspec (.yaml) file

Run:
sbatch submit_nextflow.sh main.nf

==========================
TODO list: 

- enable workflow notification -> sends an email upon completion:
    nextflow run <pipeline name> -N <recipient address>

- incorporate seqspec

==========================
Help:

To print the output of processes:
println DEMULTIPLEX.out.splitted_files

------------------------------------------------------------------------------
*/

// ================     CHANNELS    ================= \\
// Read initial sample IDs from accession list
Channel.fromPath("${params.accessions}")
    .splitText()
    .map { it.trim() }
    .filter { it }
    .set { sample_ids }


// =================     IMPORTS    ================= \\
// Import sub-workflows
include { parse_workflow } from './workflows/parse_workflow'
include { bd_rhapsody_workflow } from './workflows/bd_rhapsody_workflow'
include { QC_mapping_workflow } from './workflows/QC_mapping_workflow'

// Import processes
include { MAPPING_STATS } from './modules/mapping_statistics'


// =================   MAIN WORKFLOW  ================ \\ 
// Selects pre-processing workflow depending on the    \\
// sequencing technique and returns a bam file         \\
// (output alignment step). This will be used in       \\
// downstream processes, which are identical for all   \\
// techniques.                                         \\

workflow {
    // Sequencing-specific analysis of data:
    //  - Parse Bioscience: Demultiplexing using groups and mapping using split-pipe
    //  - BD Rhapsody: Demultiplexing using groups
    if (params.seqTech == 'parse_biosciences') {     
        data_output = parse_workflow(sample_ids)
    } else if (params.seqTech == 'bd_rhapsody') {
        data_output = bd_rhapsody_workflow(sample_ids)
    } else {
        error "Invalid sequencing technology specified. Use 'parse_biosciences' or 'bd_rhapsody'."
    }
    
    // MultiQC and mapping steps using STARsolo including gene extension and testing different configurations
    qc_output = QC_mapping_workflow(data_output)

    // Collect all outputs into a single channel and create trigger to summarize mapping steps into excel format
    all_outputs = data_output.mix(qc_output)
    mapping_stats_trigger = all_outputs.collect().map { it -> true }
    MAPPING_STATS(mapping_stats_trigger) 

    // Filtering raw matrices of ambient RNA and detecting doublets
    // filtering_workflow()
}


// ============  RUNTIME INFORMATION  ============ \\ 
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

