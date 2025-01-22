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

// ====================  CHANNELS ===================== \\
// Read initial sample IDs from accession list
Channel.fromPath("${params.accessions}")
    .splitText()
    .map { it.trim() }
    .filter { it }
    .set { sample_ids }


// ================= IMPORTS ================ \\
include { parse_workflow } from './workflows/parse_workflow'
include { bd_rhapsody_workflow } from './workflows/bd_rhapsody_workflow'

// Import processes
include { FASTQC } from './modules/fastqc'
include { MULTIQC } from './modules/multiqc'

include { GENINDEX_STARSOLO as GENINDEX_STARSOLO_N } from './modules/genindex_starsolo'
include { GENINDEX_STARSOLO as GENINDEX_STARSOLO_CR } from './modules/genindex_starsolo'
include { GENINDEX_STARSOLO as REINDEX_STARSOLO_N } from './modules/genindex_starsolo'
include { GENINDEX_STARSOLO as REINDEX_STARSOLO_CR } from './modules/genindex_starsolo'

include { MAPPING_STARSOLO as MAPPING_STARSOLO_N } from './modules/mapping_starsolo'
include { MAPPING_STARSOLO as MAPPING_STARSOLO_CR } from './modules/mapping_starsolo'
include { MAPPING_STARSOLO as MAPPING_STARSOLO_CRGE } from './modules/mapping_starsolo'
include { MAPPING_STARSOLO as REMAPPING_STARSOLO_N } from './modules/mapping_starsolo'
include { MAPPING_STARSOLO as REMAPPING_STARSOLO_CR } from './modules/mapping_starsolo'

include { INDEX_BAM as INDEX_BAM_N } from './modules/index_bam'
include { INDEX_BAM as INDEX_BAM_CR } from './modules/index_bam'
include { INDEX_BAM as INDEX_BAM_CRGE } from './modules/index_bam'
include { INDEX_BAM as INDEX_BAM_NGE } from './modules/index_bam'

include { SATURATION as SATURATION_N } from './modules/saturation'
include { SATURATION as SATURATION_CR } from './modules/saturation'
include { SATURATION as SATURATION_CRGE } from './modules/saturation'
include { SATURATION as SATURATION_NGE } from './modules/saturation'

include { GENE_EXT } from './modules/gene_ext'
include { MAPPING_STATS } from './modules/mapping_statistics'


// =================  MAIN WORKFLOW  ================= \\ 
// Selects pre-processing workflow depending on the    \\
// sequencing technique and returns a bam file         \\
// (output alignment step). This will be used in       \\
// downstream processes, which are identical for all   \\
// techniques.                                         \\


workflow {
    if (params.seqTech == 'parse_biosciences') {     
        parse_workflow(sample_ids)
        data_output = parse_workflow.out
    
    } else if (params.seqTech == 'bd_rhapsody') {
        bd_rhapsody_workflow(sample_ids)
        data_output = bd_rhapsody_workflow.out
    
    } else {
        error "Invalid sequencing technology specified. Use 'parse_biosciences' or 'bd_rhapsody'."
    }
    
    // Quality Control
    FASTQC(data_output)
    MULTIQC(FASTQC.out.collect())

    // Mapping STARsolo
    // GENINDEX_STARSOLO_N(params.ref_star_gtf, file(params.star_config_mkref_N), 'N')
    // GENINDEX_STARSOLO_CR(params.ref_star_gtf, file(params.star_config_mkref_CR), 'CR')

    // Mapping: standard configuration
    // MAPPING_STARSOLO_N(data_output, GENINDEX_STARSOLO_N.out, file(params.star_config_ED), params.barcodeDir, 'N')
    // INDEX_BAM_N(MAPPING_STARSOLO_N.out)
    // SATURATION_N(MAPPING_STARSOLO_N.out, INDEX_BAM_N.out)
    // CALC_MT_RRNA_N(MAPPING_STARSOLO_N.out)

    // GENE_EXT(MAPPING_STARSOLO_N.out) //INDEX_BAM_N.out
    // REINDEX_STARSOLO_N(GENE_EXT.out, file(params.star_config_mkref_N), 'N')
    // REMAPPING_STARSOLO_N(data_output, REINDEX_STARSOLO_N.out, file(params.star_config_ED), params.barcodeDir, 'remappedNGE')
   
    // INDEX_BAM_NGE(REMAPPING_STARSOLO_N.out)
    // SATURATION_NGE(REMAPPING_STARSOLO_N.out, INDEX_BAM_NGE.out)

    // // Mapping: CR-like + Gene extension
    // MAPPING_STARSOLO_CR(data_output, GENINDEX_STARSOLO_CR.out, file(params.star_config_CRED), params.barcodeDemux, 'CR')
    // INDEX_BAM_CR(MAPPING_STARSOLO_CR.out)
    // SATURATION_CR(MAPPING_STARSOLO_CR.out, INDEX_BAM_CR.out)
    // CALC_MT_RRNA_CR(MAPPING_STARSOLO_CR.out)

    // GENE_EXT_CR(MAPPING_STARSOLO_CR.out)
    // REINDEX_STARSOLO_CR(GENE_EXT_CR.out, file(params.star_config_mkref_CR), 'CR')
    // REMAPPING_STARSOLO_CR(data_output, REINDEX_STARSOLO_CR.out, file(params.star_config_CRED), params.barcodeDemux, 'remappedCRGE')
 
    // INDEX_BAM_CRGE(REMAPPING_STARSOLO_CR.out)
    // SATURATION_CRGE(REMAPPING_STARSOLO_CR.out, INDEX_BAM_CRGE.out)

    // Downstream processes (continuing with config 1)
    // MAPPING_STATS() 
    // DOUBLET_DET(mapping_output)
}


// ==================  RUNTIME INFORMATION  =================== \\ 
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

