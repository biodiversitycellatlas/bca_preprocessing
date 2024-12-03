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

==========================
TODO list: 

- enable workflow notification -> sends an email upon completion:
    nextflow run <pipeline name> -N <recipient address>
- publishDir: , overwrite: false 
- save files somewhere outside work directory, now all files are symlink to /work/

- check why 1st splitted fastq files are empty for Nvec
- incorporate seqspec
- split_parse_data: save unused reads to file and inspect

==========================
Help:

To print the output of processes:
println DEMULTIPLEX.out.splitted_files

==========================
Last run:
Exp: 240810_ParseBio_Nvec_Tcas
- Nvec_Integr: 
- Tcas: download_data, parse index, starsolo ref -- 2fee931e-063f-4b8a-9ea2-6aa39c3641a3

Exp: 241106_BD_Rhapsody_Nvec
- Nvec: started (a373b83b-b4de-4e25-908b-75506019f904)

Currently running:
Exp: 240810_ParseBio_Nvec_Tcas
- Nvec_Integr: 
- Tcas: resume from parse & starsolo index

Exp: 241106_BD_Rhapsody_Nvec
- Nvec: from start 

------------------------------------------------------------------------------
*/

// ============== Define CUSTOM parameters ============== \\

// PARSE -- NVEC !!
// params.experiment = "240810_ParseBio_Nvec_Tcas"
// params.species = "Nvec"
// params.seqTech = "parse_biosciences"
// params.dataDir = "/users/asebe/bvanwaardenburg/git/bca_preprocessing/data"
// params.resDir = "${params.dataDir}/${params.experiment}/${params.species}_testIntegr"
// params.ref_star_gtf = "${params.resDir}/genome/Nvec_v5_merged_annotation_sort.gtf"
// params.ref_parse_gtf = "${params.resDir}/genome/Nvec_v4_merged_annotation_parse_sort.gtf"
// params.ref_fasta = "${params.resDir}/genome/Nvec_vc1.1_gDNA_mtDNA.fasta"

// params.barcodeDir = "${params.baseDir}/seq_techniques/${params.seqTech}/bc_data_n26_R1_v3_4 \
//                      ${params.baseDir}/seq_techniques/${params.seqTech}/bc_data_v1 \
//                      ${params.baseDir}/seq_techniques/${params.seqTech}/bc_data_R3_v3"


// PARSE -- TCAS !!
params.outDir = "240810_ParseBio_Nvec_Tcas/Tcas"
params.seqTech = "parse_biosciences"
params.dataDir = "/users/asebe/bvanwaardenburg/git/data"
params.accessions = "/users/asebe/bvanwaardenburg/git/data/accession_lists/Tcas_parse_biosciences_accessions.txt" // optional
params.ref_star_gtf = "/users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Tcas/genome/genomic.gtf"
params.ref_parse_gtf = "/users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Tcas/genome/genomic.gtf" // optional, otherwise state same file
params.ref_fasta = "/users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Tcas/genome/GCF_031307605.1_icTriCast1.1_genomic.fna"

// should be optional, if not given check for sequencing technique and automatically use them
params.barcodeDemux = "./seq_techniques/${params.seqTech}/bc_data_n26_R1_v3_4.csv"
params.barcodeDir = "./seq_techniques/${params.seqTech}/bc_data_n26_R1_v3_4 \
                     ./seq_techniques/${params.seqTech}/bc_data_v1 \
                     ./seq_techniques/${params.seqTech}/bc_data_R3_v3"

// BD RHAPSODY -- NVEC !!
// params.experiment = "241106_BD_Rhapsody_Nvec"
// params.species = "Nvec"
// params.seqTech = "bd_rhapsody"
// params.resDir = "${params.dataDir}/${params.experiment}"
// params.ref_star_gtf = "${params.resDir}/genome/Nvec_v5_merged_annotation_sort.gtf"
// params.ref_parse_gtf = "${params.resDir}/genome/Nvec_v4_merged_annotation_parse_sort.gtf"
// params.ref_fasta = "${params.resDir}/genome/Nvec_vc1.1_gDNA_mtDNA.fasta"

// params.barcodeDir = "./seq_techniques/${params.seqTech}/BD_CLS1.txt \
//                      ./seq_techniques/${params.seqTech}/BD_CLS2.txt \
//                      ./seq_techniques/${params.seqTech}/BD_CLS3.txt"



// ============== DEFINE BASE PARAMETERS ============== \\
params.resDir = "${params.dataDir}/${params.outDir}_Integr" 

params.star_config = "./seq_techniques/${params.seqTech}/config_${params.seqTech}_starsolo.txt"
params.star_config_CR = "./seq_techniques/${params.seqTech}/config_${params.seqTech}_starsolo_CR.txt"
params.star_config_CRED = "./seq_techniques/${params.seqTech}/config_${params.seqTech}_starsolo_CRED.txt"

params.barcodeDoublet = params.barcodeDir.split()[0]  // should be a tsv file (or 1 column with barcodes)


// ====================  CHANNELS ===================== \\
// Read initial sample IDs from accession list
Channel.fromPath("${params.accessions}")
    .splitText()
    .map { it.trim() }
    .filter { it }
    .set { sample_ids }


// ================= IMPORT WORKFLOWS ================ \\
include { parse_workflow } from './workflows/parse_workflow'
include { bd_rhapsody_workflow } from './workflows/bd_rhapsody_workflow'

// Import processes
include { FASTQC } from './modules/fastqc'
include { MULTIQC } from './modules/multiqc'
include { GENINDEX_STARSOLO } from './modules/genindex_starsolo'
include { REINDEX_STARSOLO } from './modules/genindex_starsolo'
include { MAPPING_STARSOLO as MAPPING_STARSOLO_N } from './modules/mapping_starsolo'
include { MAPPING_STARSOLO as MAPPING_STARSOLO_CR } from './modules/mapping_starsolo'
include { MAPPING_STARSOLO as MAPPING_STARSOLO_CRGE } from './modules/mapping_starsolo'
include { MAPPING_STARSOLO as REMAPPING_STARSOLO } from './modules/mapping_starsolo'
include { SATURATION as SATURATION_N } from './modules/saturation'
include { SATURATION as SATURATION_CR } from './modules/saturation'
include { SATURATION as SATURATION_CRGE } from './modules/saturation'
include { GENE_EXT } from './modules/gene_ext'
include { INDEX_BAM } from './modules/index_bam'
include { DOUBLET_DET } from './modules/doublet_det'


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
    GENINDEX_STARSOLO(params.ref_star_gtf)

    // Mapping: standard configuration
    MAPPING_STARSOLO_N(data_output, GENINDEX_STARSOLO.out, file(params.star_config), 'N')
    SATURATION_N(MAPPING_STARSOLO_N.out)

    // Mapping: CR-like
    MAPPING_STARSOLO_CR(data_output, GENINDEX_STARSOLO.out, file(params.star_config_CRED), 'CR')
    SATURATION_CR(MAPPING_STARSOLO_CR.out)

    // Mapping: CR-like + Gene extension
    MAPPING_STARSOLO_CRGE(data_output, GENINDEX_STARSOLO.out, file(params.star_config_CR), 'CRGE')
    GENE_EXT(MAPPING_STARSOLO_CRGE.out)
    REINDEX_STARSOLO(GENE_EXT.out)
    REMAPPING_STARSOLO(data_output, REINDEX_STARSOLO.out, file(params.star_config_CRED), 'remappedCRGE')
    remapped_output = REMAPPING_STARSOLO.out
    SATURATION_CRGE(remapped_output)

    // Downstream processes (continuing with config 3)
    INDEX_BAM(remapped_output)
    DOUBLET_DET(remapped_output)
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

