#!/usr/bin/env nextflow

/*
==============================================================================
BCA Pre-processing Pipeline
==============================================================================
This pipeline handles the analysis of single-cell RNA sequencing data, including
data downloading, quality control, genome indexing, demultiplexing, mapping, and
saturation analysis.

Usage:
nextflow run preprocs_pipeline.nf

Pre-requisites:
- Accession list file with sample IDs (e.g., accession_lists/species_accessions.txt)
- Annotation files (FASTA & GTF/GFF) in the genome directory
- Barcodes (if using Parse Biosciences technology)

Optional:
- seqspec (.yaml) file

Folder structure: 
- code/
    - integrated_pipe/
        - .nf
        - .config
        - . 
        - process_configs/
- data/
    - accession_lists/
        - {species}_accessions.txt
    - {species}/
        - fasta/
        - genome/
            - .fasta
            - .gtf OR .gff
        - barcodes/
            - .csv


==========================
TODO list: 

- enable workflow notification -> sends an email upon completion:
    nextflow run <pipeline name> -N <recipient address>
- publishDir: , overwrite: false 
- save files somewhere outside work directory, now all files are symlink to /work/

- check why 1st splitted fastq files are empty for Nvec
- error calculating sjdboverhang 

==========================
Help:

To print the output of processes:
println splitting_parse.out.splitted_files

==========================
Last run:
Exp: 240810_ParseBio_Nvec_Tcas
- Nvec: -
- Tcas: downloading -> multiqc

Currently running:
Exp: 240810_ParseBio_Nvec_Tcas
- Nvec: up until multiqc from start - gzip in split_parse
- Tcas: (( next up )) resume: trying to fix mkref parse

------------------------------------------------------------------------------
*/

// Define parameters

// PARSE -- NVEC !!
// params.dataDir = "/users/asebe/bvanwaardenburg/git/bca_preprocessing/data"
// params.codeDir = "/users/asebe/bvanwaardenburg/git/bca_preprocessing/code/integrated_pipe"
// params.barcodeDir = "/users/asebe/bvanwaardenburg/ParseBiosciences-Pipeline.1.3.1/splitpipe/barcodes/"
// params.experiment = "240810_ParseBio_Nvec_Tcas"
// params.species = "Nvec"
// params.seqTech = "parse"
// params.resDir = "${params.dataDir}/${params.experiment}/${params.species}_testIntegr"
// params.ref_star_gtf = "${params.resDir}/genome/Nvec_v5_merged_annotation_sort.gtf"
// params.ref_parse_gtf = "${params.resDir}/genome/Nvec_v4_merged_annotation_parse_sort.gtf"
// params.ref_fasta = "${params.resDir}/genome/Nvec_vc1.1_gDNA_mtDNA.fasta"

// PARSE -- TCAS !!
params.dataDir = "/users/asebe/bvanwaardenburg/git/bca_preprocessing/data"
params.codeDir = "/users/asebe/bvanwaardenburg/git/bca_preprocessing/code/integrated_pipe"
params.barcodeDir = "/users/asebe/bvanwaardenburg/ParseBiosciences-Pipeline.1.3.1/splitpipe/barcodes/"
params.experiment = "240810_ParseBio_Nvec_Tcas"
params.species = "Tcas"
params.seqTech = "parse"
params.resDir = "${params.dataDir}/${params.experiment}/${params.species}"
params.ref_star_gtf = "${params.resDir}/genome/genomic.gtf"
params.ref_parse_gtf = "${params.resDir}/genome/genomic.gtf"
params.ref_fasta = "${params.resDir}/genome/GCF_031307605.1_icTriCast1.1_genomic.fna"


// Define channels
// Read initial sample IDs from accession list
Channel.fromPath("${params.dataDir}/accession_lists/${params.species}_accessions.txt")
    .splitText()
    .map { it.trim() }
    .filter { it }
    .set { sample_ids }

// Workflow definition
workflow {
    download_data(sample_ids)           // Download data
    splitting_parse(download_data.out)  // Demultiplexing
    
    /* Function that retrieves the output from each splitting_parse process 
       (1 parse sample -> multiple splitted files) and flattens it, 
       so the next process (fastqc) only gets one file name at a time. */
    def single_fastqs = splitting_parse.out.splitted_files
        .flatten()
        .map { it.toString() }
    
    fastqc(single_fastqs)               // Quality control with FastQC 
    multiqc(fastqc.out.collect())       // Generate MultiQC report - collect(): requires the fastqc step to be finished before proceeding
    
    // genome_index_starsolo()          // Generate genome index (only once)
    ref_genome_parse()               // Generate reference genome (only once)
    
//     // Function to handle splitting_parse output
//     def paired_fastqs = splitting_parse.out.splitted_files
//         .flatten()
//         .groupTuple(2) // Group the flattened output into pairs

//     mapping_STARsolo(splitting_parse.out, genome_index_starsolo.out)    // Mapping with STARsolo
//     mapping_PB(splitting_parse.out, ref_genome_parse.out)               // Mapping with Parse Biosciences pipeline
//     saturation(mapping_PB.out)                                          // Saturation analysis
}


// ===================== Processes ===================== //

// =================  DOWNLOADING  ================= \\ 
// Check if the user provided a fastq file with the  \\
// raw data (located in a /fastq/ folder) or if it   \\
// still needs to be downloaded. In that case, an    \\
// accession file must be provided with IDs from the \\
// SRA archive. If the data folder is given, an      \\
// accession list will be created for looping over   \\
// the samples during the workflow.                  \\
process download_data {
    tag "${sample_id}"
    debug true
    
    input:
    val sample_id

    output:
    tuple val(sample_id), path("${sample_id}_*.fastq.gz")

    script:
    """
    # Check if FASTQ files already exist
    if [ -f "${params.resDir}/fastq/${sample_id}_R1_001.fastq.gz" ]; then
        echo "FASTQ files for sample ${sample_id} already exist."
        ln -s "${params.resDir}/fastq/${sample_id}_R1_001.fastq.gz" .
        ln -s "${params.resDir}/fastq/${sample_id}_R2_001.fastq.gz" .
    else
        echo "Downloading data for sample ${sample_id}"
        # prefetch ${sample_id}
        # fastq-dump --split-files --gzip --outdir /${params.resDir}/fastq/ ${sample_id}
        # ln -s "/${params.resDir}/fastq/" .
    fi
    """
}

// ==================  DEMULTIPLEXING  ================== \\ 
// To split the fastq files of each library into separate \\
// fastq files for each fixation method, a script for     \\
// demultiplexing the reads is called. From a txt file    \\
// 'sample_wells', wells associated with each fixation    \\
// method are given. Based on the provided barcode file,  \\
// the samples are linked to the barcodes, and splitted   \\
// into n seperate fastq's. This step is repeated for     \\
// all libraries. 
// TODO: The unassigned reads are then saved in a separate \\
// file for manual inspection.                            \\
// only do this for the parse data                        \\
process splitting_parse {
    tag "${sample_id}"
    label "big_cpus"
    publishDir "${params.resDir}/demux_fastq", mode: 'symlink'
    debug true
    maxForks 4
    
    input:
    tuple val(sample_id), path(fastq_files)

    output:
    path("*.fastq.gz"), emit: splitted_files

    script:
    """
    # Verify the input files
    echo "Processing sample ${sample_id}"
    echo "Input files: ${fastq_files}"
    ls -lh ${fastq_files}

    # Run the demux code
    bash ${params.codeDir}/split_parse_data.sh ${params.resDir} ${sample_id} "${fastq_files}" ${params.barcodeDir} 
    """
}

// ==================  FASTQC  ================== \\ 
// Generate Quality Control reports using FASTQC  \\
process fastqc {
    publishDir "${params.resDir}/fastqc", mode: 'symlink'
    debug true

    input:
    path(single_fastq)

    output:
    path("*_fastqc.*")

    script:
    """
    echo "Running FastQC for ${single_fastq}"
    echo "Path: ${single_fastq}"

    fastqc ${single_fastq} 
    """
}

// ==================  MULTI-QC  ================== \\ 
// A way of combining seperate FASTQC reports into  \\
// a single analysis. Provides an overview of the   \\
// data, including checks which of the files passed \\
// the quality metrics.                             \\
process multiqc {
    publishDir "${params.resDir}/fastqc", mode: 'symlink'
    debug true
    
    input:
    path(fastqc_reports)

    output:
    file("multiqc_report.html")

    script:
    """
    echo "Running MultiQC for: ${fastqc_reports.join(' ')}"
    multiqc ${fastqc_reports.join(' ')} --outdir . --filename multiqc_report.html
    """
}

// =============  GENOME INDEX STARSOLO  ============= \\ 
// Creates the Genome Indeces required for mapping.    \\
// Two folders will be created within /genome/,        \\
// genome_index: contains genome index for STARsolo    \\
// parse_refgenome: see next process                   \\
process genome_index_starsolo {   
    publishDir "${params.resDir}/genome/genome_index", mode: 'symlink'
    debug true

    output:
    path("*")     
       
    script:
    """
    # Retrieve the first accession number
    first_accs=\$(head -1 ${params.dataDir}/accession_lists/${params.species}_accessions.txt)
    first_fastq="${params.resDir}/fastq/\${first_accs}_1.fastq.gz"  

    # Calculate SJDB overhang using the first read from the first fastq file
    sjdb_overhang=\$(zcat "\${first_fastq}" 2>/dev/null | awk 'NR==2 {print length(\$0)-1; exit}' || echo "") 

    echo "Generating genome index with STAR"
    STAR --runMode genomeGenerate \\
        --genomeFastaFiles ${params.ref_fasta} \\
        --sjdbGTFfile ${params.ref_star_gtf} \\
        --sjdbOverhang \${sjdb_overhang} \\
        --genomeSAindexNbases 12
    """
}

// =============  REF GENOME PARSE PIPELINE  ============= \\ 
// Creates the Genome Indeces required for mapping.        \\
// parse_refgenome: contains the reference genome for      \\
// the Parse Biosciences (split-pipe) pipeline.            \\
process ref_genome_parse {
    publishDir "${params.resDir}/genome/parse_refgenome", mode: 'symlink'
    debug true

    output:
    path("*")

    script:
    """
    split-pipe -m mkref \\
        --genome_name ${params.species} \\
        --genes ${params.ref_parse_gtf} \\
        --fasta ${params.ref_fasta} \\
        --output_dir .
    """
}


// ==================  MAPPING STARSOLO  ================== \\ 
// Mapping using STARsolo
process mapping_STARsolo {    
    publishDir "${params.resDir}/mapping_STARsolo/${sample_id}", mode: 'symlink'
    debug true

    input:
    val sample_id
    val(sample_id), file(demux_fastq_files) 
    file genome_index_files 
    file config_file from "${params.codeDir}/config_${params.seqTech}_starsolo.txt"

    output:
    val(sample_id), file("${sample_id}*")

    script:
    """
    echo "Mapping sample ${sample_id} with STARsolo"

    # Read configuration file
    CONFIG_OPTIONS=\$(cat ${config_file})

    # Mapping step and generating count matrix using STAR
    STAR \\
        --genomeDir genome_index_files \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${sample_id}_ \\
        --readFilesIn ${demux_fastq_files.join(' ')} \\
        --soloCBwhitelist ${params.barcodeDir} \\
        \${CONFIG_OPTIONS}

    # Create the index file
    samtools index ${sample_id}_Aligned.sortedByCoord.out.bam
    """
}

// ==============  MAPPING PARSE PIPELINE  ============== \\ 
// Mapping using ParseBiosciences pipeline
process mapping_PB {
    tag "${sample_id}"
    publishDir "${params.resDir}/mapping_parsepipe/${sample_id}", mode: 'symlink'
    debug true

    input:
    val sample_id
    path demux_fastq_files
    path parse_refgenome_files

    output:
    path("*")

    script:
    """
    echo "Mapping sample ${sample_id} with Parse Biosciences pipeline"

    # Generate parameter file
    echo "post_min_map_frac 0.01" > config_${params.seqTech}_parsepipe.txt
        
    split-pipe -m "all" \\
        --chemistry "v3" \\
        --input ${demux_fastq_files.join(' ')} \\
        --reference parse_refgenome_files \\
        --output_dir . \\
        --parfile config_${params.seqTech}_parsepipe.txt
    """
}

// ========================  SATURATION  ========================= \\ 
// Creates saturation plots using the tool: 10x_saturate           \\
process saturation {
    publishDir "${params.resDir}/saturation_plots", mode: 'symlink'
    debug true

    input:
    path mapping_files

    output:
    file("saturation_plot.png")

    script:
    """
    echo "Running saturation analysis for sample ${sample_id}"
    10x_saturate \\
        --input ${mapping_files.join(' ')} \\
        --output saturation_${sample_id}.png
    """
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