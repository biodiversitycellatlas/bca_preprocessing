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
- assign demux_sample_id

- FIX: Output and Input Files Outside Working Directory:
    Nextflow expects processes to produce outputs in their working 
    directories ($PWD during execution). By specifying output files 
    located outside the working directory (e.g., ${params.resDir}/fastq/), 
    Nextflow may not be able to track these files properly, leading to 
    issues with process execution, reproducibility, and caching

------------------------------------------------------------------------------
*/

// Define parameters
params.dataDir = "/users/asebe/bvanwaardenburg/git/bca_preprocessing/data"
params.codeDir = "/users/asebe/bvanwaardenburg/git/bca_preprocessing/code/integerated_pipe"
params.barcodeDir = "/users/asebe/bvanwaardenburg/ParseBiosciences-Pipeline.1.3.1/splitpipe/barcodes/"
params.species = 'nvec'
params.seqTech = "parse"
params.resDir = "${params.dataDir}/${params.species}"

// Define channels
// Read initial sample IDs from accession list
Channel.fromPath("${params.dataDir}/accession_lists/${params.species}_accessions.txt")
    .splitText()
    .map { it.trim() }
    .filter { it }
    .set { sample_ids }

// Workflow definition
workflow {
    download_data(sample_ids)                           // Download data
    splitting_parse(download_data.out)                  // Demultiplexing
    fastqc(splitting_parse.out)                         // Quality control with FastQC
    multiqc(fastqc.out)                                 // Generate MultiQC report
    genome_index_starsolo()                             // Generate genome index (only once)
    ref_genome_parse()                                  // Generate reference genome (only once)
    mapping_STARsolo(splitting_parse.out, genome_index_starsolo.out)      // Mapping with STARsolo
    mapping_PB(splitting_parse.out, ref_genome_parse.out) // Mapping with Parse Biosciences pipeline
    saturation(mapping_PB.out)                          // Saturation analysis
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
    
    input:
    val sample_id

    output:
    set val(sample_id), file("${sample_id}_*.fastq.gz") into fastq_files

    publishDir "${params.resDir}/fastq", mode: 'move'  // Use 'move' instead of 'symlink'

    script:
    """
    # Check if FASTQ files already exist
    if [ -f "${params.resDir}/fastq/${sample_id}_1.fastq.gz" ]; then
        echo "FASTQ files for sample ${sample_id} already exist."
        ln -s "${params.resDir}/fastq/${sample_id}_1.fastq.gz" .
        ln -s "${params.resDir}/fastq/${sample_id}_2.fastq.gz" .
    else
        echo "Downloading data for sample ${sample_id}"
        prefetch ${sample_id}
        fastq-dump --split-files --gzip --outdir . ${sample_id}
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
process splitting_parse {
    tag "${sample_id}"
    
    input:
    val sample_id
    set val(sample_id), file(fastq_files) from fastq_files.filter{ it[0]==sample_id }

    output:
    set val(demux_sample_id), file("*.fastq.gz") into demuxed_fastqs

    publishDir "${params.resDir}/demux_fastq", mode: 'move'

    script:
    """
    # Run the demux code without copying files
    bash ${params.codeDir}/split_parse_data.sh ${sample_id} \\
        ${params.resDir}/fastq \\
        ${params.barcodeDir} \\
        .

    # Extract demux_sample_id from the output files
    demux_sample_id=\$(ls *.fastq.gz | cut -d'_' -f1 | sort | uniq)
    echo "\${demux_sample_id}" > demux_sample_ids.txt
    """
}

// ==================  FASTQC  ================== \\ 
// Generate Quality Control reports using FASTQC  \\
process fastqc {
    tag "${demux_sample_id}"
    
    input:
    set val(demux_sample_id), file(demux_fastq_files) from demuxed_fastqs.groupTuple()

    output:
    set val(demux_sample_id), file("*_fastqc.*") into fastqc_reports

    publishDir "${params.resDir}/fastqc", mode: 'move'

    script:
    """
    echo "Running FastQC for sample ${demux_sample_id}"
    fastqc ${demux_fastq_files.join(' ')} --outdir .
    """
}

// ==================  MULTI-QC  ================== \\ 
// A way of combining seperate FASTQC reports into  \\
// a single analysis. Provides an overview of the   \\
// data, including checks which of the files passed \\
// the quality metrics.                             \\
process multiqc {
    input:
    file fastqc_reports from fastqc_reports.collect()

    output:
    file("multiqc_report.html") into multiqc_report

    publishDir "${params.resDir}/fastqc", mode: 'move'
    
    script:
    """
    echo "Running MultiQC"
    multiqc . --filename multiqc_report.html
    """
}

// =============  GENOME INDEX STARSOLO  ============= \\ 
// Creates the Genome Indeces required for mapping.    \\
// Two folders will be created within /genome/,        \\
// genome_index: contains genome index for STARsolo    \\
// parse_refgenome: see next process                   \\
process genome_index_starsolo {
    publishDir "${params.resDir}/genome/genome_index", mode: 'move'
        
    output:
    file("genome_index/*") into genome_index_files
       
    script:
    """
    # Retrieve the first accession number
    first_accs=\$(head -1 ${params.dataDir}/accession_lists/${params.species}_accessions.txt)
    first_fastq="${params.resDir}/fastq/\${first_accs}_1.fastq.gz"  

    # Calculate SJDB overhang using the first read from the first fastq file
    sjdb_overhang=\$(zcat "\${first_fastq}" 2>/dev/null | awk 'NR==2 {print length(\$0)-1; exit}' || echo "") 

    echo "Generating genome index with STAR"
    mkdir -p genome_index
    STAR --runMode genomeGenerate \\
        --genomeDir genome_index \\
        --genomeFastaFiles ${params.resDir}/genome/Nvec_vc1.1_gDNA_mtDNA.fasta \\
        --sjdbGTFfile ${params.resDir}/genome/Nvec_v5_merged_annotation_sort.gtf \\
        --sjdbOverhang \${sjdb_overhang} \\
        --genomeSAindexNbases 12
    """
}

// =============  REF GENOME PARSE PIPELINE  ============= \\ 
// Creates the Genome Indeces required for mapping.        \\
// parse_refgenome: contains the reference genome for      \\
// the Parse Biosciences (split-pipe) pipeline.            \\
process ref_genome_parse {
    publishDir "${params.resDir}/genome/parse_refgenome", mode: 'move'
        
    output:
    file("parse_refgenome/*") into parse_refgenome_files
       
    script:
    """
    ref_gtf="${params.resDir}/genome/Nvec_v4_merged_annotation_parse_sort.gtf"
    ref_fasta="${params.resDir}/genome/Nvec_vc1.1_gDNA_mtDNA.fasta"    

    mkdir -p parse_refgenome

    split-pipe -m mkref \\
        --genome_name ${params.species} \\
        --genes \${ref_gtf} \\
        --fasta \${ref_fasta} \\
        --output_dir parse_refgenome
    """
}

// ==================  MAPPING STARSOLO  ================== \\ 
// Mapping using STARsolo
process mapping_STARsolo {
    tag "${demux_sample_id}"
    
    input:
    val demux_sample_id
    set val(demux_sample_id), file(demux_fastq_files) from demuxed_fastqs.groupTuple()
    file genome_index_files from genome_index_files.collect()
    file config_file from "${params.codeDir}/config_${params.seqTech}_starsolo.txt"

    output:
    set val(demux_sample_id), file("${demux_sample_id}*") into mapping_results

    publishDir "${params.resDir}/mapping_STARsolo/${demux_sample_id}", mode: 'move'

    script:
    """
    echo "Mapping sample ${demux_sample_id} with STARsolo"

    # Read configuration file
    CONFIG_OPTIONS=\$(cat ${config_file})

    # Mapping step and generating count matrix using STAR
    STAR \\
        --genomeDir genome_index_files \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${demux_sample_id}_ \\
        --readFilesIn ${demux_fastq_files.join(' ')} \\
        --soloCBwhitelist ${params.barcodeDir} \\
        \${CONFIG_OPTIONS}

    # Create the index file
    samtools index ${demux_sample_id}_Aligned.sortedByCoord.out.bam
    """
}

// ==============  MAPPING PARSE PIPELINE  ============== \\ 
// Mapping using ParseBiosciences pipeline
process mapping_PB {
    tag "${demux_sample_id}"
    
    input:
    val demux_sample_id
    set val(demux_sample_id), file(demux_fastq_files) from demuxed_fastqs.groupTuple()
    file parse_refgenome_files from parse_refgenome_files.collect()

    output:
    set val(demux_sample_id), file("*") into mapping_PB_results

    publishDir "${params.resDir}/mapping_parsepipe/${demux_sample_id}", mode: 'move'

    script:
    """
    echo "Mapping sample ${demux_sample_id} with Parse Biosciences pipeline"

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
    tag "${demux_sample_id}"
    
    input:
    val demux_sample_id
    file(mapping_files) from mapping_PB_results.filter{ it[0]==demux_sample_id }.collect()

    output:
    set val(demux_sample_id), file("saturation_${demux_sample_id}.png") into saturation_plots

    publishDir "${params.resDir}/saturation_plots", mode: 'move'

    script:
    """
    echo "Running saturation analysis for sample ${demux_sample_id}"
    10x_saturate \\
        --input ${mapping_files.join(' ')} \\
        --output saturation_${demux_sample_id}.png
    """
}

