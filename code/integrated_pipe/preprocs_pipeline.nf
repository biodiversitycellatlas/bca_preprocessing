#!/usr/bin/env nextflow

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
// params.codeDir = "/users/asebe/bvanwaardenburg/git/bca_preprocessing/code/integrated_pipe"
// params.barcodeDir = "${params.baseDir}/bca_preprocessing/seq_techniques/${params.seqTech}/bc_data_n26_R1_v3_4 \
//                      ${params.baseDir}/bca_preprocessing/seq_techniques/${params.seqTech}/bc_data_v1 \
//                      ${params.baseDir}/bca_preprocessing/seq_techniques/${params.seqTech}/bc_data_R3_v3"
// params.resDir = "${params.dataDir}/${params.experiment}/${params.species}_testIntegr"
// params.ref_star_gtf = "${params.resDir}/genome/Nvec_v5_merged_annotation_sort.gtf"
// params.ref_parse_gtf = "${params.resDir}/genome/Nvec_v4_merged_annotation_parse_sort.gtf"
// params.ref_fasta = "${params.resDir}/genome/Nvec_vc1.1_gDNA_mtDNA.fasta"

// PARSE -- TCAS !!
params.experiment = "240810_ParseBio_Nvec_Tcas"
params.species = "Tcas"
params.seqTech = "parse_biosciences"
params.baseDir = "/users/asebe/bvanwaardenburg/git/"
params.dataDir = "/users/asebe/bvanwaardenburg/git/bca_preprocessing/data"
params.codeDir = "/users/asebe/bvanwaardenburg/git/bca_preprocessing/code/integrated_pipe"
params.barcodeDemux = "${params.baseDir}/bca_preprocessing/seq_techniques/${params.seqTech}/bc_data_n26_R1_v3_4.csv"
params.barcodeDir = "${params.baseDir}/bca_preprocessing/seq_techniques/${params.seqTech}/bc_data_n26_R1_v3_4 \
                     ${params.baseDir}/bca_preprocessing/seq_techniques/${params.seqTech}/bc_data_v1 \
                     ${params.baseDir}/bca_preprocessing/seq_techniques/${params.seqTech}/bc_data_R3_v3"
params.resDir = "${params.dataDir}/${params.experiment}/${params.species}"
params.ref_star_gtf = "${params.resDir}/genome/genomic.gtf"
params.ref_parse_gtf = "${params.resDir}/genome/genomic.gtf"
params.ref_fasta = "${params.resDir}/genome/GCF_031307605.1_icTriCast1.1_genomic.fna"

// BD RHAPSODY -- NVEC !!
// params.experiment = "241106_BD_Rhapsody_Nvec"
// params.species = "Nvec"
// params.seqTech = "bd_rhapsody"
// params.baseDir = "/users/asebe/bvanwaardenburg/git/"
// params.dataDir = "/users/asebe/bvanwaardenburg/git/bca_preprocessing/data"
// params.codeDir = "/users/asebe/bvanwaardenburg/git/bca_preprocessing/code/integrated_pipe"
// params.barcodeDir = "${params.baseDir}/bca_preprocessing/seq_techniques/${params.seqTech}/BD_CLS1.txt \
//                         ${params.baseDir}/bca_preprocessing/seq_techniques/${params.seqTech}/BD_CLS2.txt \
//                         ${params.baseDir}/bca_preprocessing/seq_techniques/${params.seqTech}/BD_CLS3.txt"
// params.resDir = "${params.dataDir}/${params.experiment}"
// params.ref_star_gtf = "${params.resDir}/genome/Nvec_v5_merged_annotation_sort.gtf"
// params.ref_parse_gtf = "${params.resDir}/genome/Nvec_v4_merged_annotation_parse_sort.gtf"
// params.ref_fasta = "${params.resDir}/genome/Nvec_vc1.1_gDNA_mtDNA.fasta"


// ============== DEFINE BASE PARAMETERS ============== \\
params.star_config = "${params.codeDir}/config_${params.seqTech}_starsolo.txt"
params.saturationPath = "/users/asebe/bvanwaardenburg/git/10x_saturate"


// ====================  CHANNELS ===================== \\
// Read initial sample IDs from accession list
Channel.fromPath("${params.dataDir}/accession_lists/${params.species}_${params.seqTech}_accessions.txt")
    .splitText()
    .map { it.trim() }
    .filter { it }
    .set { sample_ids }


// ====================  WORKFLOWS ===================== \\

// =================  SUBWORKFLOWS  ================= \\ 
// The subworkflows are called from the main workflow \\
// (see below) and are created to handle different    \\
// sequencing techniques.                             \\
// Available:                                         \\
// - parse_workflow: Parse Biosciences                \\
// - bd_rhapsody_workflow: BD Rhapsody                \\

// ===============  Parse Biosciences  ==============  \\ 
// Based on the sample wells, the fastq sequences      \\
// must be split based on the first barcode round.     \\
// Performs both custom- and commercial pre-processing \\
// (using the Parse Biosciences pipeline v1.3.1)       \\
// enabling comparison between the methods and         \\
// validation of steps.       
workflow parse_workflow {
    take:
        sample_ids
    main:
        DOWNLOAD_DATA(sample_ids)
        DEMULTIPLEX(DOWNLOAD_DATA.out)
        FASTQC(DEMULTIPLEX.out.splitted_files)
        MULTIQC(FASTQC.out.collect())
        
        GENINDEX_STARSOLO()
        REFGEN_PARSEBIO()

        // def mapping_input = DEMULTIPLEX.out.splitted_files
        // .map { sample_id, fastq_files -> 
        //     [sample_id, fastq_files.join(' ')]
        // }
        star_config = file(params.star_config)
        MAPPING_STARSOLO(DEMULTIPLEX.out.splitted_files, GENINDEX_STARSOLO.out, star_config)
        MAPPING_PARSEBIO(DEMULTIPLEX.out.splitted_files, REFGEN_PARSEBIO.out)
        INDEX_BAM(MAPPING_STARSOLO.out)
    emit:
        MAPPING_STARSOLO.out
}

// ===================  BD Rhapsody ================  \\ 
// A basic approach, where first quality control is   \\
// performed and mapping using the (complete) fastq   \\
// sequences.                                         \\
workflow bd_rhapsody_workflow {
    take:
        sample_ids
    main:
        DOWNLOAD_DATA(sample_ids)
        FASTQC(DOWNLOAD_DATA.out)
        MULTIQC(FASTQC.out.collect())
        GENINDEX_STARSOLO()
        
        // def mapping_input = DOWNLOAD_DATA.out
        // .map { sample_id, fastq_files -> 
        //     [sample_id, fastq_files.join(' ')]
        // }
        star_config = file(params.star_config)
        MAPPING_STARSOLO(DOWNLOAD_DATA.out, GENINDEX_STARSOLO.out, star_config)
        INDEX_BAM(MAPPING_STARSOLO.out)
    emit:
        MAPPING_STARSOLO.out
}

// =================  MAIN WORKFLOW  ================= \\ 
// Selects pre-processing workflow depending on the    \\
// sequencing technique and returns a bam file         \\
// (output alignment step). This will be used in       \\
// downstream processes, which are identical for all   \\
// techniques.                                         \\
workflow {
    // 
    if (params.seqTech == 'parse_biosciences') {
        parse_workflow(sample_ids)
        mapping_output = parse_workflow.out
    } else if (params.seqTech == 'bd_rhapsody') {
        bd_rhapsody_workflow(sample_ids)
        mapping_output = bd_rhapsody_workflow.out
    } else {
        error "Invalid sequencing technology specified. Use 'parse' or 'starsolo'."
    }

    // Downstream processes
    SATURATION(mapping_output)
    GENE_EXT(mapping_output)
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
process DOWNLOAD_DATA {
    tag "${sample_id}"
    debug true
    
    input:
    val sample_id

    output:
    tuple val(sample_id), path("${sample_id}*R{1,2}*.fastq.gz")

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
process DEMULTIPLEX {
    publishDir "${params.resDir}/demux_fastq", mode: 'symlink'
    tag "${sample_id}"
    label 'big_mem'
    debug true
    
    input:
    tuple val(sample_id), path(fastq_files)

    output:
    tuple val(sample_id), path("*.fastq.gz"), emit: splitted_files

    script:
    def fastq_list = fastq_files instanceof List ? fastq_files : [fastq_files]
    def r1_fastq = fastq_list.find { it.name.contains('_R1') }
    def r2_fastq = fastq_list.find { it.name.contains('_R2') }
    """
    echo "\n\n==================  splitting  =================="
    echo "Processing sample: ${sample_id}"
    echo "First barcode path: ${params.barcodeDemux}"
    echo "FQ 1: ${r1_fastq ?: 'Not provided'}"
    echo "FQ 2: ${r2_fastq ?: 'Not provided'}"

    # Run the demux code
    bash ${params.codeDir}/split_parse_data.sh ${params.resDir} ${sample_id} ${r1_fastq} ${r2_fastq} ${params.barcodeDemux}
    """
}

// ==================  FASTQC  ================== \\ 
// Generate Quality Control reports using FASTQC  \\
process FASTQC {
    publishDir "${params.resDir}/fastqc", mode: 'symlink'
    tag "${fastq_files}"
    debug true

    input:
    tuple val(sample_id), path(fastq_files)

    output:
    path("*_fastqc.*")

    script:
    """
    echo "\n\n==================  FASTQC  =================="
    echo "Running FASTQC for ${fastq_files}"
    echo "Path: ${fastq_files}"

    fastqc ${fastq_files} 
    """
}

// ==================  MULTI-QC  ================== \\ 
// A way of combining seperate FASTQC reports into  \\
// a single analysis. Provides an overview of the   \\
// data, including checks which of the files passed \\
// the quality metrics.                             \\
process MULTIQC {
    publishDir "${params.resDir}/fastqc", mode: 'symlink'
    tag "all"
    debug true
    
    input:
    path('*')

    output:
    file("multiqc_report.html")

    script:
    """
    echo "\n\n==================  Multi qc  =================="
    echo "Running MULTIQC"
    multiqc .
    """
}

// =============  GENOME INDEX STARSOLO  ============= \\ 
// Creates the Genome Indeces required for mapping.    \\
// Two folders will be created within /genome/,        \\
// genome_index: contains genome index for STARsolo    \\
// parse_refgenome: see next process                   \\
process GENINDEX_STARSOLO {   
    publishDir "${params.resDir}/genome/genome_index", mode: 'symlink'
    debug true

    output:
    path("*")     
       
    script:
    """
    echo "\n\n==================  GENOME INDEX STARSOLO  =================="
    # Retrieve the first accession number
    first_accs=\$(head -1 ${params.dataDir}/accession_lists/${params.species}_${params.seqTech}_accessions.txt)
    first_fastq="${params.resDir}/fastq/\${first_accs}*1*.fastq.gz"  

    echo "\${first_accs}"
    echo "\${first_fastq}"

    # Calculate SJDB overhang using the first read from the first fastq file
    sjdb_overhang=\$(zcat \${first_fastq} 2>/dev/null | awk 'NR==2 {print length(\$0)-1; exit}' || echo "") 

    echo "\${sjdb_overhang}"

    echo "Generating genome index with STAR"
    STAR --runMode genomeGenerate \\
        --genomeFastaFiles ${params.ref_fasta} \\
        --sjdbGTFfile ${params.ref_star_gtf} \\
        --sjdbOverhang "\${sjdb_overhang}" \\
        --genomeSAindexNbases 12
    """
}

// =============  REF GENOME PARSE PIPELINE  ============= \\ 
// Creates the Genome Indeces required for mapping.        \\
// parse_refgenome: contains the reference genome for the  \\
// Parse Biosciences (split-pipe) pipeline.                \\
process REFGEN_PARSEBIO {
    publishDir "${params.resDir}/genome/parse_refgenome", mode: 'symlink'
    debug true

    output:
    path("*")

    script:
    """
    echo "\n\n==================  REF GENOME PARSE PIPELINE  =================="
    split-pipe -m mkref \\
        --genome_name ${params.species} \\
        --genes ${params.ref_parse_gtf} \\
        --fasta ${params.ref_fasta} \\
        --output_dir .
    """
}


// ==============  MAPPING STARSOLO  =============== \\ 
// Mapping using STAR (version 2.7.11b).             \\
// The functions specify read 1 and read 2 from      \\
// the channel, and otherwise state if not provided. \\
// The settings for STAR are set in the variable     \\
// config_file, and is specific per sequencing tech. \\
process MAPPING_STARSOLO {    
    publishDir "${params.resDir}/mapping_STARsolo/${sample_id}", mode: 'symlink'
    debug true
    label 'big_mem'
    tag "${fastq_files}"

    input:
    tuple val(sample_id), path(fastq_files)
    path genome_index_files
    path star_config
    
    output:
    path("*")

    script:
    def fastq_list = fastq_files instanceof List ? fastq_files : [fastq_files]
    def r1_fastq = fastq_list.find { it.name.contains('_R1') }
    def r2_fastq = fastq_list.find { it.name.contains('_R2') }
    """
    echo "\n\n==============  MAPPING STARSOLO  ================"
    echo "Mapping sample ${sample_id} with STARsolo"
    echo "FQ 1: ${r1_fastq ?: 'Not provided'}"
    echo "FQ 2: ${r2_fastq ?: 'Not provided'}"
    echo "Genome index directory: ${genome_index_files}"

    # Read configuration file
    config_file=\$(cat ${star_config})

    # Mapping step and generating count matrix using STAR
    STAR \\
        --genomeDir ${genome_index_files.toRealPath()} \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${sample_id}_ \\
        --readFilesIn ${r1_fastq} ${r2_fastq} \\
        --soloCBwhitelist ${params.barcodeDir} \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMattributes CR UR CB UB \\
        \${config_file} 
    """
}

// ============  MAPPING PARSE PIPELINE  =========== \\ 
// Mapping using ParseBiosciences pipeline 1.3.1     \\
// The functions specify read 1 and read 2 from      \\
// the channel, and otherwise state if not provided. \\
// The settings are set in the variable config_file, \\
// and is specific per sequencing tech.              \\
// As all previously created reference files are     \\
// given as one string, these are added to a new     \\
// directory.                                        \\
process MAPPING_PARSEBIO {
    tag "${fastq_files}"
    publishDir "${params.resDir}/mapping_parsepipe/${sample_id}", mode: 'symlink'
    debug true
    label 'big_mem'

    input:
    tuple val(sample_id), path(fastq_files)
    path parse_refgenome_files

    output:
    path("*")

    script:
    def fastq_list = fastq_files instanceof List ? fastq_files : [fastq_files]
    def r1_fastq = fastq_list.find { it.name.contains('_R1') }
    def r2_fastq = fastq_list.find { it.name.contains('_R2') }
    """
    echo "\n\n=============  MAPPING PARSE BIOSCIENCES  ================"
    echo "Mapping sample ${sample_id} with Parse Biosciences pipeline"
    echo "FQ 1: ${r1_fastq ?: 'Not provided'}"
    echo "FQ 2: ${r2_fastq ?: 'Not provided'}"
    echo "Genome index directory: ${parse_refgenome_files}"

    # Generate parameter file
    echo "post_min_map_frac 0.01" > config_${params.seqTech}_parsepipe.txt

    # Create directory for the genome index files
    mkdir -p genome_index
    
    # Move all genome index files to the new directory
    mv ${parse_refgenome_files} genome_index/
        
    split-pipe -m all \\
        --chemistry v3 \\
        --fq1 ${r1_fastq} \\
        --fq2 ${r2_fastq} \\
        --genome_dir genome_index \\
        --output_dir . \\
        --parfile config_${params.seqTech}_parsepipe.txt
    """
}

// =================  INDEX BAM FILES  ================== \\ 
process INDEX_BAM {
    publishDir "${params.resDir}/mapping_parsepipe/${sample_id}", mode: 'symlink'
    debug true

    input:
    path(mapping_files)

    output:
    file("Aligned.sortedByCoord.out.bam.bai")

    script:
    """
    echo "\n\n==================  INDEX BAM FILES  =================="
    echo "Processing files: ${mapping_files}"

    samtools index "Aligned.sortedByCoord.out.bam"
    """
}

// ====================  SATURATION  ===================== \\ 
// Creates saturation plots using the tool: 10x_saturate   \\
// which is an external package linked using the github    \\
// submodule function.                                     \\
process SATURATION {
    publishDir "${params.resDir}/saturation_plots", mode: 'symlink'
    debug true

    input:
    path(mapping_files)

    script:
    """
    echo "\n\n==================  SATURATION  =================="
    # Verify the input files
    echo "Processing bam file: ${mapping_files}"

    # Run the script
    bash ${params.codeDir}/calculate_saturation.sh ${mapping_files.toRealPath()} ${params.saturationPath} ${mapping_files}
    """
}

// ==================  GENE EXTENSION  =================== \\ 
process GENE_EXT {
    publishDir "${params.resDir}/gene_ext", mode: 'symlink'
    conda '/users/asebe/bvanwaardenburg/miniconda3/envs/geneext'
    debug true

    input:
    path(mapping_files)

    output:
    file("result.gtf")
    
    script:
    """
    echo "\n\n==================  GENE EXTENSION  =================="
    echo "BAM filepath: ${mapping_files}"
    echo "Original GTF: ${params.ref_star_gtf}"
    
    python ${params.baseDir}/ext_programs/GeneExt/geneext.py \\
        -g ${params.ref_star_gtf} \\
        -b ${mapping_files}/*.bam \\
        -o result.gtf \\
        --peak_perc 0
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
