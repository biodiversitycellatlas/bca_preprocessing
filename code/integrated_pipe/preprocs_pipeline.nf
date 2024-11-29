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
- incorporate seqspec
- split_parse_data: save unused reads to file and inspect

==========================
Help:

To print the output of processes:
println splitting_parse.out.splitted_files

==========================
Last run:
Exp: 240810_ParseBio_Nvec_Tcas
- Nvec: downloading -> starsolo index (66772b21-a2e5-45e8-acc0-577bc23aab14)
- Tcas: downloading -> mapping starsolo (prob) (45b71e5-9a04-443e-bcd2-ec69b0f8fa81) --- made mistake, go again?

Exp: 241106_BD_Rhapsody_Nvec
- Nvec: 

Currently running:
Exp: 240810_ParseBio_Nvec_Tcas
- Nvec: 
- Tcas: fixing mapping starsolo

Exp: 241106_BD_Rhapsody_Nvec
- Nvec: from start

------------------------------------------------------------------------------
*/

// ============== Define CUSTOM parameters ============== \\

// PARSE -- NVEC !!
// params.dataDir = "/users/asebe/bvanwaardenburg/git/bca_preprocessing/data"
// params.codeDir = "/users/asebe/bvanwaardenburg/git/bca_preprocessing/code/integrated_pipe"
// params.barcodeDir = "/users/asebe/bvanwaardenburg/ParseBiosciences-Pipeline.1.3.1/splitpipe/barcodes/"
// params.experiment = "240810_ParseBio_Nvec_Tcas"
// params.species = "Nvec"
// params.seqTech = "parse_biosciences"
// params.resDir = "${params.dataDir}/${params.experiment}/${params.species}_testIntegr"
// params.ref_star_gtf = "${params.resDir}/genome/Nvec_v5_merged_annotation_sort.gtf"
// params.ref_parse_gtf = "${params.resDir}/genome/Nvec_v4_merged_annotation_parse_sort.gtf"
// params.ref_fasta = "${params.resDir}/genome/Nvec_vc1.1_gDNA_mtDNA.fasta"

// PARSE -- TCAS !!
// params.dataDir = "/users/asebe/bvanwaardenburg/git/bca_preprocessing/data"
// params.codeDir = "/users/asebe/bvanwaardenburg/git/bca_preprocessing/code/integrated_pipe"
// params.barcodeDir = "/users/asebe/bvanwaardenburg/ParseBiosciences-Pipeline.1.3.1/splitpipe/barcodes/"
// params.experiment = "240810_ParseBio_Nvec_Tcas"
// params.species = "Tcas"
// params.seqTech = "parse_biosciences"
// params.resDir = "${params.dataDir}/${params.experiment}/${params.species}"
// params.ref_star_gtf = "${params.resDir}/genome/genomic.gtf"
// params.ref_parse_gtf = "${params.resDir}/genome/genomic.gtf"
// params.ref_fasta = "${params.resDir}/genome/GCF_031307605.1_icTriCast1.1_genomic.fna"

// BD RHAPSODY -- NVEC !!
params.baseDir = "/users/asebe/bvanwaardenburg/git/"
params.dataDir = "/users/asebe/bvanwaardenburg/git/bca_preprocessing/data"
params.codeDir = "/users/asebe/bvanwaardenburg/git/bca_preprocessing/code/integrated_pipe"
params.experiment = "241106_BD_Rhapsody_Nvec"
params.species = "Nvec"
params.seqTech = "bd_rhapsody"
params.barcodeDir = "/users/asebe/bvanwaardenburg/bca_preprocessing/seq_techniques/${params.seqTech}"
params.resDir = "${params.dataDir}/${params.experiment}"
params.ref_star_gtf = "${params.resDir}/genome/Nvec_v5_merged_annotation_sort.gtf"
params.ref_parse_gtf = "${params.resDir}/genome/Nvec_v4_merged_annotation_parse_sort.gtf"
params.ref_fasta = "${params.resDir}/genome/Nvec_vc1.1_gDNA_mtDNA.fasta"


// ============== Define BASE parameters ============== \\
params.star_config = "${params.codeDir}/config_${params.seqTech}_starsolo.txt"
params.saturationPath = "/users/asebe/bvanwaardenburg/git/10x_saturate"


// ====================  Channels ===================== \\
// Read initial sample IDs from accession list
Channel.fromPath("${params.dataDir}/accession_lists/${params.species}_${params.seqTech}_accessions.txt")
    .splitText()
    .map { it.trim() }
    .filter { it }
    .set { sample_ids }


// ====================  Workflow ===================== \\
// Define subworkflows
workflow parse_workflow {
    take:
        sample_ids
    main:
        download_data(sample_ids)
        splitting_parse(download_data.out)
        fastqc(splitting_parse.out.splitted_files)
        multiqc(fastqc.out.collect())
        
        genome_index_starsolo()
        ref_genome_parse()

        def mapping_input = splitting_parse.out.splitted_files
        .map { sample_id, fastq_files -> 
            [sample_id, fastq_files.join(' ')]
        }
        star_config = file(params.star_config)
        mapping_STARsolo(mapping_input, genome_index_starsolo.out, star_config)
        mapping_PB(mapping_input, ref_genome_parse.out)
    emit:
        mapping_PB.out
}

workflow bd_rhapsody_workflow {
    take:
        sample_ids
    main:
        download_data(sample_ids)
        fastqc(download_data.out)
        multiqc(fastqc.out.collect())
        genome_index_starsolo()
        
        def mapping_input = download_data.out
        .map { sample_id, fastq_files -> 
            [sample_id, fastq_files.join(' ')]
        }
        star_config = file(params.star_config)
        mapping_STARsolo(download_data.out, genome_index_starsolo.out, star_config)
    emit:
        mapping_STARsolo.out
}

// Main workflow
workflow {
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
    // saturation(mapping_output)
    // gene_ext(mapping_output)
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
process splitting_parse {
    publishDir "${params.resDir}/demux_fastq", mode: 'symlink'
    tag "${sample_id}"
    label 'big_cpus'
    debug true
    
    input:
    tuple val(sample_id), path(fastq_files)

    output:
    tuple val(sample_id), path("${sample_id}_*_R{1,2}.fastq.gz"), emit: splitted_files

    script:
    """
    echo "\n\n==================  splitting  =================="
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
    tag "${fastq_files}"
    debug true

    input:
    tuple val(sample_id), path(fastq_files)

    output:
    path("*_fastqc.*")

    script:
    """
    echo "\n\n==================  fastqc  =================="
    echo "Running FastQC for ${fastq_files}"
    echo "Path: ${fastq_files}"

    fastqc ${fastq_files} 
    """
}

// ==================  MULTI-QC  ================== \\ 
// A way of combining seperate FASTQC reports into  \\
// a single analysis. Provides an overview of the   \\
// data, including checks which of the files passed \\
// the quality metrics.                             \\
process multiqc {
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
    echo "Running MultiQC"
    multiqc .
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
// parse_refgenome: contains the reference genome for      \\
// the Parse Biosciences (split-pipe) pipeline.            \\
process ref_genome_parse {
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


// ==================  MAPPING STARSOLO  ================== \\ 
// Mapping using STARsolo
process mapping_STARsolo {    
    publishDir "${params.resDir}/mapping_STARsolo/${sample_id}", mode: 'symlink'
    debug true
    label 'big_mem'
    tag "${fastq_files}"

    input:
    tuple val(sample_id), val(fastq_files)
    path genome_index_files
    path star_config
    
    output:
    path("*")

    script:
    // def (r1_fastq, r2_fastq) = fastq_files.tokenize()
    """
    echo "\n\n==================  MAPPING STARSOLO  =================="
    echo "Mapping sample ${sample_id} with STARsolo"
    echo "Input files: ${fastq_files}"
    echo "Genome index directory: ${genome_index_files}"

    # Read configuration file
    config_file=\$(cat ${star_config})

    # Mapping step and generating count matrix using STAR
    STAR \\
        --genomeDir ${genome_index_files.toRealPath()} \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${sample_id}_ \\
        --readFilesIn (${fastq_files}) \\
        --soloCBwhitelist ${params.barcodeDir} \\
        --outSAMtype BAM SortedByCoordinate \\
        \${config_file} 
    """
}

// ==============  MAPPING PARSE PIPELINE  ============== \\ 
// Mapping using ParseBiosciences pipeline
process mapping_PB {
    tag "${fastq_files}"
    publishDir "${params.resDir}/mapping_parsepipe/${sample_id}", mode: 'symlink'
    debug true
    label 'big_mem'

    input:
    tuple val(sample_id), val(fastq_files)
    path parse_refgenome_files

    output:
    path("*")

    script:
    // def (r1_fastq, r2_fastq) = fastq_files.tokenize()
    """
    echo "\n\n==================  MAPPING Parse Biosciences  =================="
    echo "Mapping sample ${sample_id} with Parse Biosciences pipeline"
    echo "Input files: ${fastq_files}"
    echo "Genome index directory: ${parse_refgenome_files}"

    # Generate parameter file
    echo "post_min_map_frac 0.01" > config_${params.seqTech}_parsepipe.txt
        
    split-pipe -m "all" \\
        --chemistry "v3" \\
        --input (${fastq_files}) \\
        --reference ${parse_refgenome_files.toRealPath()} \\
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

    script:
    """
    echo "\n\n==================  SATURATION  =================="
    # Verify the input files
    echo "Processing bam file: ${sample_id}"

    # Run the script
    bash ${params.codeDir}/calculate_saturation.sh ${sample_id.toRealPath()} ${params.saturationPath} ${sample_id}
    """
}

// ========================  Gene Extension  ========================= \\ 
process gene_ext {
    publishDir "${params.resDir}/gene_ext", mode: 'symlink'
    conda '/users/asebe/bvanwaardenburg/miniconda3/envs/geneext'
    debug true

    input:
    path bam_file

    output:
    file("result.gtf")
    
    script:
    """
    echo "\n\n==================  Gene Extension  =================="
    cd ${params.baseDir}

    python geneext.py \\
        -g test_data/annotation.gtf \\
        -b test_data/alignments.bam \\
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