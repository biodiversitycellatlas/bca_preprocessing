// ==============  MAPPING STARSOLO  =============== \\ 
// Mapping using STAR (version 2.7.11b).             \\
// The functions specify read 1 and read 2 from      \\
// the channel, and otherwise state if not provided. \\
// The settings for STAR are set in the variable     \\
// config_file, and is specific per sequencing tech. \\

process MAPPING_STARSOLO_SINGLE { 
    publishDir "${params.resDir}/mapping_STARsolo/mapping_STARsolo_single/${sample_id}", mode: 'copy'
    tag "${sample_id}_STARsolo_single"

    input:
    tuple val(sample_id), path(fastq_files)
    path genome_index_files
    
    output:
    tuple val(sample_id), path("*")

    script:
    def fastq_list = fastq_files instanceof List ? fastq_files : [fastq_files]
    def r1_fastq = fastq_list.find { it.name.contains('_R1') }
    def r2_fastq = fastq_list.find { it.name.contains('_R2') }

    // If seqTech is "bd_rhapsody", then cDNA = R2 and CB/UMI = R1, else by default cDNA = R1 and CB/UMI = R2
    def cDNA_read
    if (params.seqTech.toLowerCase().contains("bd_rhapsody")) {
        cDNA_read = r2_fastq
    } else {
        cDNA_read = r1_fastq
    }

    """
    echo "\n\n==============  MAPPING STARSOLO SINGLE ================"
    echo "Mapping sample ${sample_id} with STARsolo"
    echo "cDNA read: ${cDNA_read ?: 'Not provided'}"
    echo "Genome index directory: ${genome_index_files}"

    # Read configuration file
    config_file=\$(cat ${params.star_config_single})

    # Mapping step and generating count matrix using STAR
    STAR \\
        --runThreadN 4 \\
        --readFilesIn ${cDNA_read} \\
        --genomeDir ${genome_index_files.toRealPath()} \\
        --readFilesCommand zcat \\
        --outSAMtype BAM SortedByCoordinate \\
        --soloFeatures Gene GeneFull \\
        --outFileNamePrefix ${sample_id}_ \\
        \${config_file} 

    # Assign reads to genes
    featureCounts -a ${params.ref_star_gtf} -o gene_assigned -R BAM *_Aligned.sortedByCoord.out.bam -T 4

    # Sort the output
    samtools sort *_Aligned.sortedByCoord.out.bam.featureCounts.bam -o assigned_sorted.bam
    samtools index assigned_sorted.bam

    # Counting molecules
    umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I assigned_sorted.bam -S counts.tsv.gz
    """
}
