process FASTQ_SPLITTER {
    publishDir "${params.output_dir}/fastqsplitter/${sample_id}", mode: 'copy'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(fastq_files)

    output:
    tuple val(sample_id), file("raw_reads_split/${sample_id}_R1_*"), file("raw_reads_split/${sample_id}_R2_*") into split_files_ch

    script:
    def outArgs_R1 = (1..params.fastq_chunks).collect { i -> 
        "-o raw_reads_split/${sample_id}_R1_${i}-of-${params.fastq_chunks}.fastq.gz" 
    }.join(" ")
    def outArgs_R2 = (1..params.fastq_chunks).collect { i -> 
        "-o raw_reads_split/${sample_id}_R2_${i}-of-${params.fastq_chunks}.fastq.gz" 
    }.join(" ")

    """
    echo "\n\n==================  TRIM FASTQs WITH FASTP  =================="
    echo "Sample ID: ${sample_id}"
    echo "Processing files: ${fastq_files}"

    mkdir -p raw_reads_split

    # Split R1 file
    fastqsplitter -i ${fastq_files[0]} ${outArgs_R1} -t 1 -c 1
    
    # Split R2 file
    fastqsplitter -i ${fastq_files[1]} ${outArgs_R2} -t 1 -c 1
    """
}
