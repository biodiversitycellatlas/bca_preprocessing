process MAPPING_ALEVIN {
    publishDir "${params.resDir}/mapping_alevin/${sample_id}", mode: 'copy'
    tag "${sample_id}"
    debug true

    input:
    tuple val(sample_id), path(fastq_files)
    path index
    path transcript_tsv

    output:
    tuple val(meta), path("*_alevin_results"), emit: alevin_results
       
    script:
    def fastq_list = fastq_files instanceof List ? fastq_files : [fastq_files]
    def r1_fastq = fastq_list.find { it.name.contains('_R1') }
    def r2_fastq = fastq_list.find { it.name.contains('_R2') }

    // If seqTech is "bd_rhapsody", then cDNA = R2 and CB/UMI = R1
    // Else by default cDNA = R1 and CB/UMI = R2
    def cDNA_read
    def CBUMI_read
    if (params.seqTech.toLowerCase().contains("bd_rhapsody")) {
        cDNA_read = r2_fastq
        CBUMI_read = r1_fastq
    } else {
        cDNA_read = r1_fastq
        CBUMI_read = r2_fastq
        protocol = "splitseqV2"
    }
    """
    echo "\n\n==================  GENOME INDEX ALEVIN =================="
    echo "Sample ID: ${sample_id}"
    echo "Transcript TSV: ${transcript_tsv}"
    echo "Reference fasta: ${params.ref_fasta}"

    # Export required vars
    export ALEVIN_FRY_HOME=.

    # simpleaf configuration
    simpleaf set-paths

    # Quantify the sample
    simpleaf quant \\
        --reads1 ${cDNA_read} \\
        --reads2 ${CBUMI_read} \\
        --index ${index} \\
        --t2g-map ${transcript_tsv} \\
        --chemistry $protocol
    """
}
