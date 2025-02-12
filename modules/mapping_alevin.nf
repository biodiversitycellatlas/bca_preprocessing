process MAPPING_ALEVIN {   
    input:
    tuple val(sample_id), path(fastq_files)
    path index
    path transcript_tsv
    path whitelist

    output:
    tuple val(meta), path("*_alevin_results"), emit: alevin_results
    path  "versions.yml"                     , emit: versions
       
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
    }
    """
    echo "\n\n==================  GENOME INDEX ALEVIN =================="
    echo "Reference fasta: ${params.ref_fasta}"

    # Define filename pattern
    reads1_pat="_R1_"
    reads2_pat="_R2_"

    # Obtain and sort filenames
    reads1="\$(find -L ${FASTQ_DIR} -name "*\$reads1_pat*" -type f | sort | awk -v OFS=, '{\$1=\$1;print}' | paste -sd, -)"
    reads2="\$(find -L ${FASTQ_DIR} -name "*\$reads2_pat*" -type f | sort | awk -v OFS=, '{\$1=\$1;print}' | paste -sd, -)"

    # Export required var
    export ALEVIN_FRY_HOME=.
    export NUMBA_CACHE_DIR=.

    # simpleaf configuration
    simpleaf set-paths

    # Quantify the sample
    simpleaf quant \\
        --reads1 \${reads1} \\
        --reads2 \${reads2} \\
        --index ${index} \\
        --t2g-map ${transcript_tsv} \\
        --chemistry "$protocol" 
    """
}
