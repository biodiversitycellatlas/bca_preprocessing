process KRAKEN {
    publishDir "${params.outdir}/kraken/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    label 'process_high'
    debug true

    input:
    path db_path_file
    tuple val(meta), path(filtered_fasta)

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/kraken2_samtools_coreutils_pigz:8961943c277652a6"

    output:
    path("*")

    script:
    """
    echo "\n\n==================  KRAKEN  =================="
    echo "Kraken db path file: ${db_path_file}"
    echo "Running KRAKEN for ${meta.id}"
    echo "FASTA file: ${filtered_fasta}"

    kraken_db_path=\$(cat ${db_path_file})

    k2 classify \\
        --threads ${task.cpus} \\
        --db \${kraken_db_path} \\
        --report ${meta.id}_kraken_taxonomy.txt \\
        --report-minimizer-data \\
        --use-names \\
        --memory-mapping \
        --log ${meta.id}_kraken.log \\
        --output ${meta.id}_kraken_output.txt \\
        ${filtered_fasta}
    """
}
