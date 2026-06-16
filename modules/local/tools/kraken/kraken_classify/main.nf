process KRAKEN {
    publishDir "${params.outdir}/kraken/", mode: 'copy'
    tag "${meta.id}"
    label 'process_high_memory'
    debug true

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/kraken2_samtools_coreutils_pigz:8961943c277652a6"

    input:
    path db_path_file
    tuple val(meta), path(filtered_fasta)

    output:
    path("${meta.id}.k2report")

    script:
    """
    echo "\n\n==================  KRAKEN  =================="
    echo "Kraken db path file: ${db_path_file}"
    echo "Running KRAKEN for ${meta.id}"
    echo "FASTA file: ${filtered_fasta}"

    kraken_db_path=\$(cat ${db_path_file})

    kraken2 \\
        --threads ${task.cpus} \\
        --db \${kraken_db_path} \\
        --report ${meta.id}.k2report \\
        --output ${meta.id}.kraken \\
        ${filtered_fasta}
    """
}
