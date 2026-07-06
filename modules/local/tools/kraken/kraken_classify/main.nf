process KRAKEN {
    publishDir "${params.outdir}/kraken/", mode: 'copy'
    tag "${meta.id}"
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/kraken2:2.17.1--9b129e69b7bcd776"

    input:
    path db_path_file
    tuple val(meta), path(filtered_fasta)

    output:
    path("${meta.id}.k2report"), emit: k2report
    path "versions.yml",         emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(kraken2 --version | head -n1 | sed 's/Kraken version //')
    END_VERSIONS
    """
}
