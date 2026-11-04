process KRAKEN {
    publishDir path: { "${params.outdir}/kraken/" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${meta.id}"

    // kraken2 loads the database's *.k2d files whole, so memory follows the database
    // actually used, downloaded or given in kraken_db_path, rather than the reads. Only
    // those files are measured: a database built in place also keeps its library/.
    // Coefficients live in params.dynamic_memory; without an entry the request is
    // process_high_memory2's 128 GB. 
    cpus   8
    time   { 45.m * task.attempt }
    memory { BcaResources.scaledMemory(
        params.dynamic_memory?.KRAKEN, BcaResources.filesEndingWith(kraken_db, '.k2d'), task.attempt, 128) }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/kraken2:2.17.1--9b129e69b7bcd776' :
        'community.wave.seqera.io/library/kraken2:2.17.1--beecc0e7e3db7c8b' }"

    input:
    path kraken_db
    tuple val(meta), path(filtered_fasta)

    output:
    path("${meta.id}.k2report"), emit: k2report
    path "versions.yml",         emit: versions

    script:
    """
    echo -e "\\n\\n==================  KRAKEN  =================="
    echo "Kraken database: ${kraken_db}"
    echo "Running KRAKEN for ${meta.id}"
    echo "FASTA file: ${filtered_fasta}"

    kraken2 \\
        --threads ${task.cpus} \\
        --db ${kraken_db} \\
        --report ${meta.id}.k2report \\
        --output ${meta.id}.kraken \\
        ${filtered_fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(kraken2 --version | head -n1 | sed 's/Kraken version //')
    END_VERSIONS
    """
}
