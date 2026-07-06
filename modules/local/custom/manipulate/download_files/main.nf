process DOWNLOAD_DATA {
    label 'process_single'

    output:
    path("bc_whitelist.txt"), emit: whitelist
    path "versions.yml",      emit: versions

    script:
    """
    # Download the whitelist file
    wget -O bc_whitelist.txt.gz ${params.bc_whitelist}

    # Unzip the whitelist file
    gunzip bc_whitelist.txt.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -n1 | sed 's/^GNU Wget //; s/ .*\$//')
    END_VERSIONS
    """
}
