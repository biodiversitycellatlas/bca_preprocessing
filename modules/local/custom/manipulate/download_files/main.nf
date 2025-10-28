process DOWNLOAD_DATA {
    label 'process_single'

    output:
    path("bc_whitelist.txt")

    script:
    """
    # Download the whitelist file
    wget -O bc_whitelist.txt.gz ${params.bc_whitelist}

    # Unzip the whitelist file
    gunzip bc_whitelist.txt.gz
    """
}
