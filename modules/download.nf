process DOWNLOAD_DATA {
    output:
    path("bc_whitelist.txt")

    script:
    bc_whitelist  = params.seqtech_parameters[params.protocol].bc_whitelist

    """
    # Download the whitelist file
    wget -O bc_whitelist.txt.gz ${bc_whitelist}
    
    # Unzip the whitelist file
    gunzip bc_whitelist.txt.gz
    """
}
