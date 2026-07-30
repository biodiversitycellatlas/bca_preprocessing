process DOWNLOAD_WHITELIST {
    label 'process_single'

    input:
    val(whitelist_urls)

    output:
    path("bc_whitelist_*.txt"), emit: whitelist
    path "versions.yml",        emit: versions

    script:
    """
    # Support one or more whitespace-separated whitelist URLs
    urls=( ${whitelist_urls} )

    for i in "\${!urls[@]}"; do
        # Zero-padded index keeps the downloads in the order the URLs were given
        out=\$(printf "bc_whitelist_%03d.txt" "\$i")
        wget -O "\${out}.download" "\${urls[\$i]}"

        # Unzip the downloaded file only if it is actually gzipped
        if gzip -t "\${out}.download" 2>/dev/null; then
            gunzip -c "\${out}.download" > "\$out"
        else
            mv "\${out}.download" "\$out"
        fi
        rm -f "\${out}.download"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -n1 | sed 's/^GNU Wget //; s/ .*\$//')
    END_VERSIONS
    """
}
