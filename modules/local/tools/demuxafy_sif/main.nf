process DOWNLOAD_DEMUXAFY_SIF {
    publishDir "${params.outdir}/containers/demuxafy", mode: 'copy', overwrite: false
    tag "demuxafy_sif"
    label 'process_single2'

    output:
    path 'Demuxafy.sif', emit: sif_file
    path 'demuxafy_sif_path.txt', emit: sif_path_file

    script:
    // dl=1 makes Dropbox serve the raw file instead of the HTML preview page
    def sif_url = 'https://www.dropbox.com/scl/fi/kykwi78vk4yifbbag5ajz/Demuxafy.sif?rlkey=5hcugu6ztpy0eik3xno63xiar&dl=1'
    def md5_url = 'https://www.dropbox.com/scl/fi/37oj9y1frzhqazl4h8s21/Demuxafy.sif.md5?rlkey=o2bn5wp9q68numlaav8gg95kh&dl=1'
    """
    echo "\n\n==================  DOWNLOAD DEMUXAFY SIF =================="

    if [ -z "${params.demuxafy_sif}" ] || [ "${params.demuxafy_sif}" = "null" ]; then
        echo "No demuxafy_sif provided. Downloading Demuxafy.sif (~7.5GB, this can take 15-30 min)..."
        wget -q -O Demuxafy.sif '${sif_url}'
        wget -q -O Demuxafy.sif.md5 '${md5_url}'

        # Guard against an HTML page being saved instead of the actual files
        if grep -qi '<html' Demuxafy.sif.md5 || head -c 1024 Demuxafy.sif | grep -qi '<html'; then
            echo "ERROR: Download returned an HTML page instead of Demuxafy.sif; check the download URLs or provide --demuxafy_sif" >&2
            exit 1
        fi

        expected_md5=\$(grep -o '[a-f0-9]\\{32\\}' Demuxafy.sif.md5 | head -1)
        actual_md5=\$(md5sum Demuxafy.sif | awk '{print \$1}')
        if [ "\$expected_md5" != "\$actual_md5" ]; then
            echo "ERROR: MD5 checksum mismatch for Demuxafy.sif (expected \$expected_md5, got \$actual_md5)" >&2
            exit 1
        fi

        # Record where the image ends up, so it can be passed to --demuxafy_sif on a later run.
        echo "${params.outdir}/containers/demuxafy/Demuxafy.sif" > demuxafy_sif_path.txt
    else
        echo "Using existing Demuxafy.sif: ${params.demuxafy_sif}"
        echo "${params.demuxafy_sif}" > demuxafy_sif_path.txt
        # If params.demuxafy_sif is provided, copy it to the work directory
        cp "${params.demuxafy_sif}" Demuxafy.sif
    fi
    """
}
