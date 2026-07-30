process DOWNLOAD_DEMUXAFY_SIF {
    publishDir "${params.outdir}/containers/demuxafy", mode: 'copy', overwrite: false
    tag "demuxafy_sif"
    label 'process_low'

    output:
    path 'demuxafy_sif_path.txt', emit: sif_path_file

    script:
    def sif_url = 'https://www.dropbox.com/scl/fi/kykwi78vk4yifbbag5ajz/Demuxafy.sif?rlkey=5hcugu6ztpy0eik3xno63xiar'
    def md5_url = 'https://www.dropbox.com/scl/fi/37oj9y1frzhqazl4h8s21/Demuxafy.sif.md5?rlkey=o2bn5wp9q68numlaav8gg95kh'
    """
    echo "\n\n==================  DOWNLOAD DEMUXAFY SIF =================="

    if [ -z "${params.demuxafy_sif}" ]; then
        echo "No demuxafy_sif provided. Downloading Demuxafy.sif (~7.5GB, this can take 15-30 min)..."
        wget -q -O Demuxafy.sif '${sif_url}'
        wget -q -O Demuxafy.sif.md5 '${md5_url}'

        expected_md5=\$(grep -o '[a-f0-9]\\{32\\}' Demuxafy.sif.md5 | head -1)
        actual_md5=\$(md5sum Demuxafy.sif | awk '{print \$1}')
        if [ "\$expected_md5" != "\$actual_md5" ]; then
            echo "ERROR: MD5 checksum mismatch for Demuxafy.sif (expected \$expected_md5, got \$actual_md5)" >&2
            exit 1
        fi

        echo "\$(pwd)/Demuxafy.sif" > demuxafy_sif_path.txt
    else
        echo "Using existing Demuxafy.sif: ${params.demuxafy_sif}"
        echo "${params.demuxafy_sif}" > demuxafy_sif_path.txt
    fi
    """
}
