process MERGE_FASTQS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fastq_cDNA_list), path(fastq_BC_UMI_list), path(fastq_indices), path(input_file)

    output:
    tuple val(meta), path("${meta.id}_merged_cDNA.fastq.gz"), path("${meta.id}_merged_BC_UMI.fastq.gz"), path("${meta.id}_merged_I*.fastq.gz", arity: '0..2'), path(input_file), emit: merged_files
    path "versions.yml", emit: versions

    script:
    """
    echo "================== MERGE FASTQs =================="
    echo "Sample ID: ${meta.id}"

    COMPRESSOR="pigz -p ${task.cpus}"
    CDNA=( ${fastq_cDNA_list} )
    BCUMI=( ${fastq_BC_UMI_list} )
    IND=( ${fastq_indices} )

    echo "Counts: cDNA=\${#CDNA[@]} BC_UMI=\${#BCUMI[@]} Index=\${#IND[@]}"

    # Inputs may be plain (.fastq/.fq) or gzipped (.fastq.gz/.fq.gz), and a sample may mix both. 
    # Compression is detected from the gzip magic bytes rather than the file name, so mislabelled extensions are handled too. 
    # Whatever comes in, the merged output is always gzipped, which is what every downstream module expects.
    is_gzipped() {
        [ "\$(od -An -N2 -tx1 "\$1" | tr -d ' \\n')" = "1f8b" ]
    }

    # Stream the plain-text contents of one or more FASTQs to stdout
    cat_fastqs() {
        for f in "\$@"; do
            if is_gzipped "\$f"; then
                pigz -dc -p ${task.cpus} "\$f"
            else
                cat "\$f"
            fi
        done
    }

    # A lone gzipped input is symlinked instead of recompressed.
    merge_fastqs() {
        local out=\$1
        shift
        if [ "\$#" -eq 1 ] && is_gzipped "\$1"; then
            ln -s "\$1" "\$out"
        else
            cat_fastqs "\$@" | \$COMPRESSOR -c > "\$out"
        fi
    }

    merge_fastqs "${meta.id}_merged_cDNA.fastq.gz" "\${CDNA[@]}"
    merge_fastqs "${meta.id}_merged_BC_UMI.fastq.gz" "\${BCUMI[@]}"

    # Index FASTQs (commonly I1 and I2)
    mkdir -p ${meta.id}_merged_indices

    if [ \${#IND[@]} -eq 0 ]; then
        echo "[3/3] No index FASTQs detected — skipping."
    else
        echo "[3/3] Handling index FASTQs..."
        I1_LIST=(\$(echo ${fastq_indices} | tr ' ' '\\n' | grep -E 'I1|_R3|_index1' || true))
        I2_LIST=(\$(echo ${fastq_indices} | tr ' ' '\\n' | grep -E 'I2|_R4|_index2' || true))

        # I1
        if [ \${#I1_LIST[@]} -gt 0 ]; then
            merge_fastqs "${meta.id}_merged_I1.fastq.gz" "\${I1_LIST[@]}"
        fi

        # I2
        if [ \${#I2_LIST[@]} -gt 0 ]; then
            merge_fastqs "${meta.id}_merged_I2.fastq.gz" "\${I2_LIST[@]}"
        fi
    fi

    echo "Merging completed for ${meta.id}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$(pigz --version 2>&1 | sed 's/pigz //')
    END_VERSIONS
    """
}
