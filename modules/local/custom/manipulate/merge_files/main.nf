process MERGE_FASTQS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fastq_cDNA_list), path(fastq_BC_UMI_list), path(fastq_indices), path(input_file)

    output:
    tuple val(meta), path("${meta.id}_merged_cDNA.fastq.gz"), path("${meta.id}_merged_BC_UMI.fastq.gz"), path("${meta.id}_merged_indices/"), path(input_file)

    script:
    """
    echo "================== MERGE FASTQs =================="
    echo "Sample ID: ${meta.id}"

    COMPRESSOR="pigz -p ${task.cpus}"
    CDNA=( ${fastq_cDNA_list} )
    BCUMI=( ${fastq_BC_UMI_list} )
    IND=( ${fastq_indices} )

    echo "Counts: cDNA=\${#CDNA[@]} BC_UMI=\${#BCUMI[@]} Index=\${#IND[@]}"

    if [ \${#CDNA[@]} -eq 1 ]; then
        ln -s "\${CDNA[0]}" "${meta.id}_merged_cDNA.fastq.gz"
    else
        zcat "\${CDNA[@]}" | \$COMPRESSOR -c > "${meta.id}_merged_cDNA.fastq.gz"
    fi

    if [ \${#BCUMI[@]} -eq 1 ]; then
        ln -s "\${BCUMI[0]}" "${meta.id}_merged_BC_UMI.fastq.gz"
    else
        zcat "\${BCUMI[@]}" | \$COMPRESSOR -c >"${meta.id}_merged_BC_UMI.fastq.gz"
    fi

    # Index FASTQs (commonly I1 and I2)
    mkdir -p ${meta.id}_merged_indices

    if [ \${#IND[@]} -eq 0 ]; then
        echo "[3/3] No index FASTQs detected — skipping."
    else
        echo "[3/3] Handling index FASTQs..."
        I1_LIST=(\$(echo ${fastq_indices} | tr ' ' '\\n' | grep -E 'I1|_R3|_index1' || true))
        I2_LIST=(\$(echo ${fastq_indices} | tr ' ' '\\n' | grep -E 'I2|_R4|_index2' || true))

        # I1
        if [ \${#I1_LIST[@]} -eq 1 ]; then
            ln -s \${I1_LIST[0]} ${meta.id}_merged_indices/${meta.id}_merged_I1.fastq.gz
        elif [ \${#I1_LIST[@]} -gt 1 ]; then
            zcat \${I1_LIST[@]} | \$COMPRESSOR -c > ${meta.id}_merged_indices/${meta.id}_merged_I1.fastq.gz
        fi

        # I2
        if [ \${#I2_LIST[@]} -eq 1 ]; then
            ln -s \${I2_LIST[0]} ${meta.id}_merged_indices/${meta.id}_merged_I2.fastq.gz
        elif [ \${#I2_LIST[@]} -gt 1 ]; then
            zcat \${I2_LIST[@]} | \$COMPRESSOR -c > ${meta.id}_merged_indices/${meta.id}_merged_I2.fastq.gz
        fi
    fi

    echo "Merging completed for ${meta.id}"
    """
}
