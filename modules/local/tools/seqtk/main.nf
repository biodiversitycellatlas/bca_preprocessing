process SUBSAMPLE_FASTQS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file)

    output:
    tuple val(meta), path("${meta.id}_subsampled_cDNA.fastq.gz"), path("${meta.id}_subsampled_BC_UMI.fastq.gz"), path("${meta.id}_subsampled_indices/"), path(input_file), emit: subsampled_files
    path "versions.yml", emit: versions

    script:
    """
    echo "================== SUBSAMPLE FASTQs =================="
    echo "Sample ID: ${meta.id}"
    echo "Subsampling to ${params.subsample_nreads} reads"

    SEED=100
    NREADS=${params.subsample_nreads}
    # COMPRESSOR="pigz -p ${task.cpus}"
    COMPRESSOR="gzip"

    echo "[1/2] Subsampling paired reads..."
    seqtk sample -s\$SEED ${fastq_cDNA} \$NREADS | \$COMPRESSOR > ${meta.id}_subsampled_cDNA.fastq.gz
    seqtk sample -s\$SEED ${fastq_BC_UMI} \$NREADS | \$COMPRESSOR > ${meta.id}_subsampled_BC_UMI.fastq.gz

    mkdir -p ${meta.id}_subsampled_indices

    if [ -n "${fastq_indices}" ]; then
        echo "[2/2] Subsampling index reads..."
        for idx in ${fastq_indices}; do
            base=\$(basename \$idx .fastq.gz)
            seqtk sample -s\$SEED \$idx \$NREADS | \$COMPRESSOR > ${meta.id}_subsampled_indices/${meta.id}_subsampled_\${base}.fastq.gz
        done
    else
        echo "No index FASTQs provided. Skipping index subsampling."
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1 | grep Version | sed 's/Version: //'))
    END_VERSIONS
    """
}
