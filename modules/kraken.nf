process KRAKEN {
    publishDir "${params.output_dir}/kraken/${meta.id}", mode: 'copy'
    tag "${meta.id}"

    input:
    path db_path_file 
    tuple val(meta), path(mapping_files)

    output:
    path("*")     

    script:
    """
    echo "\n\n==================  KRAKEN  =================="
    echo "Conda environment: \$CONDA_DEFAULT_ENV"
    echo "Running KRAKEN for ${meta}"
    echo "Path: ${mapping_files}"

    kraken_db_path=\$(cat ${db_path_file})
    echo "Using Kraken2 DB at: \$kraken_db_path"

    # Saving unmapped reads to a fasta file
    samtools view -f 0x4 -b *_Aligned.sortedByCoord.out.bam | samtools fasta - > ${meta.id}_unmapped.fasta

    # Run Kraken2
    k2 classify \\
        --threads 8 \\
        --db \${kraken_db_path} \\
        --report ${meta.id}_kraken_taxonomy.txt \\
        --report-minimizer-data \\
        --use-names \\
        --memory-mapping \\
        --log ${meta.id}_kraken.log \\
        --output ${meta.id}_kraken_output.txt \\
        ${meta.id}_unmapped.fasta

    """
}
