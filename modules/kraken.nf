process KRAKEN {
    publishDir "${params.output_dir}/kraken/${sample_id}", mode: 'copy'
    tag "${sample_id}"

    input:
    path db_path_file 
    tuple val(sample_id), path(mapping_files)

    output:
    path("*")     

    script:
    """
    echo "\n\n==================  KRAKEN  =================="
    echo "Conda environment: \$CONDA_DEFAULT_ENV"
    echo "Running KRAKEN for ${sample_id}"
    echo "Path: ${mapping_files}"

    kraken_db_path=\$(cat ${db_path_file})
    echo "Using Kraken2 DB at: \$kraken_db_path"

    # Saving unmapped reads to a fasta file
    samtools view -f 0x4 -b *_Aligned.sortedByCoord.out.bam | samtools fasta - > ${sample_id}_unmapped.fasta

    # Run Kraken2
    k2 classify \\
        --threads 8 \\
        --db \${kraken_db_path} \\
        --report ${sample_id}_kraken_taxonomy.txt \\
        --report-minimizer-data \\
        --use-names \\
        --memory-mapping \\
        --log ${sample_id}_kraken.log \\
        --output ${sample_id}_kraken_output.txt \\
        ${sample_id}_unmapped.fasta

    """
}
