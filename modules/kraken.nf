process KRAKEN {
    publishDir "${params.resDir}/kraken/${sample_id}", mode: 'copy'
    tag "${sample_id}"

    input:
    path db_path_file 
    tuple val(sample_id), path(mapping_files)

    script:
    """
    echo "\n\n==================  KRAKEN  =================="
    echo "Running KRAKEN for ${sample_id}"
    echo "Path: ${mapping_files}"

    kraken_db_path=\$(cat ${db_path_file})
    echo "Using Kraken2 DB at: \$kraken_db_path"

    # Saving unmapped reads to a new bam file, where the 0x4 flag specifies unmapped reads
    samtools view -f 0x4 *_Aligned.sortedByCoord.out.bam > unmapped.sortedByCoord.out.bam

    # Saving unmapped reads to a fasta file
    samtools view -f 0x4 *_Aligned.sortedByCoord.out.bam | samtools fasta - > unmapped.fasta

    # Run Kraken2
    k2 classify \
        --threads 8 \
        --db \${kraken_db_path} \
        --report kraken_taxonomy.txt \
        --report-minimizer-data \
        --use-names \
        --memory-mapping \
        --log kraken.log \
        --output kraken_output.txt \
        ${sample_id}_unmapped.fasta

    """
}
