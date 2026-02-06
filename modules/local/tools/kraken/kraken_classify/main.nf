process KRAKEN {
    publishDir "${params.outdir}/kraken/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    label 'process_high_memory'
    debug true

    input:
    path db_path_file
    tuple val(meta), path(mapping_files)

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/kraken2_samtools_coreutils_pigz:8961943c277652a6"

    output:
    path("*")

    script:
    """
    echo "\n\n==================  KRAKEN  =================="
    echo "Kraken db path file: ${db_path_file}"
    echo "Running KRAKEN for ${meta.id}"
    echo "Mapping files: ${mapping_files}"
    echo "BAM file: \$(ls *_Aligned.sortedByCoord.out.bam)"

    kraken_db_path=\$(cat ${db_path_file})
    echo "Using Kraken2 DB at: \$kraken_db_path"

    # Saving unmapped reads to a fasta file
    # TODO: Move this step to a separate process or mapping_starsolo process
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
