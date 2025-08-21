process KRAKEN {
    publishDir "${params.output_dir}/kraken/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    label 'process_high_memory'

    input:
    path db_path_file 
    tuple val(meta), path(mapping_files)
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/29/29ed8f68315625eca61a3de9fcb7b8739fe8da23f5779eda3792b9d276aa3b8f/data' :
        'community.wave.seqera.io/library/kraken2_coreutils_pigz:45764814c4bb5bf3' }"

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
