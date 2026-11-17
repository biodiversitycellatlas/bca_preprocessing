process PARSEBIO_PIPELINE_MKREF {
    label 'process_medium'

    conda "${params.splitpipe_conda_env}"

    input:
    path fasta
    path gtf

    output:
    path "ref_splitpipe", emit: reference
    path "versions.yml",  emit: versions

    script:
    """
    echo -e "\\n\\n==================  REF GENOME PARSE PIPELINE  =================="
    echo "Conda environment: \$CONDA_DEFAULT_ENV"

    # Check if a specific reference GTF file is provided, otherwise use the default
    GTF_FILE="${gtf}"

    # Create reference genome index files, in their own directory (split-pipe accepts an existing one)
    mkdir -p ref_splitpipe
    split-pipe -m mkref \\
        --genome_name ref_splitpipe \\
        --genes \$GTF_FILE \\
        --fasta ${fasta} \\
        --output_dir ref_splitpipe

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        split-pipe: \$(split-pipe --version 2>&1 | sed -n '1p')
    END_VERSIONS
    """
}
