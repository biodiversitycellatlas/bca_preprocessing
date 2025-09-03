process PARSEBIO_PIPELINE {
    publishDir "${params.outdir}/ParseBio_pipeline/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    label 'process_high'

    conda "${params.splitpipe_conda_env}"
    
    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI)
    path parse_refgenome_files

    output:
    path("*")

    script:
    def kitskip_arg   = task.ext.args ?: ''          // If ext.args is defined assign it to kitskip_arg
    """
    echo "\n\n=============  PARSE BIOSCIENCES PIPELINE  ================"
    echo "Mapping sample ${meta.id} with Parse Biosciences pipeline"
    echo "FASTQ cDNA: ${fastq_cDNA}"
    echo "FASTQ BC & UMI: ${fastq_BC_UMI}"
    echo "Genome index directory: ${parse_refgenome_files}"
    echo "Conda environment: \$CONDA_DEFAULT_ENV"

    # Generate parameter file
    echo "post_min_map_frac 0.01" > config_${params.protocol}_splitpipe.txt

    # Create directory for the genome index files
    mkdir -p genome_index
    
    # Move all genome index files to the new directory
    mv ${parse_refgenome_files} genome_index/
        
    split-pipe -m all \\
        --chemistry v3 \\
        --fq1 ${fastq_cDNA} \\
        --fq2 ${fastq_BC_UMI} \\
        --nthreads 16 \\
        --genome_dir genome_index \\
        --outdir . \\
        --parfile config_${params.protocol}_splitpipe.txt \\
        --no_keep_going
    """
}
