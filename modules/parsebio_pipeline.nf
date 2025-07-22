// ============  MAPPING PARSE PIPELINE  =========== \\ 
// Mapping using ParseBiosciences pipeline 1.3.1     \\
// The functions specify read 1 and read 2 from      \\
// the channel, and otherwise state if not provided. \\
// The settings are set in the variable config_file, \\
// and is specific per sequencing tech.              \\
// As all previously created reference files are     \\
// given as one string, these are added to a new     \\
// directory.                                        \\

process PARSEBIO_PIPELINE {
    publishDir "${params.output_dir}/ParseBio_pipeline/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    
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
        --output_dir . \\
        --parfile config_${params.protocol}_splitpipe.txt \\
        --no_keep_going
    """
}
