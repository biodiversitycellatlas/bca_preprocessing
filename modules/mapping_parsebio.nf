// ============  MAPPING PARSE PIPELINE  =========== \\ 
// Mapping using ParseBiosciences pipeline 1.3.1     \\
// The functions specify read 1 and read 2 from      \\
// the channel, and otherwise state if not provided. \\
// The settings are set in the variable config_file, \\
// and is specific per sequencing tech.              \\
// As all previously created reference files are     \\
// given as one string, these are added to a new     \\
// directory.                                        \\

process MAPPING_PARSEBIO {
    publishDir "${params.resDir}/mapping_parsepipe/${sample_id}", mode: 'copy', overwrite: false
    tag "${fastq_files}"
    
    input:
    tuple val(sample_id), path(fastq_files)
    path parse_refgenome_files

    output:
    path("*")

    script:
    def fastq_list = fastq_files instanceof List ? fastq_files : [fastq_files]
    def r1_fastq = fastq_list.find { it.name.contains('_R1') }
    def r2_fastq = fastq_list.find { it.name.contains('_R2') }
    """
    echo "\n\n=============  MAPPING PARSE BIOSCIENCES  ================"
    echo "Mapping sample ${sample_id} with Parse Biosciences pipeline"
    echo "FQ 1: ${r1_fastq ?: 'Not provided'}"
    echo "FQ 2: ${r2_fastq ?: 'Not provided'}"
    echo "Genome index directory: ${parse_refgenome_files}"

    # Generate parameter file
    echo "post_min_map_frac 0.01" > config_${params.seqTech}_parsepipe.txt

    # Create directory for the genome index files
    mkdir -p genome_index
    
    # Move all genome index files to the new directory
    mv ${parse_refgenome_files} genome_index/
        
    split-pipe -m all \\
        --chemistry v3 \\
        --fq1 ${r1_fastq} \\
        --fq2 ${r2_fastq} \\
        --nthreads 16 \\
        --genome_dir genome_index \\
        --output_dir . \\
        --parfile config_${params.seqTech}_parsepipe.txt
    """
}
