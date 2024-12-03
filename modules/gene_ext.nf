// ==================  GENE EXTENSION  =================== \\ 

process GENE_EXT {
    publishDir "${params.resDir}/gene_ext_${config_name}", mode: 'symlink'
    conda '/users/asebe/bvanwaardenburg/miniconda3/envs/geneext'
    debug true

    input:
    tuple val(sample_id), val(config_name), path(mapping_files)

    output:
    tuple val(sample_id), val(config_name), path("result.gtf")
    
    script:
    """
    echo "\n\n==================  GENE EXTENSION ${config_name} =================="
    echo "BAM filepath: ${mapping_files}"
    echo "Original GTF: ${params.ref_star_gtf}"
    
    python ../ext_programs/GeneExt/geneext.py \\
        -g ${params.ref_star_gtf} \\
        -b ${mapping_files}/*.bam \\
        -o result.gtf \\
        --peak_perc 0
    """
}
