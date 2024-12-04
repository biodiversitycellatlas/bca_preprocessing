// ==================  GENE EXTENSION  =================== \\ 

process GENE_EXT {
    publishDir "${params.resDir}/gene_ext_${config_name}", mode: 'symlink'
    conda '/users/asebe/bvanwaardenburg/miniconda3/envs/geneext'
    tag "${sample_id}"
    debug true

    input:
    tuple val(sample_id), val(config_name), path(mapping_files)

    output:
    path("*")
    
    script:
    """
    echo "\n\n==================  GENE EXTENSION ${config_name} =================="
    echo "Mapping files: ${mapping_files}"
    echo "Original GTF: ${params.ref_star_gtf}"

    bam_file=\$(ls *Aligned.sortedByCoord.out.bam | head -n 1)
    echo "BAM file: \${bam_file}"
    
    python ${params.baseDir}/ext_programs/GeneExt/geneext.py \\
        -g ${params.ref_star_gtf} \\
        -b \${bam_file} \\
        -o result.gtf \\
        --peak_perc 0
    """
}
