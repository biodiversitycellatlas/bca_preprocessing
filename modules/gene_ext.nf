// ==================  GENE EXTENSION  =================== \\ 

process GENE_EXT {
    publishDir "${params.resDir}/gene_ext_${config_name}", mode: 'symlink'
    tag "${sample_id}"
    debug true

    input:
    tuple val(sample_id), val(config_name), path(mapping_files)
    tuple val(sample_id2), val(config_name2), path(indexed_bam)

    output:
    path("result.gtf")
    
    script:
    """
    echo "\n\n==================  GENE EXTENSION ${config_name} =================="
    echo "Mapping files: ${mapping_files}"
    echo "Original GTF: ${params.ref_star_gff}"

    bam_file=\$(ls ${mapping_files}/*Aligned.sortedByCoord.out.bam | head -n 1)
    echo "BAM file: \${bam_file}"
    
    # rm -r ${params.baseDir}/ext_programs/GeneExt/tmp
    # rm ${params.baseDir}/ext_programs/GeneExt/result*
    
    python ${params.baseDir}/ext_programs/GeneExt/geneext.py \\
        -g ${params.ref_star_gff} \\
        -b \${bam_file} \\
        -o result.gtf \\
        --peak_perc 0
    """
}
