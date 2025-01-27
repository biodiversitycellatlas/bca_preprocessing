// ==================  GENE EXTENSION  =================== \\ 

process GENE_EXT {
    publishDir "${params.resDir}/gene_ext/gene_ext_${config_name}/${sample_id}", mode: 'copy'
    tag "${sample_id}"
    debug true

    input:
    tuple val(sample_id), val(config_name), path(mapping_files)
    file(bam_index)

    output:
    path("${sample_id}_geneext.gtf")
    
    script:
    """
    echo "\n\n==================  GENE EXTENSION ${config_name} =================="
    echo "Sample ID: ${sample_id}"
    echo "Mapping files: ${mapping_files}"
    echo "BAM index: ${bam_index}"
    echo "Original GTF: ${params.ref_star_gtf}"
    echo "Conda environment: \$CONDA_DEFAULT_ENV"

    if [ -d "tmp" ]; then rm -r tmp; fi
    
    python ${params.baseDir}/ext_programs/GeneExt/geneext.py \\
        -g ${params.ref_star_gtf} \\
        -b Aligned.sortedByCoord.out.bam \\
        -o ${sample_id}_geneext.gtf \\
        --force \\
        -j 4 \\
        --verbose 0 \\
        --keep_intermediate_files
    """
}
