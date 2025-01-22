// ==================  GENE EXTENSION  =================== \\ 

process GENE_EXT {
    publishDir "${params.resDir}/genome/gene_ext/${sample_id}_${config_name}", mode: 'copy'
    tag "${sample_id}"
    debug true

    input:
    tuple val(sample_id), val(config_name), path(mapping_files)

    output:
    path("result.gtf")
    
    script:
    """
    echo "\n\n==================  GENE EXTENSION ${config_name} =================="
    echo "Sample ID: ${sample_id}"
    echo "Mapping files: ${mapping_files}"
    echo "Original GTF: ${params.ref_star_gff}"

    conda list 

    which python
    python -c "import gffutils; print('gffutils found!')"

    # bam_file=\$(ls Aligned.sortedByCoord.out.bam | head -n 1)
    # echo "BAM file: \${bam_file}"
    
    if [ -d "tmp" ]; then rm -Rf tmp; fi
    if [ -d "${params.baseDir}/ext_programs/GeneExt/tmp" ]; then rm -Rf ${params.baseDir}/ext_programs/GeneExt/tmp; fi
    if [ -d "${params.baseDir}/ext_programs/GeneExt/result*" ]; then rm -Rf ${params.baseDir}/ext_programs/GeneExt/result*; fi
    
    python ${params.baseDir}/ext_programs/GeneExt/geneext.py \\
        -g /users/asebe/bvanwaardenburg/git/data/240810_ParseBio_Nvec_Tcas/Nvec_BCA001_BCA002_TestREF/genome/Nvec_v5_merged_annotation_sort.gtf \\ 
        -b Aligned.sortedByCoord.out.bam \\
        -o result.gtf \\
        --peak_perc 0
    """
}
