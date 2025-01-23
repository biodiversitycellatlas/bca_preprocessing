// ==================  GENE EXTENSION  =================== \\ 

process GENE_EXT {
    publishDir "${params.resDir}/genome/gene_ext/${sample_id}_${config_name}", mode: 'copy'
    tag "${sample_id}"
    debug true

    input:
    tuple val(sample_id), val(config_name), path(mapping_files)

    output:
    path("${sample_id}_geneext.gtf")
    
    script:
    """
    echo "\n\n==================  GENE EXTENSION ${config_name} =================="
    echo "Sample ID: ${sample_id}"
    echo "Mapping files: ${mapping_files}"
    echo "Original GTF: ${params.ref_star_gtf}"
    echo "Conda environment: \$CONDA_DEFAULT_ENV"

    if [ -d "tmp" ]; then rm -r tmp; fi
    if [ -d "${params.baseDir}/ext_programs/GeneExt/tmp" ]; then rm -r ${params.baseDir}/ext_programs/GeneExt/tmp; fi
    # if [ -d "${params.baseDir}/ext_programs/GeneExt/result*" ]; then rm -r ${params.baseDir}/ext_programs/GeneExt/result*; fi
    
    python ${params.baseDir}/ext_programs/GeneExt/geneext.py \\
        -g ${params.ref_star_gtf} \\
        -b Aligned.sortedByCoord.out.bam \\
        -o ${sample_id}_geneext.gtf \\
        -j 40

    # Wait for output file to exist
    while [ ! -f ${sample_id}_geneext.gtf ]; do
        echo "Waiting for output file..."
        sleep 10
    done
    """
}
