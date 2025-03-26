process GENE_EXT {
    publishDir "${params.resDir}/gene_ext/geneext_${config}/${sample_id}", mode: 'copy'
    tag "${sample_id}"
    debug true

    input:
    tuple val(sample_id), path(mapping_files)
    file(bam_index)
    val(config)

    output:
    path("${sample_id}*")
    
    script:
    """
    echo "\n\n==================  GENE EXTENSION ${config} =================="
    echo "Sample ID: ${sample_id}"
    echo "Mapping files: ${mapping_files}"
    echo "BAM index: ${bam_index}"
    echo "Original GTF: ${params.ref_star_gtf}"
    echo "Conda environment: \$CONDA_DEFAULT_ENV"

    if [ -d "tmp" ]; then rm -r tmp; fi
    
    if [ ${params.annot_type} == "GFF" ];
    then
        gtf_output="${sample_id}_${config}_geneext.gff"
    else
        gtf_output="${sample_id}_${config}_geneext.gtf"
    fi
    echo \${gtf_output}
    bam_file=\$(ls *_Aligned.sortedByCoord.out.bam | head -n 1)
    
    python ${params.baseDir}/submodules/GeneExt/geneext.py \\
        -g ${params.ref_star_gtf} \\
        -b \${bam_file} \\
        -o \${gtf_output} \\
        -j 4 

    """
}
