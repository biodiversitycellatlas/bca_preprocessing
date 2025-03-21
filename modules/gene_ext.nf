process GENE_EXT {
    publishDir "${params.resDir}/gene_ext/gene_ext_${config_name}/${sample_id}", mode: 'copy'
    tag "${sample_id}"
    debug true

    input:
    tuple val(sample_id), val(config_name), path(mapping_files)
    file(bam_index)

    output:
    path("${sample_id}*")
    
    script:
    """
    echo "\n\n==================  GENE EXTENSION ${config_name} =================="
    echo "Sample ID: ${sample_id}"
    echo "Mapping files: ${mapping_files}"
    echo "BAM index: ${bam_index}"
    echo "Original GTF: ${params.ref_star_gtf}"
    echo "Conda environment: \$CONDA_DEFAULT_ENV"

    if [ -d "tmp" ]; then rm -r tmp; fi
    
    if [ ${params.annot_type} == "GFF" ];
    then
        gtf_input="reference.gff"
        gtf_output="${sample_id}_${config_name}_geneext.gff"
        cp ${params.ref_star_gtf} \${gtf_input}
    else
        gtf_input="reference.gtf"
        gtf_output="${sample_id}_${config_name}_geneext.gtf"
        cp ${params.ref_star_gtf} \${gtf_input}
    fi
    echo \${gtf_output}
    bam_file=\$(ls ${sample_id}_${config_name}_Aligned.sortedByCoord.out.bam | head -n 1)
    
    python ${params.baseDir}/submodules/GeneExt/geneext.py \\
        -g \${gtf_input} \\
        -b \${bam_file} \\
        -o \${gtf_output} \\
        --force \\
        -j 4 
    """
}
