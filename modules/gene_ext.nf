process GENE_EXT {
    publishDir "${params.output_dir}/gene_ext/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    debug true

    input:
    tuple val(meta), path(mapping_files)
    file(bam_index)

    output:
    path("${meta.id}*.geneext.g{tf, ff}")
    
    script:
    """
    echo "\n\n==================  GENE EXTENSION =================="
    echo "Sample ID: ${meta}"
    echo "Mapping files: ${mapping_files}"
    echo "BAM index: ${bam_index}"
    echo "Original GTF: ${params.ref_gtf}"
    echo "Conda environment: \$CONDA_DEFAULT_ENV"

    if [ -d "tmp" ]; then rm -r tmp; fi
    
    if [ ${params.annot_type} == "GFF" ];
    then
        gtf_output="${meta.id}_geneext.gff"
    else
        gtf_output="${meta.id}_geneext.gtf"
    fi
    echo \${gtf_output}
    bam_file=\$(ls *_Aligned.sortedByCoord.out.bam | head -n 1)
    
    python ${launchDir}/submodules/GeneExt/geneext.py \\
        -g ${params.ref_gtf} \\
        -b \${bam_file} \\
        -o \${gtf_output} \\
        -j 4 

    """
}
