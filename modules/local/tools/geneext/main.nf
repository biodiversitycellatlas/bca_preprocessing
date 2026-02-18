process GENE_EXT {
    publishDir "${params.outdir}/gene_ext/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    label 'process_medium'

    // conda "${projectDir}/submodules/GeneExt/environment.yaml"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(mapping_files)
    file(bam_index)

    output:
    path("${meta.id}*_geneext.g{tf, ff}")

    script:
    """
    echo "\n\n==================  GENE EXTENSION =================="
    echo "Sample ID: ${meta}"
    echo "BAM index: ${bam_index}"
    echo "Original GTF: ${params.ref_gtf}"

    # Remove temporary directory if it exists
    if [ -d "tmp" ]; then rm -r tmp; fi

    # Extract file extension
    extension=\$(echo "${params.ref_gtf}" | awk -F. '{print \$NF}')

    if [ \$extension == "gff" ];
    then
        gtf_output="${meta.id}_geneext.gff"
    else
        gtf_output="${meta.id}_geneext.gtf"
    fi
    echo \${gtf_output}
    bam_file=\$(ls *_Aligned.sortedByCoord.out.bam | head -n 1)

    # Run GeneExt
    python ${projectDir}/submodules/GeneExt/geneext.py \\
        -g ${params.ref_gtf} \\
        -b \${bam_file} \\
        -o \${gtf_output} \\
        -j 4

    """
}
