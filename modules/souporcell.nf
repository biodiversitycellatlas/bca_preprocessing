process SOUPORCELL {
    publishDir "${params.output_dir}/souporcell/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    debug true
    
    input:
        tuple val(meta), path(mapping_files)
        path bc_whitelist

    output:
        path "clusters.tsv"
        path "cluster_genotypes.vcf"
        path "ambient_rna.txt"

    script:
    """
    echo "\n\n===============  DOUBLET DETECTION  ==============="
    echo "Sample ID: ${meta}"
    echo "Mapping files: ${mapping_files}"

    bam_file=\$(ls ${meta.id}_Aligned.sortedByCoord.out.bam | head -n 1)
    echo "BAM file: \${bam_file}"

    souporcell_pipeline.py \\
        -i \${bam_file} \\
        -b ${bc_whitelist} \\
        -f ${params.ref_fasta} \\
        -t ${task.cpus} \\
        -o ./ \\
        --sample ${meta.id}
    """
}
