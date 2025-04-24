process DOUBLET_DET {
    publishDir "${params.output_dir}/souporcell/${meta.id}", mode: 'symlink'
    tag "${meta.id}"
    debug true

    input:
    tuple val(meta), path(mapping_files)

    script:
    """
    echo "\n\n===============  DOUBLET DETECTION  ==============="
    echo "Sample ID: ${meta}"
    echo "Mapping files: ${mapping_files}"

    bc_file=\$(ls ${mapping_files}/Solo.out/GeneFull/raw/barcodes.tsv | head -n 1)
    bam_file=\$(ls ${mapping_files}/${meta.id}_Aligned.sortedByCoord.out.bam | head -n 1)

    echo "Barcodes file: \${bc_file}"
    echo "BAM file: \${bam_file}"
    
    ${launchDir}/submodules/souporcell/souporcell_pipeline.py \\
        --bam \${bam_file} \\
        --barcodes \${bc_file} \\
        --fasta ${params.ref_fasta} \\
        --threads 40 \\
        --out_dir . \\
        --clusters 4
    """
}
