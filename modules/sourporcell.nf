process DOUBLET_DET {
    publishDir "${params.resDir}/souporcell/${sample_id}", mode: 'symlink'
    tag "${sample_id}"
    debug true

    input:
    tuple val(sample_id), path(mapping_files)

    script:
    """
    echo "\n\n===============  DOUBLET DETECTION  ==============="
    echo "Sample ID: ${sample_id}"
    echo "Mapping files: ${mapping_files}"

    bc_file=\$(ls ${mapping_files}/Solo.out/GeneFull/raw/barcodes.tsv | head -n 1)
    bam_file=\$(ls ${mapping_files}/${sample_id}_Aligned.sortedByCoord.out.bam | head -n 1)

    echo "Barcodes file: \${bc_file}"
    echo "BAM file: \${bam_file}"
    
    ${params.baseDir}/submodules/souporcell/souporcell_pipeline.py \\
        --bam \${bam_file} \\
        --barcodes \${bc_file} \\
        --fasta ${params.ref_fasta} \\
        --threads 40 \\
        --out_dir . \\
        --clusters 4
    """
}
