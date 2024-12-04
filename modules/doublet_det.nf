// ==================  DOUBLET DETECTION  =================== \\ 

process DOUBLET_DET {
    publishDir "${params.resDir}/doublets", mode: 'symlink'
    conda '${params.codeDir}/ext_programs/souporcell/souporcell_env.yaml'
    debug true

    input:
    path(mapping_files)

    script:
    """
    echo "\n\n===============  DOUBLET DETECTION  ==============="
    echo "Mapping files: ${mapping_files}"

    bc_file=\$(ls *_Solo.out/GeneFull/raw/barcodes.tsv | head -n 1)
    bam_file=\$(ls *Aligned.sortedByCoord.out.bam | head -n 1)
    
    echo "Barcodes file: \${bc_file}"
    echo "BAM file: \${bam_file}"
    
    ${params.codeDir}/ext_programs/souporcell/souporcell_pipeline.py \\
        --bam \${bam_file} \\
        --barcodes \${bc_file} \\
        --fasta ${params.ref_fasta} \\
        --threads 1 \\
        --out_dir . \\
        --clusters 4
    """
}
