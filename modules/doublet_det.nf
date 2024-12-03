// ==================  DOUBLET DETECTION  =================== \\ 

process DOUBLET_DET {
    publishDir "${params.resDir}/doublets", mode: 'symlink'
    conda '../ext_programs/souporcell/souporcell_env.yaml'
    debug true

    input:
    path(mapping_files)

    script:
    """
    echo "\n\n===============  DOUBLET DETECTION  ==============="
    echo "BAM filepath: ${mapping_files}"
    echo "BAM file: ${mapping_files}/*.bam"
    
    ../ext_programs/souporcell/souporcell_pipeline.py \\
        --bam ${mapping_files}/*.bam \\
        --barcodes ${params.barcodeDoublet} \\
        --fasta ${params.ref_fasta} \\
        --threads 1 \\
        --out_dir . \\
        --clusters 4
    """
}
