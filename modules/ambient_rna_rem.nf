// ==================  DOUBLET DETECTION  =================== \\ 

process AMBIENT_RNA_RM {
    publishDir "${params.resDir}/ambient_rna_rm/${sample_id}", mode: 'symlink'
    tag "${sample_id}"
    debug true

    input:
    tuple val(sample_id), val(config_name), path(mapping_files)

    script:
    """
    echo "\n\n===============  DOUBLET DETECTION ${config_name}  ==============="
    echo "Sample ID: ${sample_id}"
    
    """
}
