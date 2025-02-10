// ==================  AMBIENT_RNA_RM =================== \\ 

process AMBIENT_RNA_RM {
    publishDir "${params.resDir}/ambient_rna_rm/${sample_id}", mode: 'copy'
    tag "${sample_id}"
    debug true

    input:
    tuple val(sample_id), val(config_name), path(mapping_files)

    script:
    """
    echo "\n\n===============  AMBIENT_RNA_RM ${config_name}  ==============="
    echo "Sample ID: ${sample_id}"
    
    """
}
