// ============= CREATE INTRONIC COUNT MATRIX ============= \\ 
// 

process CREATE_INTRON_MTX {   
    input:
    tuple val(sample_id), val(config_name), path(mapping_files)  
       
    script:
    """
    echo "\n\n==================  CREATE INTRONIC COUNT MATRIX - ${config_name} =================="
    echo "Sample ID: ${sample_id}"
    echo "Config name: ${config_name}"
    echo "Mapping files: ${mapping_files}"

    mapping_dir="${params.resDir}/mapping_STARsolo/mapping_STARsolo_${config_name}/${sample_id}"
    sbatch ${params.baseDir}/scripts/create_intron_mtx.sh \${mapping_dir}
    """
}
