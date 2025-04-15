process CREATE_INTRON_MTX {   
    input:
    tuple val(sample_id), path(mapping_files)  
       
    script:
    """
    echo "\n\n==================  CREATE INTRONIC COUNT MATRIX  =================="
    echo "Sample ID: ${sample_id}"
    echo "Mapping files: ${mapping_files}"

    mapping_dir="${params.output_dir}/mapping_STARsolo/${sample_id}"
    sbatch ${params.code_dir}/bin/create_intron_mtx.sh \${mapping_dir}
    """
}
