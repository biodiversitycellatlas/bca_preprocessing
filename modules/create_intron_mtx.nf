process CREATE_INTRON_MTX {   
    input:
    tuple val(meta), path(mapping_files)  
       
    script:
    """
    echo "\n\n==================  CREATE INTRONIC COUNT MATRIX  =================="
    echo "Sample ID: ${meta}"
    echo "Mapping files: ${mapping_files}"

    mapping_dir="${params.output_dir}/mapping_STARsolo/${meta.id}"
    sbatch ${launchDir}/bin/create_intron_mtx.sh \${mapping_dir}
    """
}
