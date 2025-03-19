process CELLBENDER {
    publishDir "${params.resDir}/Cellbender/${sample_id}", mode: 'copy'
    tag "${sample_id}"
    debug true

    input:
    tuple val(sample_id), val(config_name), path(mapping_files)

    script:
    """
    echo "\n\n===============  Ambient RNA removal ${config_name}  ==============="
    echo "Sample ID: ${sample_id}"
    echo "Mapping files: ${mapping_files}"
    
    # Copy features file as cellbender expects the file to be named genes.tsv
    cp ${mapping_files}/features.tsv ${mapping_files}/genes.tsv

    cellbender remove-background \
        --input ${mapping_files}/Solo.out/Gene/raw/ \
        --output cellbender_output.h5 
    """
}