process CELLBENDER {
    publishDir "${params.resDir}/cellbender/${sample_id}", mode: 'copy'
    tag "${sample_id}"
    debug true

    input:
    tuple val(sample_id), path(mapping_files)

    script:
    """
    echo "\n\n===============  Ambient RNA removal  ==============="
    echo "Sample ID: ${sample_id}"
    echo "Mapping files: ${mapping_files}"
    
    # Copy features file as cellbender expects the file to be named genes.tsv
    cp ${mapping_files}/features.tsv ${mapping_files}/genes.tsv

    cellbender remove-background \\
        --input ${mapping_files}/Solo.out/Gene/raw/ \\
        --output cellbender_output.h5 \\
        --epochs 150 \\
        --expected-cells 2500 \\
        --fpr 0.01 \\
        --cuda
    """
}