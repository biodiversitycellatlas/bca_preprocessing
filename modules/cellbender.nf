process CELLBENDER {
    publishDir "${params.resDir}/cellbender/${sample_id}", mode: 'copy'
    tag "${sample_id}"
    debug true

    input:
    tuple val(sample_id), path(mapping_files)

    script:
    """
    echo "\n\n===============  Ambient RNA removal  ==============="
    echo "Conda environment: \$CONDA_DEFAULT_ENV"
    echo "Sample ID: ${sample_id}"
    echo "Mapping files: ${mapping_files}"
    
    matrix_path=\$(echo ./*_Solo.out/GeneFull_Ex50pAS/raw)

    # Copy features file as cellbender expects the file to be named genes.tsv
    cp \${matrix_path}/features.tsv \${matrix_path}/genes.tsv

    cellbender remove-background \\
        --input \${matrix_path} \\
        --output cellbender_output.h5 \\
        --epochs 150 \\
        --expected-cells 2500 \\
        --fpr 0.01
    """
}