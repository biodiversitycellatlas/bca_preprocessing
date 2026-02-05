process CELLBENDER {
    publishDir "${params.outdir}/cellbender/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    label 'process_high_memory'

    container "oras://community.wave.seqera.io/library/cellbender:0.3.2--4f86af6695399b4f"

    input:
    tuple val(meta), path(mapping_files)

    output:
    path("cellbender_output*")

    script:
    """
    echo "\n\n===============  Ambient RNA removal  ==============="
    echo "Sample ID: ${meta}"
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
