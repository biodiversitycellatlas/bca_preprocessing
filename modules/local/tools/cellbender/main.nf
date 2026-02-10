process CELLBENDER {
    publishDir "${params.outdir}/cellbender/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    label 'process_high'
    debug true

    container "oras://community.wave.seqera.io/library/cellbender:0.3.2--4f86af6695399b4f"

    input:
    tuple val(meta), path(starsolo_genefull50_raw)

    output:
    path("cellbender_output*")

    script:
    """
    echo "\n\n===============  Ambient RNA removal  ==============="
    echo "Sample ID: ${meta}"
    echo "STARsolo files: ${starsolo_genefull50_raw}"
    echo "Number of expected cells: ${meta.expected_cells}"

    # Copy features file as cellbender expects the file to be named genes.tsv
    cp ${starsolo_genefull50_raw}/features.tsv ${starsolo_genefull50_raw}/genes.tsv

    cellbender remove-background \\
        --cpu-threads ${task.cpus} \\
        --input ${starsolo_genefull50_raw} \\
        --output cellbender_output.h5 \\
        --epochs 150 \\
        --expected-cells ${meta.expected_cells} \\
        --fpr 0.01 \\
        --low-count-threshold 50
    """
}
