process CELLBENDER {
    publishDir "${params.output_dir}/cellbender/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    label 'process_high_memory'
    debug true

    conda "${moduleDir}/environment.yml"
    container "${ task.ext.use_gpu ? 'us.gcr.io/broad-dsde-methods/cellbender:0.3.2' :
        workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/eb/ebcf140f995f79fcad5c17783622e000550ff6f171771f9fc4233484ee6f63cf/data':
        'community.wave.seqera.io/library/cellbender_webcolors:156d413fdfc16cdb' }"

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