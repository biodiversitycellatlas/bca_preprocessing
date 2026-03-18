process CELLSWEEP {
    publishDir "${params.outdir}/cellsweep/${meta.id}/${mapping_method}", mode: 'copy'
    tag "${meta.id} | ${mapping_method}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(raw_h5ad), val(mapping_method), val(datatype)

    output:
    path("${meta.id}_${mapping_method}_cs_filtered.h5ad"),      emit: cs_filtered_h5ad
    path("${meta.id}_${mapping_method}_cs_full.h5ad"),          emit: cs_full_h5ad
    path("*.png"),                                              emit: cs_images
    path("*.csv"),                                              emit: cs_top_genes

    script:
    """
    echo "\n\n==================  CellSweep =================="
    echo "Meta: ${meta}"
    echo "Raw h5ad: ${raw_h5ad}"
    echo "Mapping method: ${mapping_method}"
    echo "Datatype: ${datatype}"

    run_cellsweep.py \\
        --input_h5ad ${raw_h5ad} \\
        --cs_filtered_h5ad ${meta.id}_${mapping_method}_cs_filtered.h5ad \\
        --cs_full_h5ad ${meta.id}_${mapping_method}_cs_full.h5ad \\
        --image_prefix ${meta.id}_${mapping_method}_ \\
        --expected_cells ${meta.expected_cells} \\
        --threads ${task.cpus}
    """
}
