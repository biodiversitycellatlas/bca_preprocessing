process METACELL2 {
    publishDir path: { "${params.outdir}/metacells/${meta.id}" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${meta.id} | ${meta.mapping_method}"

    // A sample with too few cells is not an error: run_metacells.py writes a 'skipped'
    // summary and exits 0, so the task is cached like any other. Only real failures
    // land here, and those drop the sample from the report rather than end the run.
    errorStrategy {
        if (task.exitStatus in ((130..145) + 104 + 175 + 255)) {
            return 'retry'
        }
        log.warn "No metacells for ${meta.id} (${meta.mapping_method}): Metacell2 failed with exit status " +
            "${task.exitStatus} (see .command.log in the task work directory). The sample is left out of " +
            "the filtering report."
        return 'ignore'
    }
    maxRetries 1

    // Metacell2 splits its work by the processor count, and the grouping can depend on
    // how it was split: the count is fixed rather than escalated, so a retry or a rerun
    // on another node gives the assignment a selection was made against. Deliberately
    // no resource label, since a withLabel setting would take precedence over these.
    cpus 4
    time { 8.h * task.attempt }

    memory { BcaResources.scaledMemory(
        params.dynamic_memory?.METACELL2, [cells_h5ad], task.attempt, 36) }

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(cells_h5ad), path(ambient_h5ad)
    path(gene_table)

    output:
    tuple val(meta), path("${meta.id}_mc2_cells.h5ad"),     emit: cells_h5ad,     optional: true
    tuple val(meta), path("${meta.id}_mc2_metacells.h5ad"), emit: metacells_h5ad, optional: true
    tuple val(meta), path("${meta.id}_mc_summary.json"),    emit: summary
    path "versions.yml",                                    emit: versions

    script:
    def ambient_arg  = ambient_h5ad ? "--ambient_h5ad ${ambient_h5ad}" : ""
    def max_umis_arg = params.mc2_max_cell_umis ? "--max_cell_umis ${params.mc2_max_cell_umis}" : ""
    def max_excl_arg = params.mc2_max_excluded_genes_fraction != null
        ? "--max_excluded_genes_fraction ${params.mc2_max_excluded_genes_fraction}" : ""
    def removed_arg  = params.perform_doublet_filtering ? "--doublets_removed" : ""
    """
    export METACELLS_PROCESSORS_COUNT=${task.cpus}
    export NUMBA_CACHE_DIR=\$PWD/.numba_cache

    run_metacells.py \\
        --cells_h5ad ${cells_h5ad} \\
        ${ambient_arg} \\
        --gene_table ${gene_table} \\
        --sample_id ${meta.id} \\
        --mapping_method ${meta.mapping_method} \\
        --prefix ${meta.id} \\
        --target_metacell_size ${params.mc2_target_metacell_size} \\
        --random_seed ${params.mc2_random_seed} \\
        --min_cell_umis ${params.mc2_min_cell_umis} \\
        ${max_umis_arg} \\
        ${max_excl_arg} \\
        --excluded_gene_patterns '${params.mc2_excluded_gene_patterns ?: ''}' \\
        --lateral_gene_patterns '${params.mc2_lateral_gene_patterns ?: ''}' \\
        --min_cells ${params.mc2_min_cells} \\
        ${removed_arg} \\
        --cpus ${task.cpus} \\
        --compression ${params.h5ad_compression ?: 'gzip'} \\
        --versions_yml versions.yml \\
        --process_name "${task.process}"
    """

    stub:
    """
    touch ${meta.id}_mc2_cells.h5ad ${meta.id}_mc2_metacells.h5ad versions.yml
    echo '{"schema":"bca_metacell_summary","schema_version":1,"sample":{"id":"${meta.id}","status":"skipped","reason":"stub run"},"metacells":[],"genes":null}' > ${meta.id}_mc_summary.json
    """
}
