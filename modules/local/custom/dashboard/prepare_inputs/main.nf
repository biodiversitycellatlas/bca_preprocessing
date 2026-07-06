process PREPARE_DASHBOARD_INPUTS {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(files)

    output:
    path "*_Summary.csv"           , emit: summary, optional: true
    path "*_CellReads.stats"       , emit: cell_stats, optional: true
    path "*_UMIperCellSorted.txt"  , emit: knee_files, optional: true
    path "*_meta_info.json"         , emit: af_meta_info, optional: true
    path "*_quant.json"            , emit: af_quant_json, optional: true
    path "*_cell_meta.tsv"         , emit: af_cell_meta, optional: true
    path "versions.yml"            , emit: versions

    script:
    """
    for f in $files; do
        # Use simple bash logic to match the file type and link it
        if [[ "\$f" == *"Summary.csv"* ]]; then ln -s "\$f" "${meta.id}_Summary.csv"; fi
        if [[ "\$f" == *"CellReads.stats"* ]]; then ln -s "\$f" "${meta.id}_CellReads.stats"; fi
        if [[ "\$f" == *"UMIperCellSorted.txt"* ]]; then ln -s "\$f" "${meta.id}_UMIperCellSorted.txt"; fi
        if [[ "\$f" == *"meta_info.json"* ]]; then ln -s "\$f" "${meta.id}_meta_info.json"; fi
        if [[ "\$f" == *"quant.json"* ]]; then ln -s "\$f" "${meta.id}_quant.json"; fi
        if [[ "\$f" == *"cell_meta.tsv"* ]]; then ln -s "\$f" "${meta.id}_cell_meta.tsv"; fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | sed 's/^GNU bash, version //; s/ .*\$//')
    END_VERSIONS
    """
}
