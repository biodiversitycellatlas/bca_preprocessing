process PREPARE_DASHBOARD_INPUTS {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(summary), path(cell_stats), path(knee_files)

    output:
    path "*_Summary.csv"           , emit: summary
    path "*_CellReads.stats"       , emit: cell_stats
    path "*_UMIperCellSorted.txt"  , emit: knee_files

    script:
    """
    # Rename files by adding the Sample ID prefix
    ln -s ${summary}    ${meta.id}_Summary.csv
    ln -s ${cell_stats} ${meta.id}_CellReads.stats
    ln -s ${knee_files} ${meta.id}_UMIperCellSorted.txt
    """
}
