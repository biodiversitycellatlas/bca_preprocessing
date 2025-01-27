// ==================  MAPPING STATISTICS  =================== \\ 

process MAPPING_STATS {
    input:
    val(trigger)
    
    output:
    file("${params.resDir}/mapping_stats.tsv")

    script:
    """
    echo "\n\n==================  MAPPING STATISTICS  =================="

    sbatch ${params.baseDir}/scripts/mapping_statistics.sh ${params.resDir}
    """
}
