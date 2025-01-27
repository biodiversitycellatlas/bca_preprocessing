// ==================  MAPPING STATISTICS  =================== \\ 

process MAPPING_STATS {
    publishDir "${params.resDir}/mapping_stats", mode: 'copy'
    debug true

    input:
    val(trigger)
    
    output:
    file("mapping_stats.tsv")

    script:
    """
    echo "\n\n==================  MAPPING STATISTICS  =================="

    sbatch ${params.baseDir}/scripts/mapping_statistics.sh ${params.dataDir}
    """
}
