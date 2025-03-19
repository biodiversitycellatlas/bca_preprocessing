process MAPPING_STATS {
    input:
    val(trigger)

    script:
    """
    sbatch ${params.baseDir}/bin/mapping_statistics.sh ${params.resDir}
    """
}
