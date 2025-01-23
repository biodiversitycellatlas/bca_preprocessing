// ==================  MAPPING STATISTICS  =================== \\ 

process MAPPING_STATS {
    publishDir "${params.resDir}/mapping_stats", mode: 'copy'
    debug true

    input:
    // tuple val(sample_id_N), val(config_name_N), path(files_N)
    // tuple val(sample_id_CR), val(config_name_CR), path(files_CR)

    output:
    file("mapping_stats.tsv")

    script:
    """
    sbatch ${params.baseDir}/scripts/mapping_statistics.sh ${params.dataDir}
    """
}
