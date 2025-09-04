process MERGE_FASTQS {
    tag { meta.id }
    label 'process_single'

    input:
    val  meta
    path cDNA_list
    path BCUMI_list

    output:
    tuple val(meta), path("${meta.id}_cDNA.fastq.gz"), path("${meta.id}_BC_UMI.fastq.gz")

    script:
    """
    # cDNA
    cat ${cDNA_list.join(' ')} > ${meta.id}_cDNA.fastq.gz
    # Barcodes + UMIs
    cat ${BCUMI_list.join(' ')} > ${meta.id}_BC_UMI.fastq.gz
    """
}
