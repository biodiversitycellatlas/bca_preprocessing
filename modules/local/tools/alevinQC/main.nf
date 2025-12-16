process ALEVIN_QC {
    publishDir "${params.outdir}/mapping_alevin/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/bioconductor-alevinqc_r-base:62474ff3de946976"

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file)
    path(alevin_fry_output)

    output:
    path("${meta.id}_run/${meta.id}_alevinFry_QC.html"), emit: alevinQC_report
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file), emit: alevinQC_samplesheet

    script:
    """
    echo "\n\n==================  ALEVIN-FRY =================="
    echo "Sample ID: ${meta.id}"
    echo "Alevin-Fry output: ${alevin_fry_output}"

    alevin_QC.R ${alevin_fry_output} ${meta.id}
    """
}
