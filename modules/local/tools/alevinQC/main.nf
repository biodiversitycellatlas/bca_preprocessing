process ALEVIN_QC {
    publishDir "${params.outdir}/mapping_alevin/${meta.id}", mode: 'copy'
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file)
    path(alevin_fry_output)

    output:
    path("*_alevinFry_QC.html"), emit: alevinQC_report

    script:
    """
    echo "\n\n==================  ALEVIN-FRY =================="
    echo "Sample ID: ${meta.id}"
    echo "Alevin-Fry output: ${alevin_fry_output}"

    Rscript alevin_QC.R ${alevin_fry_output} ${meta.id}
    """
}
