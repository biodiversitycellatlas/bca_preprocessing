process ALEVIN_QC {
    publishDir path: { "${params.outdir}/mapping_alevin/${meta.id}" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${meta.id}"
    // alevinQC runs on one thread, but its memory follows the permit list (over 12 GB for
    // a 6.8 M-barcode one). No resource label, since a withLabel setting would take
    // precedence over cpus: memory and time are process_medium's.
    cpus   1
    memory { 36.GB * task.attempt }
    time   { 6.h   * task.attempt }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-alevinqc_r-base:62474ff3de946976' :
        'community.wave.seqera.io/library/bioconductor-alevinqc_r-base:6d2d752a32e8c96a' }"

    input:
    tuple val(meta), path(alevin_fry_output)

    output:
    tuple val(meta), path("${meta.id}_run/${meta.id}_alevinFry_QC.html"), emit: alevinQC_report
    path "versions.yml", emit: versions

    script:
    """
    echo -e "\\n\\n==================  ALEVIN-FRY =================="
    echo "Sample ID: ${meta.id}"
    echo "Alevin-Fry output: ${alevin_fry_output}"

    alevin_QC.R ${alevin_fry_output} ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript --version 2>&1 | sed 's/^.*version //; s/ .*\$//')
    END_VERSIONS
    """
}
