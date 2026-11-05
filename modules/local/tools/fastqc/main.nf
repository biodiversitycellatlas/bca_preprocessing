process FASTQC {
    publishDir path: { "${params.outdir}/fastqc" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${fastq_cDNA}, ${fastq_BC_UMI}"

    // Time follows the FASTQs file size
    cpus   { 1 * task.attempt }
    memory { 4.GB * task.attempt }
    time   { BcaResources.scaledTime(
        params.dynamic_time?.FASTQC, [fastq_cDNA, fastq_BC_UMI], task.attempt, 2) }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/fastqc:0.12.1--104d26ddd9519960' :
        'community.wave.seqera.io/library/fastqc:0.12.1--9971ea336a9eddae' }"

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file)

    output:
    path "*_fastqc.{zip,html}", emit: fastqc_results
    path "versions.yml",        emit: versions

    script:
    """
    fastqc ${fastq_cDNA}
    fastqc ${fastq_BC_UMI}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/FastQC v//')
    END_VERSIONS
    """
}
