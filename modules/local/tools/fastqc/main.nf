process FASTQC {
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    tag "${fastq_cDNA}, ${fastq_BC_UMI}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/fastqc:0.12.1--104d26ddd9519960"

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
