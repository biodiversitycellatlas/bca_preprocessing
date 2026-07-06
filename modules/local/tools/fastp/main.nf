process FASTP {
    publishDir "${params.outdir}/fastp", mode: 'copy'
    tag "${meta.id}"
    label 'process_low'
    debug true

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/fastp:0.24.0--0397de619771c7ae"

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file)

    output:
    tuple val(meta), path("trimmed_${fastq_cDNA}"), path("trimmed_${fastq_BC_UMI}"), path(fastq_indices), path(input_file), emit: trimmed_files
    path "versions.yml", emit: versions

    script:
    // Retrieve fastp settings from custom parameters if set, otherwise from conf/seqtech_parameters.config
    def fastp_length_required = params.fastp_length_required ?: params.seqtech_parameters[params.protocol].fastp_length_required
    def fastp_qualified_quality_phred = params.fastp_qualified_quality_phred ?: params.seqtech_parameters[params.protocol].fastp_qualified_quality_phred

    """
    echo "\n\n==================  TRIM FASTQs WITH FASTP  =================="
    echo "Metadata: ${meta}"
    echo "FASTQ cDNA: ${fastq_cDNA}"
    echo "FASTQ BC & UMI: ${fastq_BC_UMI}"
    echo "-- Settings:"
    echo "-- fastp nextflow.config: qualified_quality_phred: ${fastp_qualified_quality_phred}"
    echo "-- fastp nextflow.config: length_required: ${fastp_length_required}"

    fastp \\
        --html fastp_${meta.id}.html \\
        --json fastp_${meta.id}.json \\
        --thread 8 \\
        --in1 ${fastq_cDNA} \\
        --in2 ${fastq_BC_UMI} \\
        --out1 trimmed_${fastq_cDNA} \\
        --out2 trimmed_${fastq_BC_UMI} \\
        --length_required ${fastp_length_required} \\
        --qualified_quality_phred ${fastp_qualified_quality_phred} \\
        --detect_adapter_for_pe \\
        --dont_eval_duplication

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed 's/fastp //')
    END_VERSIONS
    """
}
