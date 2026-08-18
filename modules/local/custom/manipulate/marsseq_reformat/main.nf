process MARSSEQ_BUILD_READS {
    publishDir "${params.outdir}/marsseq_reformat/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/python:3.14.2--0bd36b5fd9edb930"

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file)

    output:
    tuple val(meta), path("${meta.id}_marsseq_cDNA.fastq.gz"), path("${meta.id}_marsseq_BC_UMI.fastq.gz"), path(fastq_indices), path(input_file), emit: reformatted_files
    tuple val(meta), path("${meta.id}_marsseq_reformat.log"), emit: stats
    path "versions.yml", emit: versions

    script:
    def marsseq_read1_design = params.marsseq_read1_design ?: params.seqtech_parameters[params.protocol].marsseq_read1_design
    def marsseq_read2_design = params.marsseq_read2_design ?: params.seqtech_parameters[params.protocol].marsseq_read2_design

    if (!marsseq_read1_design || !marsseq_read2_design) {
        error "No MARS-seq read design defined for protocol '${params.protocol}'. Set 'marsseq_read1_design' and 'marsseq_read2_design' in the configuration file."
    }
    """
    echo "\n\n==================  MARS-seq: Rebuilding reads  =================="
    echo "Processing sample: ${meta}"
    echo "Fastq files: ${fastq_cDNA}, ${fastq_BC_UMI}"
    echo "Read 1 design: ${marsseq_read1_design}"
    echo "Read 2 design: ${marsseq_read2_design}"

    # Strip the ignore bases from read 1 and prepend its batch barcode to the read 2 barcode & UMI
    marsseq_build_reads.py \\
        --sample_id ${meta.id} \\
        --fq1 ${fastq_cDNA} \\
        --fq2 ${fastq_BC_UMI} \\
        --read1-design ${marsseq_read1_design} \\
        --read2-design ${marsseq_read2_design} \\
        --output .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
