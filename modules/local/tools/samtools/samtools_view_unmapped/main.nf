process SAMTOOLS_VIEW_UNMAPPED {
    tag "${meta.id}"
    label 'process_low2'

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/samtools:1.22.1--9a10f06c24cdf05f"

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("${meta.id}_unmapped.fasta"),     emit : filtered_unmapped_fasta
    path "versions.yml",                                    emit : versions

    script:
    """
    echo "\n\n==================  SAMTOOLS VIEW UNMAPPED =================="
    echo "Metadata: ${meta}"
    echo "BAM file: ${bam_file}"

    # Extract unmapped reads (include only unmapped with -f 4)
    samtools view --threads ${task.cpus} -f 0x4 -b ${bam_file} | samtools fasta - > ${meta.id}_unmapped.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed '1!d; s/samtools //')
    END_VERSIONS
    """
}
