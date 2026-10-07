process OAK_DEMUX {
    publishDir path: { "${params.outdir}/demultiplex/${meta.id}" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file)

    output:
    tuple val(meta), path("${meta.id}_S1_L001_R2_001.fastq.gz"), path("${meta.id}_S1_L001_R1_001.fastq.gz"), path("${meta.id}_S1_L001_I*_001.fastq.gz", arity: '0..2'), path(input_file), emit: demux_files
    tuple val(meta), path("${meta.id}_oak_demux.log"), emit: stats
    path "versions.yml", emit: versions

    script:
    // Index FASTQs as named by MERGE_FASTQS (<id>_merged_I1/I2.fastq.gz)
    def idx_files = [fastq_indices].flatten().findAll { it }
    def fastq_I1  = idx_files.find { it.name ==~ /.*_I1(_001)?\.f(ast)?q(\.gz)?/ }
    def fastq_I2  = idx_files.find { it.name ==~ /.*_I2(_001)?\.f(ast)?q(\.gz)?/ }
    def p5 = meta.p5?.toString()?.trim()
    def p7 = meta.p7?.toString()?.trim()
    def max_mm = params.oak_demux_max_mismatches != null ? params.oak_demux_max_mismatches : (params.seqtech_parameters[params.protocol].oak_demux_max_mismatches ?: 1)

    if (!p5) {
        error "OAK sample '${meta.id}' has no i5 barcode. Fill in the 'p5' column of the samplesheet."
    }
    if (!fastq_I2) {
        error "OAK sample '${meta.id}' has no I2 index FASTQ, which is needed to demultiplex on i5. Point 'fastq_indices' to the index reads, e.g. /path/Undetermined_S0_I*."
    }
    if (p7 && !fastq_I1) {
        error "OAK sample '${meta.id}' has an i7 barcode but no I1 index FASTQ. Point 'fastq_indices' to both index reads, e.g. /path/Undetermined_S0_I*."
    }

    // With an i7 barcode, demultiplex on i7 first and feed its output to the i5 step
    def i7_step = p7 ? """
    echo "Demultiplexing on i7 index: ${p7}"
    oakseq_custom_demux_i5_i7.sh \\
        --barcode ${p7} \\
        --index-type i7 \\
        --r1 ${fastq_BC_UMI} \\
        --r2 ${fastq_cDNA} \\
        --i1 ${fastq_I1} \\
        --i2 ${fastq_I2} \\
        --out i7/${meta.id} \\
        --max-mm ${max_mm} | tee i7_demux.out

    echo -e "after_i7\\t${p7}\\t\$(grep 'Reads matching' i7_demux.out | grep -oE '[0-9]+\$')" >> ${meta.id}_oak_demux.log

    R1_in="i7/${meta.id}_R1_001.fastq.gz"
    R2_in="i7/${meta.id}_R2_001.fastq.gz"
    I1_arg="--i1 i7/${meta.id}_I1_001.fastq.gz"
    I2_in="i7/${meta.id}_I2_001.fastq.gz"
    """ : """
    echo "No i7 index provided, skipping i7 demultiplexing."
    R1_in="${fastq_BC_UMI}"
    R2_in="${fastq_cDNA}"
    I1_arg="${fastq_I1 ? "--i1 ${fastq_I1}" : ''}"
    I2_in="${fastq_I2}"
    """
    """
    # Fail when the demultiplexing script fails, not only when tee does
    set -o pipefail

    echo -e "\\n\\n==================  Demultiplex OAK data  =================="
    echo "Processing sample: ${meta}"
    echo "FASTQ cDNA: ${fastq_cDNA}"
    echo "FASTQ BC & UMI: ${fastq_BC_UMI}"
    echo "FASTQ indices: ${idx_files.join(' ')}"
    echo "Maximum index mismatches: ${max_mm}"

    echo -e "step\\tbarcode\\treads" > ${meta.id}_oak_demux.log
    echo -e "input\\t-\\t\$(( \$(gzip -dc ${fastq_I2} | wc -l) / 4 ))" >> ${meta.id}_oak_demux.log

    ${i7_step}

    # Always demultiplex on i5; the output follows Cell Ranger's <sample>_S1_L001_<read>_001 naming
    echo "Demultiplexing on i5 index: ${p5}"
    oakseq_custom_demux_i5_i7.sh \\
        --barcode ${p5} \\
        --index-type i5 \\
        --r1 \${R1_in} \\
        --r2 \${R2_in} \\
        \${I1_arg} \\
        --i2 \${I2_in} \\
        --out ${meta.id}_S1_L001 \\
        --max-mm ${max_mm} | tee i5_demux.out

    echo -e "after_i5\\t${p5}\\t\$(grep 'Reads matching' i5_demux.out | grep -oE '[0-9]+\$')" >> ${meta.id}_oak_demux.log

    # Intermediate i7 reads are no longer needed
    rm -rf i7

    echo "Demultiplexing completed."
    cat ${meta.id}_oak_demux.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1 | grep Version | sed 's/Version: //'))
    END_VERSIONS
    """
}
