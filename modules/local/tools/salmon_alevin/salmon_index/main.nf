process SALMON_INDEX {
    publishDir path: { "${params.outdir}/genome" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }

    // Memory tracks the size of the splici transcriptome being indexed. Coefficients live in
    // params.dynamic_memory; without an entry the request is process_high's 64 GB.
    cpus   8
    time   { BcaResources.attemptTime(task, '2h') }
    memory { BcaResources.attemptMemory(task, BcaResources.scaledMemory(
        params.dynamic_memory?.SALMON_INDEX, [splici_index_reference], 1, 64)) }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/salmon:1.10.3--726401738a281398' :
        'community.wave.seqera.io/library/salmon:1.10.3--6ea6d56a4b0c26b8' }"

    input:
    path(splici_index_reference)
    val(index_suffix)
    val(kmer)                   // k-mer length (odd, at most 31): mapping_alevin.nf resolves it from salmon_index_k

    output:
    path("salmon_index${index_suffix}/"), emit: salmon_index
    path "versions.yml",   emit: versions


    script:
    """
    echo -e "\\n\\n==================  SALMON INDEX =================="
    echo "k-mer length: ${kmer}"
    # Build reference index
    salmon index \\
        -t ${splici_index_reference}/*.fa \\
        -i ./salmon_index${index_suffix} \\
        -k ${kmer} \\
        -p ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(salmon --version | sed 's/salmon //')
    END_VERSIONS
    """
}
