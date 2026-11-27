process SPLITSEQ_WELLS {
    publishDir path: { "${params.outdir}/demultiplex/splitseq_whitelists" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }

    // Reads a few hundred barcodes, once per run
    cpus   1
    memory 1.GB
    time   15.m

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python:3.14.2--0bd36b5fd9edb930' :
        'community.wave.seqera.io/library/python:3.14.2--0562a0df3245213a' }"

    input:
    // 'bc_whitelist': one list for all three rounds, or three lists in the order bc1 (RT),
    // bc2, bc3. Each segment in a folder of its own, so a list given twice does not clash.
    path(whitelists, stageAs: 'segment?/*')
    path(wells, stageAs: 'given_wells/*')       // [] unless 'splitseq_wells' is set

    output:
    path "splitseq_bc1.txt",   emit: bc1
    path "splitseq_bc2.txt",   emit: bc2
    path "splitseq_bc3.txt",   emit: bc3
    path "splitseq_wells.csv", emit: wells
    path "versions.yml",       emit: versions

    script:
    def wells_arg = wells ? "--wells ${wells}" : ''
    """
    echo -e "\\n\\n==================  SPLiT-seq: barcode lists and round-1 wells  =================="

    splitseq_wells.py ${whitelists} ${wells_arg} --output .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
