process SEQSPEC_LAYOUT {
    tag "${seqspec.name}"
    publishDir path: { "${params.outdir}/seqspec" }, mode: 'copy', pattern: '{seqspec_layout.env,onlists.tsv}'

    // Reads one YAML file; downloading and decompressing a remote onlist
    // is what takes the time and the memory.
    cpus   1
    memory { 4.GB * task.attempt }
    time   { 30.m * task.attempt }

    // The environment and container of STARSOLO_ALIGN, which carry seqspec
    conda "${moduleDir}/../../../tools/star/starsolo_align/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/htslib_samtools_seqspec_star_pruned:cef769c7e3b03dd0' :
        'community.wave.seqera.io/library/htslib_samtools_seqspec_star_pruned:583718e0c3d7bddd' }"

    input:
    path(seqspec, stageAs: 'seqspec/*')
    val fetch_onlists                            // download the remote onlists (bc_whitelist unset)

    output:
    path "seqspec_layout.env",                   emit: layout
    path "onlists.tsv",                          emit: onlists
    path "onlist_*.txt",                         emit: onlist_files, optional: true
    path "versions.yml",                         emit: versions

    script:
    """
    echo -e "\\n\\n==============  SEQSPEC LAYOUT  ================"
    echo "seqspec file: ${seqspec}"

    seqspec_layout.py \\
        --spec ${seqspec} \\
        --modality rna \\
        ${fetch_onlists ? '--fetch-onlists' : ''} \\
        --out-env seqspec_layout.env \\
        --out-onlists onlists.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqspec: \$(seqspec --version 2>&1 | sed 's/^seqspec //')
    END_VERSIONS
    """
}
