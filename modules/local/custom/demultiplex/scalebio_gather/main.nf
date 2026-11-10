process SCALEBIO_GATHER {
    publishDir path: { "${params.outdir}/demultiplex/scalebio/${meta.id}" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python:3.14.2--0bd36b5fd9edb930' :
        'community.wave.seqera.io/library/python:3.14.2--0562a0df3245213a' }"

    input:
    // A sample's reads from every chunk (lane) of every library it was in, and those chunks' metrics
    tuple val(meta), path(cdna_parts, stageAs: 'parts/cdna/*'), path(bc_parts, stageAs: 'parts/bc/*'), path(metrics, stageAs: 'metrics/*'), path(input_file)

    output:
    tuple val(meta), path("${meta.id}_scalebio_cDNA.fastq.gz"), path("${meta.id}_scalebio_BC_UMI.fastq.gz"), path(input_file), emit: reads
    tuple val(meta), path("${meta.id}.scalebio_demux.json"), emit: stats
    path "versions.yml", emit: versions

    script:
    """
    echo -e "\\n\\n==================  ScaleBio: gather ${meta.id}  =================="

    # The parts are named <sample>.<library>_<lane>_{cDNA,BC_UMI}.fastq.gz, so sorting by
    # name puts both reads' parts in the same order; a run of gzip members is one gzip file.
    gather() {
        local out=\$1 dir=\$2
        local parts=( \$(ls -1 "\$dir" | LC_ALL=C sort) )
        if [ "\${#parts[@]}" -eq 1 ]; then
            ln -s "\$dir/\${parts[0]}" "\$out"
        else
            for p in "\${parts[@]}"; do cat "\$dir/\$p"; done > "\$out"
        fi
    }
    gather "${meta.id}_scalebio_cDNA.fastq.gz" parts/cdna
    gather "${meta.id}_scalebio_BC_UMI.fastq.gz" parts/bc

    scalebio_demux.py gather \\
        --sample ${meta.id} \\
        --metrics metrics/* \\
        --out ${meta.id}.scalebio_demux.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
