process SCDBLFINDER {
    publishDir "${params.outdir}/doublet_filtering/${meta.id}/scdblfinder", mode: 'copy'
    tag "${meta.id} | ${meta.mapping_method}"
    label 'process_medium'

    container { demuxafy_sif }

    input:
    tuple val(meta), path(tenx_dir)
    val demuxafy_sif

    output:
    tuple val(meta), path("${meta.id}_scDblFinder_results.tsv"), emit: scdblfinder_results
    path("${meta.id}_scDblFinder_summary.tsv"),                  emit: scdblfinder_summary, optional: true

    script:
    """
    echo "\n\n==================  scDblFinder =================="
    echo "Meta: ${meta}"
    echo "10x dir: ${tenx_dir}"
    echo "Datatype: ${meta.datatype}"

    # scDblFinder.R doesn't create its own output directory before writing into it.
    mkdir -p scdblfinder_out

    scDblFinder.R -o scdblfinder_out -t ${tenx_dir} ""

    mv scdblfinder_out/scDblFinder_doublets_singlets.tsv ${meta.id}_scDblFinder_results.tsv
    if [ -f scdblfinder_out/scDblFinder_doublet_summary.tsv ]; then
        mv scdblfinder_out/scDblFinder_doublet_summary.tsv ${meta.id}_scDblFinder_summary.tsv
    fi
    """
}
