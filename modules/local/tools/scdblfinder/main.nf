process SCDBLFINDER {
    publishDir "${params.outdir}/doublet_filtering/${meta.id}/scdblfinder", mode: 'copy'
    tag "${meta.id} | ${mapping_method}"
    label 'process_medium'

    container { demuxafy_sif }

    input:
    tuple val(meta), path(tenx_dir), val(mapping_method), val(datatype), path(filtered_barcodes)
    val demuxafy_sif

    output:
    tuple val(meta), path("${meta.id}_scDblFinder_results.tsv"), val(mapping_method), val(datatype), emit: scdblfinder_results
    path("${meta.id}_scDblFinder_summary.tsv"),                                                      emit: scdblfinder_summary, optional: true

    script:
    // STARsolo's raw output includes near-empty droplets, which inflates the cell count enough
    // that scDblFinder's automatic doublet-rate estimate (dbr) falls outside its required (0,1)
    // range. Restrict doublet calling to real cells when a filtered list is given.
    def filter_arg = filtered_barcodes.name != 'NO_FILE' ? "-b ${filtered_barcodes}" : ""
    """
    echo "\n\n==================  scDblFinder =================="
    echo "Meta: ${meta}"
    echo "10x dir: ${tenx_dir}"
    echo "Filtered barcodes: ${filtered_barcodes}"

    # scDblFinder.R doesn't create its own output directory before writing into it.
    mkdir -p scdblfinder_out

    # scDblFinder.R -o scdblfinder_out -t ${tenx_dir} ${filter_arg}
    scDblFinder.R -o scdblfinder_out -t ${tenx_dir} ""

    mv scdblfinder_out/scDblFinder_doublets_singlets.tsv ${meta.id}_scDblFinder_results.tsv
    if [ -f scdblfinder_out/scDblFinder_doublet_summary.tsv ]; then
        mv scdblfinder_out/scDblFinder_doublet_summary.tsv ${meta.id}_scDblFinder_summary.tsv
    fi
    """
}
