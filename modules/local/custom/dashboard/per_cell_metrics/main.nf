process PERCELL_METRICS {
    publishDir "${params.outdir}/summary_results/per-cell_metrics", mode: 'copy'
    tag "${meta.id}"
    label 'process_single_long2'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta),   path(star_bam)
    tuple val(_meta2), path(star_solodir)
    tuple val(_meta3), path(star_log)
    tuple val(_meta4), path(secondderiv_cutoff)
    tuple val(_meta5), path(filtered_matrix_dir)
    path (ref_gtf)

    output:
    path("*_metrics.json"), emit: percell_json
    path("*.png"),  emit: percell_imgs
    path "versions.yml",    emit: versions

    script:
    def cfg_name = "GeneFull_Ex50pAS"
    def cutoff_file = secondderiv_cutoff ?: ''
    def cell_barcodes = filtered_matrix_dir ? "${filtered_matrix_dir}/barcodes.tsv" : ''
    def mt_contigs = (params.mt_contig ?: '').toString().trim()
    if (!mt_contigs) {
        error "PERCELL_METRICS requires params.mt_contig to name at least one mitochondrial contig"
    }
    """
    # The second-derivative cutoff is the one the matrices were filtered on, so it takes
    # precedence over STARsolo's nUMImin whenever that method produced a cutoff file.
    cell_thres=""
    if [ -n "${cutoff_file}" ] && [ -s "${cutoff_file}" ]; then
        cell_thres=\$(cat ${cutoff_file})
        echo "Using second-derivative UMI threshold: \$cell_thres"
    else
        # Extract the UMI threshold using awk
        cell_thres=\$(awk -v cfg="${cfg_name}" '
            /Starting Solo post-map for/ { if (\$0 ~ cfg) in_block=1; else in_block=0 }
            in_block && /cellFiltering/ {
                if (match(\$0, /nUMImin=([0-9]+)/, m)) print m[1]
            }
        ' ${star_log})
        echo "Using STARsolo nUMImin threshold: \$cell_thres"
    fi

    # Fallback to 0 if threshold wasn't found
    if [ -z "\$cell_thres" ]; then
        cell_thres=0
    fi

    # The barcode list is exact; the threshold is only a fallback, since it compares an UMI cutoff against per-barcode read counts
    cell_bc_arg=""
    if [ -s "${cell_barcodes}" ]; then
        cell_bc_arg="--cell-barcodes ${cell_barcodes}"
    else
        echo "No filtered barcodes.tsv staged; splitting cells on \$cell_thres reads"
    fi

    per-cell_images.py \\
        --solo-output ${star_solodir} \\
        --bam ${star_bam} \\
        --mt-contig ${mt_contigs} \\
        --gtf ${ref_gtf} \\
        --outdir . \\
        --min-reads \${cell_thres} \\
        \${cell_bc_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
