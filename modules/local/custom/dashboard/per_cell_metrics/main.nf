process PERCELL_METRICS {
    publishDir "${params.outdir}/summary_results/per-cell_metrics", mode: 'copy'
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(star_bam)
    tuple val(meta), path(star_solodir)
    tuple val(meta), path(star_log)
    tuple val(meta), path(secondderiv_cutoff)
    path (ref_gtf)

    output:
    path("*_metrics.json"), emit: percell_json
    path("*.png"),  emit: percell_imgs
    path "versions.yml",    emit: versions

    script:
    def cfg_name = "GeneFull_Ex50pAS"
    def cutoff_file = secondderiv_cutoff ?: ''
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

    per-cell_images.py \\
        --solo-output ${star_solodir} \\
        --bam ${star_bam} \\
        --mt-contig ${params.mt_contig} \\
        --gtf ${ref_gtf} \\
        --outdir . \\
        --min-reads \${cell_thres}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
