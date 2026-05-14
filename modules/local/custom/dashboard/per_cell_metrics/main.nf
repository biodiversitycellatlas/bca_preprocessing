process PERCELL_METRICS {
    publishDir "${params.outdir}/summary_results/per-cell_metrics", mode: 'copy'
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(star_bam)
    tuple val(meta), path(star_solodir)
    tuple val(meta), path(star_log)
    path (ref_gtf)

    output:
    path("*_metrics.json"), emit: percell_json
    path("*.png"),  emit: percell_imgs

    script:
    def cfg_name = "GeneFull_Ex50pAS"
    """
    # Extract the UMI threshold using awk
    starsolo_thres=\$(awk -v cfg="${cfg_name}" '
        /Starting Solo post-map for/ { if (\$0 ~ cfg) in_block=1; else in_block=0 }
        in_block && /cellFiltering/ {
            if (match(\$0, /nUMImin=([0-9]+)/, m)) print m[1]
        }
    ' ${star_log})

    # Fallback to 0 if threshold wasn't found
    if [ -z "\$starsolo_thres" ]; then
        starsolo_thres=0
    fi

    per-cell_images.py \\
        --solo-output ${star_solodir} \\
        --bam ${star_bam} \\
        --mt-contig ${params.mt_contig} \\
        --gtf ${ref_gtf} \\
        --outdir . \\
        --min-reads \${starsolo_thres}
    """
}
