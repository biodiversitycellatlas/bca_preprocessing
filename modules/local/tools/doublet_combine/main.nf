process COMBINE_DOUBLET_RESULTS {
    publishDir path: { "${params.outdir}/doublet_filtering/${meta.id}/combined" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${meta.id} | ${meta.mapping_method}"
    label 'process_single2'
    // As the two callers: a failure leaves the sample without doublet annotation
    // instead of ending the run
    label 'error_optional'

    // Runs in scDblFinder's environment, which has R, argparse and data.table
    conda "${moduleDir}/../scdblfinder/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-base_bioconductor-scdblfinder_r-xgboost_r-argparse_pruned:60e56785b2c0e166' :
        'community.wave.seqera.io/library/r-base_bioconductor-scdblfinder_r-xgboost_r-argparse_pruned:1a632f336ca0f064' }"

    input:
    // Either results file may be empty ([]) when that caller could not run on the
    // matrix; the consensus then holds the other caller's calls alone
    tuple val(meta), path(scrublet_results), path(scdblfinder_results)

    output:
    tuple val(meta), path("${meta.id}_combined_doublets_w_combined_assignments.tsv"), emit: combined_results
    path("${meta.id}_combined_doublets_summary.tsv"),                                 emit: combined_summary, optional: true

    script:
    def scrublet_arg    = scrublet_results    ? "-r ${scrublet_results}"    : ''
    def scdblfinder_arg = scdblfinder_results ? "-n ${scdblfinder_results}" : ''
    """
    echo -e "\\n\\n==================  Combine doublet results =================="
    echo "Meta: ${meta}"
    echo "Consensus method: ${params.doublet_consensus_method}"

    demuxafy_combine_results.R \\
        -o ${meta.id}_combined_doublets.tsv \\
        ${scrublet_arg} \\
        ${scdblfinder_arg} \\
        -m "${params.doublet_consensus_method}"
    """
}
