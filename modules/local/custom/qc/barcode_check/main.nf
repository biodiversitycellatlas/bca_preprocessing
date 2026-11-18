process BARCODE_CHECK {
    tag "${meta.id}"

    // Reads the start of one FASTQ and the whitelist: seconds, and the largest whitelist fits well within this.
    cpus   1
    memory { 4.GB * task.attempt }
    time   { 30.m * task.attempt }

    // The environment and container of STARSOLO_ALIGN: pigz and awk for the check, seqspec
    // for the layout of protocol 'seqspec'
    conda "${moduleDir}/../../../tools/star/starsolo_align/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/htslib_samtools_seqspec_star_pruned:cef769c7e3b03dd0' :
        'community.wave.seqera.io/library/htslib_samtools_seqspec_star_pruned:583718e0c3d7bddd' }"

    input:
    tuple val(meta), path(fastq_BC_UMI)
    path(bc_whitelist, stageAs: 'whitelist_?/*')  // in barcode-segment order; empty without a whitelist
    val solo_type_string                          // the --soloType ... options STARsolo gets for this protocol
    path(seqspec, stageAs: 'seqspec/*')           // the seqspec file for protocol 'seqspec', empty otherwise
    val min_barcode_match                         // % of reads on the whitelist below which the sample fails
    val check_reads                               // number of reads checked, from the start of the file

    output:
    tuple val(meta), path("${meta.id}.barcode_check.tsv"), emit: verdict
    path "versions.yml",                                   emit: versions

    script:
    def whitelists = (bc_whitelist instanceof List ? bc_whitelist : [bc_whitelist]).findAll { wl -> wl }
    """
    echo -e "\\n\\n==============  BARCODE CHECK  ================"
    echo "Sample: ${meta.id}"
    echo "FASTQ BC & UMI: ${fastq_BC_UMI}"
    echo "Barcode whitelist: ${whitelists ? whitelists.join(' ') : 'None'}"

    if [[ -n "${seqspec}" ]]; then
        # Barcode structure information from seqspec file, as STARSOLO_ALIGN reads it
        layout="--soloType CB_UMI_Simple \$(seqspec index -m rna -t starsolo -s file ${seqspec})"
    else
        layout="${solo_type_string}"
    fi
    echo "Layout: \${layout}"

    # The verdict (pass, fail or skipped) is in the TSV: a failed sample is stopped or
    # dropped by the workflow (lib/BarcodeCheck.groovy), not by this task's exit code.
    barcode_check.sh \\
        --sample ${meta.id} \\
        --fastq ${fastq_BC_UMI} \\
        --layout "\${layout}" \\
        --min ${min_barcode_match} \\
        --reads ${check_reads} \\
        --out ${meta.id}.barcode_check.tsv \\
        ${whitelists.join(' ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk --version 2>/dev/null | head -n1 | sed 's/^GNU Awk //; s/,.*\$//')
    END_VERSIONS
    """
}
