process PREPARE_ADDFEATURE {
    publishDir path: { "${params.outdir}/genome" }, mode: 'copy', saveAs: { filename -> filename in ['addfeature_report.tsv', 'rrna_regions.tsv'] ? filename : null }
    label 'process_single_long2'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python:3.14.2--0bd36b5fd9edb930' :
        'community.wave.seqera.io/library/python:3.14.2--0562a0df3245213a' }"
    conda "${moduleDir}/environment.yml"

    input:
    path base_fasta
    path add_fasta
    path add_gtf

    output:
    path "addfeature.fasta",      emit: fasta
    path "rrna_regions.tsv",      emit: rrna_regions
    path "addfeature_report.tsv", emit: report
    path "versions.yml",          emit: versions

    script:
    // Either added file may be absent: the GTF alone gives rRNA loci on the genome,
    // the FASTA alone extra contigs that count as rRNA over their full length
    def fasta_arg = add_fasta ? "--add-fasta ${add_fasta}" : ""
    def gtf_arg   = add_gtf ? "--add-gtf ${add_gtf}" : ""
    """
    prepare_addfeature.py \\
        --base-fasta ${base_fasta} \\
        ${fasta_arg} \\
        ${gtf_arg} \\
        --out-fasta addfeature.fasta \\
        --out-regions rrna_regions.tsv \\
        --out-report addfeature_report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch addfeature.fasta rrna_regions.tsv addfeature_report.tsv versions.yml
    """
}
