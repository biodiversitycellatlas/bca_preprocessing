process MERGE_REF_FASTA {
    input:
        path base_fasta
        path add_fasta

    output:
        path "ref.fasta"

    script:
    def do_merge = add_fasta ? true : false
    """
    if [ "$do_merge" = "true" ]; then
        cat ${base_fasta} ${add_fasta} > ref.fasta
    else
        cp ${base_fasta} ref.fasta
    fi
    """
}
