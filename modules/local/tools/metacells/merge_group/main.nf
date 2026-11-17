process MERGE_GROUP_H5AD {
    publishDir path: { "${params.outdir}/metacells/${meta.id}" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${meta.id} | ${meta.mapping_method} | ${meta.members.size()} samples"

    // Memory tracks the members' cell-called objects, which are all held to concatenate
    // them. Of each raw (or full) h5ad only obs, var and the called rows of CellSweep's
    // layer are read, so they raise the floor (resident_factor in the entry) rather than
    // adding to the input, as in METACELL2.
    cpus   1
    time   { 2.h * task.attempt }
    memory { BcaResources.scaledMemory(
        params.dynamic_memory?.MERGE_GROUP_H5AD, [cells_h5ad], task.attempt, 16, ambient_h5ad) }

    // The image METACELL2 runs in: the merge imports run_metacells.py's CellSweep lift
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_h5py_metacells_anndata_pruned:101f1effbc94defd' :
        'community.wave.seqera.io/library/python_h5py_metacells_anndata_pruned:fe586013035ea72f' }"

    input:
    // members: one [analytical id, sample, cells file name, ambient file name or ''] per member
    tuple val(meta), val(members), path(cells_h5ad), path(ambient_h5ad)

    output:
    tuple val(meta), path("${meta.id}_group_cells.h5ad"), emit: h5ad
    tuple val(meta), path("${meta.id}_group_members.tsv"), emit: members
    path "versions.yml",                                  emit: versions

    script:
    def member_args = members.collect { m ->
        "--member '${m[0]}' '${m[1]}' '${m[2]}'" + (m[3] ? " '${m[3]}'" : "")
    }.join(' \\\n        ')
    """
    merge_group_h5ad.py \\
        --group_id ${meta.id} \\
        ${member_args} \\
        --output ${meta.id}_group_cells.h5ad \\
        --members_tsv ${meta.id}_group_members.tsv \\
        --compression ${params.h5ad_compression ?: 'gzip'} \\
        --versions_yml versions.yml \\
        --process_name "${task.process}"
    """

    stub:
    """
    touch ${meta.id}_group_cells.h5ad versions.yml
    printf 'id\\tsample\\tn_cells\\n' > ${meta.id}_group_members.tsv
    """
}
