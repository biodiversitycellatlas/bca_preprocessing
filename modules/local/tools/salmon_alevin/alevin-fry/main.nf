process ALEVIN_FRY {
    publishDir "${params.outdir}/mapping_alevin/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    label 'process_high'


    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/alevin-fry_piscem_salmon_simpleaf_pruned:c71cfb476b003414"

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file)
    path(bc_whitelist)
    path(splici_index_reference)
    path(salmon_index)

    output:
    tuple val(meta), path("${meta.id}_*"), emit: mapping_files
    tuple val(meta), path("${meta.id}_run/aux_info/meta_info.json"),    emit: af_meta_info
    tuple val(meta), path("${meta.id}_counts/quant.json"),              emit: af_quant_json
    tuple val(meta), path("${meta.id}_counts/alevin"),                  emit: af_mtx
    tuple val(meta), path("${meta.id}_counts/cell_meta.tsv"),           emit: af_cell_meta, optional: true
    path "versions.yml",                                                emit: versions

    script:
    // Retrieve settings from custom parameters if set, otherwise from conf/seqtech_parameters.config
    def bc_geom = params.alevin_bc_geometry ?: params.seqtech_parameters[params.protocol].alevin_bc_geometry
    def umi_geom = params.alevin_umi_geometry ?: params.seqtech_parameters[params.protocol].alevin_umi_geometry
    def read_geom = params.alevin_read_geometry ?: params.seqtech_parameters[params.protocol].alevin_read_geometry

    if (!bc_geom || !umi_geom || !read_geom) {
        error "No alevin geometry defined for protocol '${params.protocol}'. Set 'alevin_bc_geometry', 'alevin_umi_geometry' and 'alevin_read_geometry' in the configuration file."
    }
    """
    echo "\n\n==================  ALEVIN-FRY =================="
    echo "Sample ID: ${meta}"
    echo "Salmon Index: ${salmon_index}"
    echo "Reference fasta: ${params.ref_fasta}"
    echo "Reference ref_gtf: ${params.ref_gtf}"
    echo "cDNA read: ${fastq_cDNA}"
    echo "CB/UMI read: ${fastq_BC_UMI}"
    echo "Geometry (bc / umi / read): ${bc_geom} / ${umi_geom} / ${read_geom}"


    echo "\n\n-------------  Salmon Alevin -------------------"
    salmon alevin \\
        -i ${salmon_index} \\
        -l A \\
        -1 ${fastq_BC_UMI} \\
        -2 ${fastq_cDNA} \\
        -p 32 \\
        --bc-geometry "${bc_geom}" \\
        --umi-geo "${umi_geom}" \\
        --read-geo "${read_geom}" \\
        -o ./${meta.id}_run \\
        --justAlign

    echo "\n\n-------------  generate permit -------------------"
    alevin-fry generate-permit-list \\
        -i ./${meta.id}_run \\
        -d both \\
        --output-dir ./${meta.id}_out_permit_knee \\
        -k

    echo "\n\n-------------  collate -------------------"
    alevin-fry collate \\
        -i ./${meta.id}_out_permit_knee \\
        -t 16 \\
        -r ./${meta.id}_run

    echo "\n\n-------------  quant -------------------"
    alevin-fry quant \\
        -m ${splici_index_reference}/*t2g_3col.tsv \\
        -i ./${meta.id}_out_permit_knee  \\
        -o ./${meta.id}_counts \\
        -t 16 \\
        -r cr-like-em \\
        --use-mtx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(salmon --version | sed 's/salmon //')
        alevin-fry: \$(alevin-fry --version | sed 's/alevin-fry //')
    END_VERSIONS
    """
}
