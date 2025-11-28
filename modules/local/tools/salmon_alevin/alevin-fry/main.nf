process ALEVIN_FRY {
    publishDir "${params.outdir}/mapping_alevin/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    label 'process_high'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/simpleaf:0.19.4--ha6fb395_0':
        'quay.io/biocontainers/simpleaf:0.19.4--ha6fb395_0' }"

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file)
    path bc_whitelist
    path(splici_index_reference)
    path(salmon_index)

    output:
    path("${meta.id}_*")

    script:
    // If protocol is "bd_rhapsody", then cDNA = R2 and CB/UMI = R1
    // Else by default cDNA = R1 and CB/UMI = R2
    if (params.protocol.toLowerCase().contains("bd_rhapsody")) {
        bc_geom = "1[0-8,13-21,26-34]"
        umi_geom = "1[35-42]"
        read_geom = "2[1-end]"

    } else if (params.protocol.toLowerCase().contains("parse_biosciences")) {
        bc_geom = "1[51-58,31-38,11-18]"
        umi_geom = "1[1-10]"
        read_geom = "2[1-end]"

    } else {
        bc_geom = "1[1-16]"
        umi_geom = "1[17-28]"
        read_geom = "2[1-end]"
    }
    """
    echo "\n\n==================  ALEVIN-FRY =================="
    echo "Sample ID: ${meta}"
    echo "Salmon Index: ${salmon_index}"
    echo "Reference fasta: ${params.ref_fasta}"
    echo "Reference ref_gtf: ${params.ref_gtf}"
    echo "cDNA read: ${fastq_cDNA}"
    echo "CB/UMI read: ${fastq_BC_UMI}"


    echo "\n\n-------------  Salmon Alevin -------------------"
    salmon alevin \\
        -i ${salmon_index} \\
        -l A \\
        -1 ${fastq_BC_UMI} \\
        -2 ${fastq_cDNA} \\
        -p 32 \\
        --bc-geometry ${bc_geom} \\
        --umi-geo ${umi_geom} \\
        --read-geo ${read_geom} \\
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
    """
}
