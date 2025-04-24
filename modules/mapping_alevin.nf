process MAPPING_ALEVIN {
    publishDir "${params.output_dir}/mapping_alevin/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    debug true

    input:
    tuple val(meta), path(fastq_files)
    path(index)

    output:
    path("*")
       
    script:
    def fastq_list = fastq_files instanceof List ? fastq_files : [fastq_files]
    def r1_fastq = fastq_list.find { it.name.contains('_R1') }
    def r2_fastq = fastq_list.find { it.name.contains('_R2') }

    // If protocol is "bd_rhapsody", then cDNA = R2 and CB/UMI = R1
    // Else by default cDNA = R1 and CB/UMI = R2
    def cDNA_read
    def CBUMI_read
    if (params.protocol.toLowerCase().contains("bd_rhapsody")) {
        cDNA_read = r2_fastq
        CBUMI_read = r1_fastq
        bc_geom = "1[0-8,13-21,26-34]"
        umi_geom = "1[35-42]"
        read_geom = "2[1-end]"

    } else if (params.protocol.toLowerCase().contains("parse_biosciences")) {
        cDNA_read = r1_fastq
        CBUMI_read = r2_fastq
        bc_geom = "2[51-58,31-38,11-18]"
        umi_geom = "2[1-10]"
        read_geom = "1[1-end]"

    } else {
        cDNA_read = r2_fastq
        CBUMI_read = r1_fastq
        bc_geom = "1[1-16]"
        umi_geom = "1[17-29]"
        read_geom = "2[1-end]"
    }
    """
    echo "\n\n==================  GENOME INDEX ALEVIN =================="
    echo "Conda environment: \$CONDA_DEFAULT_ENV"
    echo "Sample ID: ${meta}"
    echo "Index: ${index}"
    echo "Reference fasta: ${params.ref_fasta}"
    echo "Reference ref_gtf: ${params.ref_gtf}"
    echo "cDNA read: ${cDNA_read}"
    echo "CB/UMI read: ${CBUMI_read}"


    echo "\n\n-------------  Salmon Alevin -------------------"
    salmon alevin \\
        -i salmon_index \\
        -l A \\
        -1 ${cDNA_read} \\
        -2 ${CBUMI_read} \\
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
    txp2gene=\$(ls ./splici_index_reference/*t2g_3col.tsv)

    alevin-fry quant \\
        -m \${txp2gene} \\
        -i ./${meta.id}_out_permit_knee  \\
        -o ./${meta.id}_counts \\
        -t 16 \\
        -r cr-like-em \\
        --use-mtx 
    """
}
