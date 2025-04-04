process MAPPING_ALEVIN {
    publishDir "${params.resDir}/mapping_alevin/${sample_id}", mode: 'copy'
    tag "${sample_id}"
    debug true

    input:
    tuple val(sample_id), path(fastq_files)
    path(index)

    output:
    path("*")
       
    script:
    def fastq_list = fastq_files instanceof List ? fastq_files : [fastq_files]
    def r1_fastq = fastq_list.find { it.name.contains('_R1') }
    def r2_fastq = fastq_list.find { it.name.contains('_R2') }

    // If seqTech is "bd_rhapsody", then cDNA = R2 and CB/UMI = R1
    // Else by default cDNA = R1 and CB/UMI = R2
    def cDNA_read
    def CBUMI_read
    if (params.seqTech.toLowerCase().contains("bd_rhapsody")) {
        cDNA_read = r2_fastq
        CBUMI_read = r1_fastq
    } else {
        cDNA_read = r1_fastq
        CBUMI_read = r2_fastq
    }
    """
    echo "\n\n==================  GENOME INDEX ALEVIN =================="
    echo "Sample ID: ${sample_id}"
    echo "Index: ${index}"
    echo "Reference fasta: ${params.ref_fasta}"
    echo "Reference ref_star_gtf: ${params.ref_star_gtf}"
    echo "cDNA read: ${cDNA_read}"
    echo "CB/UMI read: ${CBUMI_read}"


    echo "\n\n-------------  Salmon Alevin -------------------"
    salmon alevin \\
        -i salmon_index \\
        -l A \\
        -1 ${cDNA_read} \\
        -2 ${CBUMI_read} \\
        -p 32 \\
        --bc-geometry 2[51-58,31-38,11-18] \\
        --umi-geo 2[1-10] \\
        --read-geo 1[1-end] \\
        -o ./${sample_id}_run \\
        --justAlign

    echo "\n\n-------------  generate permit -------------------"
    alevin-fry generate-permit-list \\
        -i ./${sample_id}_run \\
        -d both \\
        --output-dir ./${sample_id}_out_permit_knee \\
        -k

    echo "\n\n-------------  collate -------------------"
    alevin-fry collate \\
        -i ./${sample_id}_out_permit_knee \\
        -t 16 \\
        -r ./${sample_id}_run  

    echo "\n\n-------------  quant -------------------"
    txp2gene=\$(ls ./splici_index_reference/*.tsv)

    alevin-fry quant \\
        -m \${txp2gene} \\
        -i ./${sample_id}_out_permit_knee  \\
        -o ./${sample_id}_counts \\
        -t 16 \\
        -r full \\
        --use-mtx 
    """
}
