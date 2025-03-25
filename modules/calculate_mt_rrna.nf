process CALC_MT_RRNA {   
    publishDir "${params.resDir}/rRNA_mtDNA/${config}", mode: 'copy'
    debug true

    input:
    tuple val(sample_id), path(mapping_files)
    file(bam_index)
    val(config)

    output:
    path("*")     
       
    script:
    """
    echo "\n\n==================  CALCULATION rRNA & mtDNA =================="

    bam_file=\$(ls ${sample_id}_Aligned.sortedByCoord.out.bam | head -n 1)

    # Create output file
    outfile="${sample_id}_${config}_mt_rrna_metrics.txt"
    echo -e "Metric\\tCount" > \$outfile
    echo -e "Sample_ID\\t${sample_id}" >> \$outfile
    echo -e "GTF_File\\t${params.ref_star_gtf}" >> \$outfile
    echo -e "rRNA_Grep\\t${params.grep_rrna}" >> \$outfile
    echo -e "MT_Contig\\t${params.mt_contig}" >> \$outfile

    # Total number of mapped reads
    mapped=\$(samtools view -F 4 \$bam_file | wc -l)
    echo -e "Mapped_Reads\\t\$mapped" >> \$outfile

    # Number of ribosomal RNA reads
    rrna=\$(samtools view \$bam_file | grep ${params.grep_rrna} | wc -l)
    echo -e "rRNA_Reads\\t\$rrna" >> \$outfile

    # Number of reads mapping to mtDNA contig
    mt=\$(samtools view \$bam_file | grep ${params.mt_contig} | wc -l)
    echo -e "mtDNA_Reads\\t\$mt" >> \$outfile

    # Multimapped reads (primary alignment)
    samtools view -h -F 256 \$bam_file | grep -E "^\\@|NH:i:[2-9]" | samtools view -b -o multimapped_primealign.bam
    mt_multi1=\$(samtools view multimapped_primealign.bam | grep ${params.mt_contig} | wc -l)
    echo -e "mtDNA_Primary_Multimapped\\t\$mt_multi1" >> \$outfile

    # Multimapped reads (all alignments)
    samtools view -h \$bam_file | grep -E "^\\@|NH:i:[2-9]" | samtools view -b -o multimapped_allalign.bam
    mt_multi2=\$(samtools view multimapped_allalign.bam | grep ${params.mt_contig} | wc -l)
    echo -e "mtDNA_All_Multimapped\\t\$mt_multi2" >> \$outfile

    """

}
