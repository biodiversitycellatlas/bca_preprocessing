process CALC_MT_RRNA {   
    publishDir "${params.resDir}/rRNA_mtDNA/${sample_id}", mode: 'copy'
    debug true

    input:
    tuple val(sample_id), path(mapping_files)
    file(bam_index)

    output:
    path("*")     
       
    script:
    """
    echo "\n\n==================  CALCULATION rRNA & mtDNA =================="
    echo "Sample ID: ${sample_id}"
    echo "GTF: ${params.ref_star_gtf}"
    echo "rRNA: ${params.grep_rrna}"
    echo "MT contig: ${params.mt_contig}"

    bam_file=\$(ls ${sample_id}_Aligned.sortedByCoord.out.bam | head -n 1)

    # Total number of mapped reads
    samtools view -F 4 \${bam_file} | wc -l

    # Number of ribosomal RNA reads
    samtools view \${bam_file} | grep ${params.grep_rrna} | wc -l

    # Number of reads mapping to mtDNA contig
    samtools view \${bam_file} | grep ${params.mt_contig} | wc -l

    # Create seperate BAM with filtered results - filter for multimapped reads only present in primary alignment
    samtools view -h -F 256 \${bam_file} | grep -E "^\\@|NH:i:[2-9]" | samtools view -b -o multimapped_primealign.bam

    # Number of multimapped reads (primary alignment) mapping to mtDNA contig
    samtools view multimapped_primealign.bam | grep ${params.mt_contig} | wc -l
    
    # Create seperate BAM with filtered results - filter for multimapped reads present in both secondary and primary alignments
    samtools view -h \${bam_file} | grep -E "^\\@|NH:i:[2-9]" | samtools view -b -o multimapped_allalign.bam

    # Number of multimapped reads (primary + secondary alignment) mapping to mtDNA contig
    samtools view multimapped_allalign.bam | grep ${params.mt_contig} | wc -l

    """
}
