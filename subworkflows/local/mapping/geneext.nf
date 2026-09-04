//
// Subworkflow with functionality specific to the workflow 'mapping_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SAMTOOLS_SUBSAMPLE                                } from '../../../modules/local/tools/samtools/samtools_subsample/main'
include { SAMTOOLS_MERGE                                    } from '../../../modules/local/tools/samtools/samtools_merge/main'
include { GENE_EXT                                          } from '../../../modules/local/tools/geneext/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO EXTEND THE GENE ANNOTATION
        GeneExt is run once on the merged alignments of every sample, since pooling gives
        enough 3' coverage to support extensions a single sample would not. Pooling also
        means a deeply sequenced sample would otherwise decide the extension for all of
        them, so each sample is capped at the same read count before the merge.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow geneext_workflow {
    take:
        ch_starsolo_bam

    main:
        // Cap every sample at the same number of reads, so no one sample is overrepresented
        // in the extended annotation. Set geneext_subsample_nreads = 0 to merge as sequenced.
        def ch_bams = ch_starsolo_bam
        if (params.geneext_subsample_nreads) {
            SAMTOOLS_SUBSAMPLE(ch_starsolo_bam)
            ch_bams = SAMTOOLS_SUBSAMPLE.out.subsampled_bam
        }

        // Extract the BAM files and collect them into a single list
        bams_to_merge = ch_bams
            .map { meta, bam -> bam }
            .collect()

        // Merge all BAMs
        SAMTOOLS_MERGE(bams_to_merge)

        // Run gene extension using GeneExt
        GENE_EXT(SAMTOOLS_MERGE.out.merged_bam, SAMTOOLS_MERGE.out.merged_bai)

    emit:
        ref_gtf         = GENE_EXT.out.gtf
        report          = GENE_EXT.out.report
        geneext_log     = GENE_EXT.out.log
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
