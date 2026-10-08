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
        // With a single sample there is nothing to balance, so its reads are kept in full.
        def ch_bams = ch_starsolo_bam
        def ch_versions = channel.empty()
        if (params.geneext_subsample_nreads != 0) {
            ch_by_count = ch_starsolo_bam
                .combine(ch_starsolo_bam.count())
                .branch { meta, bam, n_samples ->
                    subsample: n_samples > 1
                        return [meta, bam]
                    keep: true
                        return [meta, bam]
                }

            SAMTOOLS_SUBSAMPLE(ch_by_count.subsample)
            ch_bams = SAMTOOLS_SUBSAMPLE.out.subsampled_bam.mix(ch_by_count.keep)
            ch_versions = ch_versions.mix(SAMTOOLS_SUBSAMPLE.out.versions)
        }

        // Extract the BAM files and collect them into a single list
        bams_to_merge = ch_bams
            .map { meta, bam -> bam }
            .collect()

        // Merge all BAMs
        SAMTOOLS_MERGE(bams_to_merge)

        // Run gene extension using GeneExt
        GENE_EXT(SAMTOOLS_MERGE.out.merged_bam, SAMTOOLS_MERGE.out.merged_bai)
        ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions, GENE_EXT.out.versions)

    emit:
        ref_gtf         = GENE_EXT.out.gtf
        report          = GENE_EXT.out.report
        geneext_log     = GENE_EXT.out.log
        versions        = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
