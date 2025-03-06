// ===============  QC & Mapping Workflow  ==============  \\ 
//      \\


// Import processes
include { FASTQC } from '../modules/fastqc'

include { GENINDEX_STARSOLO } from '../modules/genindex_starsolo'

include { MAPPING_STARSOLO_SINGLE as MAPPING_STARSOLO_SINGLE } from '../modules/mapping_starsolo_single'
include { MAPPING_STARSOLO_SINGLE as REMAPPING_STARSOLO_SINGLE } from '../modules/mapping_starsolo_single'

include { MAPPING_STARSOLO_PAIRED as MAPPING_STARSOLO_PAIRED } from '../modules/mapping_starsolo_paired'
include { MAPPING_STARSOLO_PAIRED as REMAPPING_STARSOLO_PAIRED } from '../modules/mapping_starsolo_paired'

include { GENINDEX_ALEVIN } from '../modules/genindex_alevin'
include { MAPPING_ALEVIN } from '../modules/mapping_alevin'

include { INDEX_BAM as INDEX_BAM_SINGLE } from '../modules/index_bam'
include { INDEX_BAM as INDEX_BAM_SINGLE_GE } from '../modules/index_bam'
include { INDEX_BAM as INDEX_BAM_PAIRED } from '../modules/index_bam'
include { INDEX_BAM as INDEX_BAM_PAIRED_GE } from '../modules/index_bam'

include { SATURATION as SATURATION_SINGLE } from '../modules/saturation'
include { SATURATION as SATURATION_SINGLE_GE } from '../modules/saturation'
include { SATURATION as SATURATION_PAIRED } from '../modules/saturation'
include { SATURATION as SATURATION_PAIRED_GE } from '../modules/saturation'

include { CALC_MT_RRNA as CALC_MT_RRNA_SINGLE } from '../modules/calculate_mt_rrna'
include { CALC_MT_RRNA as CALC_MT_RRNA_SINGLE_GE } from '../modules/calculate_mt_rrna'
include { CALC_MT_RRNA as CALC_MT_RRNA_PAIRED } from '../modules/calculate_mt_rrna'
include { CALC_MT_RRNA as CALC_MT_RRNA_PAIRED_GE } from '../modules/calculate_mt_rrna'

include { GENE_EXT as GENE_EXT_SINGLE } from '../modules/gene_ext'
include { GENE_EXT as GENE_EXT_PAIRED } from '../modules/gene_ext'


workflow QC_mapping_workflow {
    take:
        data_output

    main:
        // Quality Control
        FASTQC(data_output)

        // Create STAR index
        GENINDEX_STARSOLO()

        // ----------- Paired end workflow -----------
        MAPPING_STARSOLO_PAIRED(data_output, GENINDEX_STARSOLO.out)
        INDEX_BAM_PAIRED(MAPPING_STARSOLO_PAIRED.out)
        // SATURATION_PAIRED(MAPPING_STARSOLO_PAIRED.out, INDEX_BAM_PAIRED.out)
        // CALC_MT_RRNA_PAIRED(MAPPING_STARSOLO_PAIRED.out, INDEX_BAM_PAIRED.out)

        // Mapping: Paired + Gene Extension
        // GENE_EXT_PAIRED(MAPPING_STARSOLO_PAIRED.out, INDEX_BAM_PAIRED.out)
        // INDEX_BAM_PAIRED_GE(GENE_EXT_PAIRED.out)
        // REMAPPING_STARSOLO_PAIRED(data_output, INDEX_BAM_PAIRED_GE.out)
        // SATURATION_PAIRED_GE(REMAPPING_STARSOLO_PAIRED.out, INDEX_BAM_PAIRED_GE.out)
        // CALC_MT_RRNA_PAIRED_GE(REMAPPING_STARSOLO_PAIRED.out, INDEX_BAM_PAIRED_GE.out)


        // ----------- Single end workflow -----------
        // MAPPING_STARSOLO_SINGLE(data_output, GENINDEX_STARSOLO.out)
        // INDEX_BAM_SINGLE(MAPPING_STARSOLO_SINGLE.out)

        // Mapping: Single + Gene Extension
        // GENE_EXT_SINGLE(MAPPING_STARSOLO_SINGLE.out, INDEX_BAM_SINGLE.out)
        // REINDEX_STARSOLO_SINGLE(GENE_EXT_SINGLE.out)
        // REMAPPING_STARSOLO_SINGLE(data_output, REINDEX_STARSOLO_SINGLE.out)        

        // Mapping: Alevin-fry
        // GENINDEX_ALEVIN()
        // MAPPING_ALEVIN(data_output) //, GENINDEX_ALEVIN.out.index, GENINDEX_ALEVIN.out.transcript_tsv

        // Generate mtx files
        // ...

    emit:
        // SATURATION_SINGLE_GE.out.collect()
        // SATURATION_PAIRED.out.collect()
        // calc_mt_rrna_n = CALC_MT_RRNA_SINGLE.out.collect()
        // calc_mt_rrna_nge = CALC_MT_RRNA_SINGLE_GE.out.collect()
        INDEX_BAM_PAIRED.out.collect()
        // calc_mt_rrna_crge = CALC_MT_RRNA_PAIRED_GE.out.collect()

        // MAPPING_ALEVIN.out.collect()
        // intron.mtx, exon.mtx, fullgenome.mtx
}
