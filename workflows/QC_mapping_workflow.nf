include { FASTQC } from '../modules/fastqc'

include { GENINDEX_STARSOLO } from '../modules/genindex_starsolo'
include { MAPPING_STARSOLO as MAPPING_STARSOLO } from '../modules/mapping_starsolo'
include { MAPPING_STARSOLO as MAPPING_STARSOLO_GENEEXT } from '../modules/mapping_starsolo'

include { GENINDEX_ALEVIN } from '../modules/genindex_alevin'
include { MAPPING_ALEVIN } from '../modules/mapping_alevin'

include { INDEX_BAM as INDEX_BAM } from '../modules/index_bam'
include { INDEX_BAM as INDEX_BAM_GENEEXT } from '../modules/index_bam'

include { SATURATION as SATURATION } from '../modules/saturation'
include { SATURATION as SATURATION_GENEEXT } from '../modules/saturation'

include { CALC_MT_RRNA as CALC_MT_RRNA } from '../modules/calculate_mt_rrna'
include { CALC_MT_RRNA as CALC_MT_RRNA_GENEEXT } from '../modules/calculate_mt_rrna'

include { GENE_EXT } from '../modules/gene_ext'


workflow QC_mapping_workflow {
    take:
        data_output

    main:
        // Quality Control
        FASTQC(data_output)

        // Create STAR index
        GENINDEX_STARSOLO()

        // Mapping: STARsolo
        MAPPING_STARSOLO(data_output, GENINDEX_STARSOLO.out)
        INDEX_BAM(MAPPING_STARSOLO.out)
        SATURATION(MAPPING_STARSOLO.out, INDEX_BAM.out)
        CALC_MT_RRNA(MAPPING_STARSOLO.out, INDEX_BAM.out)

        // Mapping: STARsolo + Gene Extension
        // GENE_EXT_PAIRED(MAPPING_STARSOLO_PAIRED.out, INDEX_BAM_PAIRED.out)
        // INDEX_BAM_PAIRED_GE(GENE_EXT_PAIRED.out)
        // REMAPPING_STARSOLO_PAIRED(data_output, INDEX_BAM_PAIRED_GE.out)
        // SATURATION_PAIRED_GE(REMAPPING_STARSOLO_PAIRED.out, INDEX_BAM_PAIRED_GE.out)
        // CALC_MT_RRNA_PAIRED_GE(REMAPPING_STARSOLO_PAIRED.out, INDEX_BAM_PAIRED_GE.out)     

        // Mapping: Alevin-fry
        // GENINDEX_ALEVIN()
        // MAPPING_ALEVIN(data_output) //, GENINDEX_ALEVIN.out.index, GENINDEX_ALEVIN.out.transcript_tsv

        // Generate mtx files
        // ...

    emit:
        MAPPING_STARSOLO.out.collect()
}
