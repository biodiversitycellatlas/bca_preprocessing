// ===============  QC & Mapping Workflow  ==============  \\ 
//      \\


// Import processes
include { FASTQC } from '../modules/fastqc'

include { GENINDEX_STARSOLO as GENINDEX_STARSOLO_N } from '../modules/genindex_starsolo'
include { GENINDEX_STARSOLO as GENINDEX_STARSOLO_CR } from '../modules/genindex_starsolo'
include { GENINDEX_STARSOLO as REINDEX_STARSOLO_N } from '../modules/genindex_starsolo'
include { GENINDEX_STARSOLO as REINDEX_STARSOLO_CR } from '../modules/genindex_starsolo'

include { MAPPING_STARSOLO as MAPPING_STARSOLO_N } from '../modules/mapping_starsolo'
include { MAPPING_STARSOLO as MAPPING_STARSOLO_CR } from '../modules/mapping_starsolo'
include { MAPPING_STARSOLO as REMAPPING_STARSOLO_N } from '../modules/mapping_starsolo'
include { MAPPING_STARSOLO as REMAPPING_STARSOLO_CR } from '../modules/mapping_starsolo'

include { INDEX_BAM as INDEX_BAM_N } from '../modules/index_bam'
include { INDEX_BAM as INDEX_BAM_NGE } from '../modules/index_bam'
include { INDEX_BAM as INDEX_BAM_CR } from '../modules/index_bam'
include { INDEX_BAM as INDEX_BAM_CRGE } from '../modules/index_bam'

include { SATURATION as SATURATION_N } from '../modules/saturation'
include { SATURATION as SATURATION_NGE } from '../modules/saturation'
include { SATURATION as SATURATION_CR } from '../modules/saturation'
include { SATURATION as SATURATION_CRGE } from '../modules/saturation'

include { CALC_MT_RRNA as CALC_MT_RRNA_N } from '../modules/calculate_mt_rrna'
include { CALC_MT_RRNA as CALC_MT_RRNA_NGE } from '../modules/calculate_mt_rrna'
include { CALC_MT_RRNA as CALC_MT_RRNA_CR } from '../modules/calculate_mt_rrna'
include { CALC_MT_RRNA as CALC_MT_RRNA_CRGE } from '../modules/calculate_mt_rrna'

include { GENE_EXT as GENE_EXT_N } from '../modules/gene_ext'
include { GENE_EXT as GENE_EXT_CR } from '../modules/gene_ext'


workflow QC_mapping_workflow {
    take:
        data_output

    main:
        // Quality Control
        FASTQC(data_output)

        // Mapping STARsolo
        GENINDEX_STARSOLO_N(params.ref_star_gtf, file(params.star_config_mkref_N), 'N')
        GENINDEX_STARSOLO_CR(params.ref_star_gtf, file(params.star_config_mkref_CR), 'CR')

        // Mapping: Normal
        MAPPING_STARSOLO_N(data_output, GENINDEX_STARSOLO_N.out, file(params.star_config_ED), params.barcodeDir, 'N')
        INDEX_BAM_N(MAPPING_STARSOLO_N.out)
        SATURATION_N(MAPPING_STARSOLO_N.out, INDEX_BAM_N.out)
        CALC_MT_RRNA_N(MAPPING_STARSOLO_N.out)
        
        // Mapping: Normal + Gene Extension
        GENE_EXT_N(MAPPING_STARSOLO_N.out, INDEX_BAM_N.out)
        REINDEX_STARSOLO_N(GENE_EXT_N.out, file(params.star_config_mkref_N), 'NGE')
        REMAPPING_STARSOLO_N(data_output, REINDEX_STARSOLO_N.out, file(params.star_config_ED), params.barcodeDir, 'NGE')
        INDEX_BAM_NGE(REMAPPING_STARSOLO_N.out)
        SATURATION_NGE(REMAPPING_STARSOLO_N.out, INDEX_BAM_NGE.out)
        CALC_MT_RRNA_NGE(REMAPPING_STARSOLO_N.out)

        // Mapping: CR-like
        MAPPING_STARSOLO_CR(data_output, GENINDEX_STARSOLO_CR.out, file(params.star_config_CRED), params.barcodeDemux, 'CR') // original: barcodeDemux
        INDEX_BAM_CR(MAPPING_STARSOLO_CR.out)
        SATURATION_CR(MAPPING_STARSOLO_CR.out, INDEX_BAM_CR.out)
        CALC_MT_RRNA_CR(MAPPING_STARSOLO_CR.out)

        // Mapping: CR-like + Gene Extension
        GENE_EXT_CR(MAPPING_STARSOLO_CR.out, INDEX_BAM_CR.out)
        REINDEX_STARSOLO_CR(GENE_EXT_CR.out, file(params.star_config_mkref_CR), 'CRGE')
        REMAPPING_STARSOLO_CR(data_output, REINDEX_STARSOLO_CR.out, file(params.star_config_CRED), params.barcodeDemux, 'CRGE')
        INDEX_BAM_CRGE(REMAPPING_STARSOLO_CR.out)
        SATURATION_CRGE(REMAPPING_STARSOLO_CR.out, INDEX_BAM_CRGE.out)
        CALC_MT_RRNA_CRGE(REMAPPING_STARSOLO_CR.out)

        // Generate mtx files
        // ...

    emit:
        // MAPPING_STARSOLO_N.out.collect()
        // MAPPING_STARSOLO_CR.out.collect()
        calc_mt_rrna_n = CALC_MT_RRNA_N.out.collect()
        calc_mt_rrna_nge = CALC_MT_RRNA_NGE.out.collect()
        calc_mt_rrna_cr = CALC_MT_RRNA_CR.out.collect()
        calc_mt_rrna_crge = CALC_MT_RRNA_CRGE.out.collect()

        // intron.mtx, exon.mtx, fullgenome.mtx
}
