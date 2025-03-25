include { FASTQC } from '../modules/fastqc'

include { GENINDEX_STARSOLO } from '../modules/genindex_starsolo'
include { GENINDEX_STARSOLO as GENINDEX_STARSOLO_GENEEXT} from '../modules/genindex_starsolo'

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

        // Mapping: STARsolo
        GENINDEX_STARSOLO(params.ref_star_gtf)
        MAPPING_STARSOLO(data_output, GENINDEX_STARSOLO.out, params.gtf_name)
        INDEX_BAM(MAPPING_STARSOLO.out)

        // Mapping: STARsolo + Gene Extension
        GENE_EXT(MAPPING_STARSOLO.out, INDEX_BAM.out, params.gtf_name)
        GENINDEX_STARSOLO_GENEEXT(GENE_EXT.out)
        MAPPING_STARSOLO_GENEEXT(data_output, GENINDEX_STARSOLO_GENEEXT.out, "${params.gtf_name}_geneext")
        INDEX_BAM_GENEEXT(MAPPING_STARSOLO_GENEEXT.out)

        // Mapping: Alevin-fry
        // GENINDEX_ALEVIN()
        // MAPPING_ALEVIN(data_output) //, GENINDEX_ALEVIN.out.index, GENINDEX_ALEVIN.out.transcript_tsv

        // Generate mtx files
        // ...

        // Calculate saturation
        SATURATION(MAPPING_STARSOLO.out, INDEX_BAM.out, params.gtf_name)
        SATURATION_GENEEXT(MAPPING_STARSOLO_GENEEXT.out, INDEX_BAM_GENEEXT.out, "${params.gtf_name}_geneext")

        // Calculate percentages mitochondrial DNA and ribosomal RNA
        CALC_MT_RRNA(MAPPING_STARSOLO.out, INDEX_BAM.out, params.gtf_name)
        CALC_MT_RRNA_GENEEXT(MAPPING_STARSOLO_GENEEXT.out, INDEX_BAM_GENEEXT.out, "${params.gtf_name}_geneext")     

    emit:
        SATURATION_GENEEXT.out.collect()
        CALC_MT_RRNA_GENEEXT.out.collect()
}
