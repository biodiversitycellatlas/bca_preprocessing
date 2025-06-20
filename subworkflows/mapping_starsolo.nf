//
// Subworkflow with functionality specific to the workflow 'mapping_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { STAR_INDEX as STAR_INDEX                          } from '../modules/genindex_starsolo'
include { STAR_INDEX as STAR_INDEX_GENEEXT                  } from '../modules/genindex_starsolo'
include { MAPPING_STARSOLO as MAPPING_STARSOLO              } from '../modules/mapping_starsolo'
include { MAPPING_STARSOLO as MAPPING_STARSOLO_GENEEXT      } from '../modules/mapping_starsolo'
include { INDEX_BAM as INDEX_BAM                            } from '../modules/index_bam'
include { INDEX_BAM as INDEX_BAM_GENEEXT                    } from '../modules/index_bam'
include { SATURATION                                        } from '../modules/saturation'
include { SATURATION_PLOT                                   } from '../modules/saturation_plot'
include { CALC_MT_RRNA as CALC_MT_RRNA                      } from '../modules/calculate_mt_rrna'
include { CALC_MT_RRNA as CALC_MT_RRNA_GENEEXT              } from '../modules/calculate_mt_rrna'
include { GENE_EXT                                          } from '../modules/gene_ext'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO RUN STARSOLO
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow mapping_starsolo_workflow {
    take:
        data_output
        bc_whitelist
        all_outputs
    
    main:
        // Mapping: STARsolo
        STAR_INDEX(data_output, params.ref_gtf)
        mapping_files = MAPPING_STARSOLO(data_output, bc_whitelist, STAR_INDEX.out)
        INDEX_BAM(MAPPING_STARSOLO.out)

        // Calculate saturation
        // SATURATION(MAPPING_STARSOLO.out, INDEX_BAM.out)
        // SATURATION_PLOT(MAPPING_STARSOLO.out, SATURATION.out)
        // all_outputs.mix(SATURATION_PLOT.out)

        // Calculate percentages mitochondrial DNA and ribosomal RNA
        if (params.perform_featurecounts) {
            CALC_MT_RRNA(MAPPING_STARSOLO.out, INDEX_BAM.out)
            all_outputs = all_outputs.mix(CALC_MT_RRNA.out)
        } else {
            log.info "Skipping mtDNA/rRNA calculation as 'perform_featurecounts' is false."
        }

        // Conditionally run Gene Extension + Remapping branch
        // params.perform_geneext (boolean) must be true
        if (params.perform_geneext) {
            GENE_EXT(MAPPING_STARSOLO.out, INDEX_BAM.out)
            STAR_INDEX_GENEEXT(data_output, GENE_EXT.out)
            MAPPING_STARSOLO_GENEEXT(data_output, bc_whitelist, STAR_INDEX_GENEEXT.out)
            INDEX_BAM_GENEEXT(MAPPING_STARSOLO_GENEEXT.out)

            mapping_files = mapping_files.mix(MAPPING_STARSOLO_GENEEXT.out)
            all_outputs = all_outputs.mix(INDEX_BAM_GENEEXT.out)
        } else {
            log.info "Skipping Gene Extension steps as 'perform_geneext' is false."
        }

    emit:
        mapping_files
        all_outputs
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/