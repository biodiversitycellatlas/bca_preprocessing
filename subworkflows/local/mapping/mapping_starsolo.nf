//
// Subworkflow with functionality specific to the workflow 'mapping_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { STARSOLO_INDEX as STARSOLO_INDEX                  } from '../../../modules/local/tools/star/starsolo_genome_generate/main'
include { STARSOLO_INDEX as STARSOLO_INDEX_GENEEXT          } from '../../../modules/local/tools/star/starsolo_genome_generate/main'
include { STARSOLO_ALIGN as STARSOLO_ALIGN                  } from '../../../modules/local/tools/star/starsolo_align/main'
include { STARSOLO_ALIGN as STARSOLO_ALIGN_GENEEXT          } from '../../../modules/local/tools/star/starsolo_align/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX                  } from '../../../modules/local/tools/samtools/samtools_index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_GENEEXT          } from '../../../modules/local/tools/samtools/samtools_index/main'
include { SATURATION_TABLE                                  } from '../../../modules/local/tools/10x_saturate/saturation_table/main'
include { SATURATION_PLOT                                   } from '../../../modules/local/tools/10x_saturate/plot_curve/main'
include { CALC_MT_RRNA as CALC_MT_RRNA                      } from '../../../modules/local/tools/featurecounts/main'
include { CALC_MT_RRNA as CALC_MT_RRNA_GENEEXT              } from '../../../modules/local/tools/featurecounts/main'
include { GENE_EXT                                          } from '../../../modules/local/tools/geneext/main'


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
        STARSOLO_INDEX(data_output, params.ref_gtf)
        mapping_files = STARSOLO_ALIGN(data_output, bc_whitelist, STARSOLO_INDEX.out)
        SAMTOOLS_INDEX(STARSOLO_ALIGN.out)

        // Calculate saturation
        SATURATION_TABLE(STARSOLO_ALIGN.out, SAMTOOLS_INDEX.out)
        SATURATION_PLOT(STARSOLO_ALIGN.out, SATURATION_TABLE.out)
        all_outputs.mix(SATURATION_PLOT.out)

        // Calculate percentages mitochondrial DNA and ribosomal RNA
        if (params.perform_featurecounts) {
            CALC_MT_RRNA(STARSOLO_ALIGN.out, SAMTOOLS_INDEX.out)
            all_outputs = all_outputs.mix(CALC_MT_RRNA.out)
        }

        // Conditionally run Gene Extension + Remapping branch
        if (params.perform_geneext) {
            GENE_EXT(STARSOLO_ALIGN.out, SAMTOOLS_INDEX.out)

            // Create STAR index with extended GTF
            STARSOLO_INDEX_GENEEXT(data_output, GENE_EXT.out)

            // Remap with STARsolo using the extended GTF
            STARSOLO_ALIGN_GENEEXT(data_output, bc_whitelist, STARSOLO_INDEX_GENEEXT.out)
            SAMTOOLS_INDEX_GENEEXT(STARSOLO_ALIGN_GENEEXT.out)

            // Add outputs from the Gene Extension branch to the output channel
            mapping_files = mapping_files.mix(STARSOLO_ALIGN_GENEEXT.out)
            all_outputs = all_outputs.mix(SAMTOOLS_INDEX_GENEEXT.out)

        } else if (params.run_method == "geneext_only") {
            GENE_EXT(STARSOLO_ALIGN.out, SAMTOOLS_INDEX.out)
            all_outputs = all_outputs.mix(GENE_EXT.out)
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
