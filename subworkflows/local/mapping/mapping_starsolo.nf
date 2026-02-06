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
include { SAMTOOLS_VIEW                                     } from '../../../modules/local/tools/samtools/samtools_view/main'
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
        // Initialize variables to ensure scope visibility
        def ch_mapping_files = Channel.empty()
        def ch_starsolo_bam  = Channel.empty()
        def ch_genefull50_raw_dir = Channel.empty()
        def ch_all_outputs   = all_outputs

        // Get the first cDNA fastq file from the samplesheet to use for STAR index generation
        def ch_first_cDNA = data_output
            .map { meta, fastq_cDNA, fastq_BC_UMI, empty_list, samplesheet -> fastq_cDNA }
            .first()

        // Check if star index is provided, if not create it
        def star_index_ch
        if (params.star_index && file(params.star_index).exists()) {
            star_index_ch = Channel.value(file(params.star_index))
        } else {
            STARSOLO_INDEX(params.ref_gtf, ch_first_cDNA)
            star_index_ch = STARSOLO_INDEX.out
        }

        // Run STARsolo alignment
        STARSOLO_ALIGN(data_output, bc_whitelist, star_index_ch)
        ch_mapping_files = STARSOLO_ALIGN.out.starsolo_files

        // Check if generate BAM is enabled/disabled in the custom configurations
        if (params.star_generateBAM) {

            // Capture the BAM output if generated
            ch_starsolo_bam = STARSOLO_ALIGN.out.bam_file
            ch_genefull50_raw_dir = STARSOLO_ALIGN.out.genefull50_raw_dir

            SAMTOOLS_INDEX(ch_starsolo_bam)

            // Calculate saturation curve if perform_10x_saturate is true
            if (params.perform_10x_saturate) {
                SAMTOOLS_VIEW(STARSOLO_ALIGN.out.starsolo_files)
                SATURATION_TABLE(STARSOLO_ALIGN.out.starsolo_files, SAMTOOLS_VIEW.out.bam_file, SAMTOOLS_VIEW.out.bam_index, SAMTOOLS_VIEW.out.mapreads)
                SATURATION_PLOT(STARSOLO_ALIGN.out.starsolo_files, SATURATION_TABLE.out)
                ch_all_outputs = ch_all_outputs.mix(SATURATION_PLOT.out)
            }

            // Calculate percentages mitochondrial DNA and ribosomal RNA
            if (params.perform_featurecounts) {
                CALC_MT_RRNA(STARSOLO_ALIGN.out.starsolo_files, SAMTOOLS_INDEX.out.bam_index)
                ch_all_outputs = ch_all_outputs.mix(CALC_MT_RRNA.out)
            }

            // Conditionally run Gene Extension + Remapping branch
            if (params.perform_geneext && params.run_method != "geneext_only") {
                GENE_EXT(STARSOLO_ALIGN.out.starsolo_files, SAMTOOLS_INDEX.out.bam_index)

                // Create STAR index with extended GTF
                STARSOLO_INDEX_GENEEXT(GENE_EXT.out, ch_first_cDNA)

                // Remap with STARsolo using the extended GTF
                STARSOLO_ALIGN_GENEEXT(data_output, bc_whitelist, STARSOLO_INDEX_GENEEXT.out)
                SAMTOOLS_INDEX_GENEEXT(STARSOLO_ALIGN_GENEEXT.out.bam_file)

                // Add outputs from the Gene Extension branch to the output channel
                ch_mapping_files = ch_mapping_files.mix(STARSOLO_ALIGN_GENEEXT.out.starsolo_files)
                ch_all_outputs = ch_all_outputs.mix(SAMTOOLS_INDEX_GENEEXT.out.bam_index)

            } else if (params.run_method == "geneext_only") {
                GENE_EXT(STARSOLO_ALIGN.out.starsolo_files, SAMTOOLS_INDEX.out.bam_index)
                ch_all_outputs = ch_all_outputs.mix(GENE_EXT.out)
            }
        }

    emit:
        mapping_files       = ch_mapping_files
        all_outputs         = ch_all_outputs
        starsolo_bam        = ch_starsolo_bam
        starsolo_genefull50_raw  = ch_genefull50_raw_dir
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
