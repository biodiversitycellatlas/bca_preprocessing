//
// Workflow with functionality specific to 'main.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { mapping_starsolo_workflow } from '../subworkflows/local/mapping/mapping_starsolo'
include { mapping_alevin_workflow   } from '../subworkflows/local/mapping/mapping_alevin'

include { FASTQC                    } from '../modules/local/tools/fastqc/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW TO RUN MAPPING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow QC_mapping_workflow {
    take:
        data_output
        bc_whitelist

    main:
        // Initialize reporting channels
        def ch_mapping_files = Channel.empty()
        def ch_starsolo_bam  = Channel.empty()
        def ch_star_solodir  = Channel.empty()
        def ch_starsolo_genefull50_raw  = Channel.empty()
        def ch_sat_imgs         = Channel.empty()
        def ch_sat_res_imgs     = Channel.empty()
        def ch_sat_logs         = Channel.empty()
        def ch_star_umi         = Channel.empty()
        def ch_star_log         = Channel.empty()
        def ch_star_final_log   = Channel.empty()
        def ch_star_summaries   = Channel.empty()
        def ch_star_cellreads  = Channel.empty()
        def ch_alevin_meta_info = Channel.empty()
        def ch_alevin_quant_json = Channel.empty()
        def ch_alevin_cell_meta = Channel.empty()
        def ch_featurecounts    = Channel.empty()
        def ch_pavian_sankey    = Channel.empty()

        // Quality Control
        FASTQC(data_output)

        // Mapping: STARsolo, Alevin-fry, or both
        if (params.mapping_software == "starsolo") {
            mapping_starsolo_workflow(data_output, bc_whitelist)

            ch_mapping_files         =  mapping_starsolo_workflow.out.mapping_files
            ch_starsolo_bam          =  mapping_starsolo_workflow.out.starsolo_bam
            ch_star_solodir          =  mapping_starsolo_workflow.out.star_solodir
            ch_starsolo_genefull50_raw  = mapping_starsolo_workflow.out.starsolo_genefull50_raw
            ch_sat_imgs              =  mapping_starsolo_workflow.out.saturation_imgs
            ch_sat_res_imgs          =  mapping_starsolo_workflow.out.saturation_residual_imgs
            ch_sat_logs              =  mapping_starsolo_workflow.out.saturation_logs
            ch_star_umi              =  mapping_starsolo_workflow.out.star_umipercell
            ch_star_log              =  mapping_starsolo_workflow.out.star_log
            ch_star_final_log        =  mapping_starsolo_workflow.out.star_final_log
            ch_star_summaries        =  mapping_starsolo_workflow.out.star_summaries
            ch_star_cellreads        =  mapping_starsolo_workflow.out.star_cellreads
            ch_featurecounts         =  mapping_starsolo_workflow.out.featurecount_txt
            ch_pavian_sankey         =  mapping_starsolo_workflow.out.pavian_sankey

        } else if (params.mapping_software == "alevin") {
            mapping_alevin_workflow(data_output, bc_whitelist)
            ch_mapping_files = mapping_alevin_workflow.out.mapping_files
            ch_alevin_meta_info = mapping_alevin_workflow.out.af_meta_info
            ch_alevin_quant_json = mapping_alevin_workflow.out.af_quant_json
            ch_alevin_cell_meta = mapping_alevin_workflow.out.af_cell_meta

        } else if (params.mapping_software == "both") {
            mapping_starsolo_workflow(data_output, bc_whitelist)
            mapping_alevin_workflow(data_output, bc_whitelist)

            ch_mapping_files         = mapping_alevin_workflow.out.mapping_files.mix(mapping_starsolo_workflow.out.mapping_files)
            ch_starsolo_bam          = mapping_starsolo_workflow.out.starsolo_bam
            ch_star_solodir          =  mapping_starsolo_workflow.out.star_solodir
            ch_starsolo_genefull50_raw  = mapping_starsolo_workflow.out.starsolo_genefull50_raw
            ch_sat_imgs              =  mapping_starsolo_workflow.out.saturation_imgs
            ch_sat_res_imgs          =  mapping_starsolo_workflow.out.saturation_residual_imgs
            ch_sat_logs              =  mapping_starsolo_workflow.out.saturation_logs
            ch_star_umi              =  mapping_starsolo_workflow.out.star_umipercell
            ch_star_log              =  mapping_starsolo_workflow.out.star_log
            ch_star_final_log        =  mapping_starsolo_workflow.out.star_final_log
            ch_star_summaries        =  mapping_starsolo_workflow.out.star_summaries
            ch_star_cellreads        =  mapping_starsolo_workflow.out.star_cellreads
            ch_featurecounts         =  mapping_starsolo_workflow.out.featurecount_txt
            ch_pavian_sankey         =  mapping_starsolo_workflow.out.pavian_sankey

        } else {
            error "Invalid mapping software specified. Use one of the following parameters: 'starsolo', 'alevin' or 'both'."
        }

    emit:
        mapping_files = ch_mapping_files
        starsolo_bam  = ch_starsolo_bam
        star_solodir  = ch_star_solodir
        starsolo_genefull50_raw  = ch_starsolo_genefull50_raw

        // REPORT REQUIRED
        saturation_imgs          = ch_sat_imgs
        saturation_residual_imgs = ch_sat_res_imgs
        saturation_logs          = ch_sat_logs
        star_umipercell          = ch_star_umi
        star_log                 = ch_star_log
        star_final_log           = ch_star_final_log
        star_summaries           = ch_star_summaries
        star_cellreads           = ch_star_cellreads
        af_meta_info             = ch_alevin_meta_info
        af_quant_json            = ch_alevin_quant_json
        af_cell_meta             = ch_alevin_cell_meta
        featurecount_txt         = ch_featurecounts
        pavian_sankey            = ch_pavian_sankey
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
