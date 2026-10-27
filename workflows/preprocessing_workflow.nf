//
// Workflow with functionality specific to 'main.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { bd_rhapsody_workflow          } from '../subworkflows/local/pre-processing/preprocs_bd_rhapsody'
include { parse_workflow                } from '../subworkflows/local/pre-processing/preprocs_parse_biosciences'
include { tenx_genomics_workflow        } from '../subworkflows/local/pre-processing/preprocs_10x'
include { sciRNAseq3_nogather_workflow  } from '../subworkflows/local/pre-processing/preprocs_sciRNAseq3_nogather'
include { marsseq_workflow              } from '../subworkflows/local/pre-processing/preprocs_marsseq'
include { seqspec_workflow              } from '../subworkflows/local/pre-processing/preprocs_seqspec'
include { oak_workflow                  } from '../subworkflows/local/pre-processing/preprocs_oak'

include { MERGE_FASTQS                  } from '../modules/local/custom/manipulate/merge_files/main'
include { DOWNLOAD_WHITELIST            } from '../modules/local/custom/manipulate/download_whitelist/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW TO RUN PRE-PROCESSING
        Sequencing-specific pre-processing of the data:
        - Parse Bioscience: Demultiplexing using groups of wells and mapping using split-pipe
        - BD Rhapsody: Removing variable bases and mapping using BD rhapsody pipeline
        - 10x Genomics: Mapping using CellRanger
        - OAK: Demultiplexing the aliquots on their i7/i5 index, then mapping using CellRanger
        - Sci-RNA-seq3: Pre-processing based on the sci-rocket pipeline
        - MARS-seq: Rebuilding the reads so the batch barcode joins the cell barcode and UMI
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow preprocessing_workflow {
    take:
        ch_samplesheet

    main:
        // Merge fastq files of duplicate sample IDs
        MERGE_FASTQS(ch_samplesheet)
        merged_samplesheet = MERGE_FASTQS.out.merged_files
        def ch_versions = MERGE_FASTQS.out.versions

        // Vendor pipeline results the mapping statistics are read from, as [pipeline, meta, file]
        def ch_vendor_results = channel.empty()

        // The barcode whitelist as a list of files, one per barcode segment in segment order,
        // and empty without one (e.g. MARS-seq). 'bc_whitelist' gives one or more URLs or
        // paths, separated by whitespace: URLs are downloaded, paths are staged into the
        // tasks, so a relative path is read against the launch directory like any parameter.
        def whitelist_tokens = params.bc_whitelist?.toString()?.trim()?.tokenize() ?: []
        def whitelist_urls   = whitelist_tokens.findAll { token -> token ==~ /(?i)^(https?|ftp):\/\/.+/ }
        if (whitelist_urls && whitelist_urls.size() != whitelist_tokens.size()) {
            error("'bc_whitelist' mixes URLs and paths ('${params.bc_whitelist}'): give either URLs or paths, so the segments stay in order.")
        }

        def resolved_whitelist
        if (whitelist_urls) {
            DOWNLOAD_WHITELIST(whitelist_urls.join(' '))
            ch_versions = ch_versions.mix(DOWNLOAD_WHITELIST.out.versions)
            resolved_whitelist = DOWNLOAD_WHITELIST.out.whitelist.map { files -> [files].flatten() }
        } else {
            def missing = whitelist_tokens.findAll { token -> !file(token).exists() }
            if (missing) {
                error("'bc_whitelist' names file(s) that do not exist: ${missing.join(', ')} (relative paths are read from the launch directory, ${workflow.launchDir}).")
            }
            resolved_whitelist = channel.value(whitelist_tokens.collect { token -> file(token) })
        }

        // Check for protocol and run appropriate pre-processing workflow
        if (params.protocol in ['parse_biosciences_WT_mini', 'parse_biosciences_WT']) {
            parse_workflow(merged_samplesheet)
            data_output_ch = parse_workflow.out.data_output
            bc_whitelist_ch  = resolved_whitelist
            ch_versions = ch_versions.mix(parse_workflow.out.versions)
            ch_vendor_results = parse_workflow.out.splitpipe_stats.map { meta, f -> ['splitpipe', meta, f] }

        } else if (params.protocol in ['bd_rhapsody_v1', 'bd_rhapsody_enhancedbeads']) {
            bd_rhapsody_workflow(merged_samplesheet)
            data_output_ch = bd_rhapsody_workflow.out.data_output
            bc_whitelist_ch  = resolved_whitelist
            ch_versions = ch_versions.mix(bd_rhapsody_workflow.out.versions)

        } else if (params.protocol in ['10xv1', '10xv2', '10xv3', '10xv4']) {
            tenx_genomics_workflow(merged_samplesheet)
            data_output_ch = tenx_genomics_workflow.out.data_output
            bc_whitelist_ch  = resolved_whitelist
            ch_versions = ch_versions.mix(tenx_genomics_workflow.out.versions)
            ch_vendor_results = tenx_genomics_workflow.out.cellranger_outs.map { meta, outs -> ['cellranger', meta, outs] }

        } else if (params.protocol == 'oak_v1') {
            oak_workflow(merged_samplesheet)
            data_output_ch = oak_workflow.out.data_output
            bc_whitelist_ch  = resolved_whitelist
            ch_versions = ch_versions.mix(oak_workflow.out.versions)
            ch_vendor_results = oak_workflow.out.cellranger_outs.map { meta, outs -> ['cellranger', meta, outs] }

        } else if (params.protocol == 'sciRNAseq3') {
            sciRNAseq3_nogather_workflow(merged_samplesheet, resolved_whitelist)
            data_output_ch   = sciRNAseq3_nogather_workflow.out.data_output
            // The mappers take one whitelist for every sample: the first sample's from sci-rocket
            bc_whitelist_ch  = sciRNAseq3_nogather_workflow.out.bc_whitelist.map { tup -> [tup].flatten() }.ifEmpty([]).first()
            ch_versions = ch_versions.mix(sciRNAseq3_nogather_workflow.out.versions)

        } else if (params.protocol in ['marsseq_v1', 'marsseq_v2']) {
            marsseq_workflow(merged_samplesheet)
            data_output_ch = marsseq_workflow.out.data_output
            ch_versions = ch_versions.mix(marsseq_workflow.out.versions)

            // MARS-seq has no barcode whitelist; an empty list makes STARsolo use '--soloCBwhitelist None'
            bc_whitelist_ch  = resolved_whitelist

        } else if (params.protocol == 'seqspec') {
            if (!params.seqspec_file || !file(params.seqspec_file).exists()) {
                error("'protocol' = 'seqspec' reads the read layout from 'seqspec_file', which is ${params.seqspec_file ? "'${params.seqspec_file}' and does not exist" : 'not set'}.")
            }
            seqspec_workflow(merged_samplesheet)
            data_output_ch = seqspec_workflow.out.data_output
            bc_whitelist_ch  = resolved_whitelist
            ch_versions = ch_versions.mix(seqspec_workflow.out.versions)

        } else {
            error """
            Invalid sequencing technology specified. Use one of the following parameters for 'protocol':
            - 'parse_biosciences_WT_mini'
            - 'parse_biosciences_WT'
            - 'bd_rhapsody_v1'
            - 'bd_rhapsody_enhancedbeads'
            - '10xv1'
            - '10xv2'
            - '10xv3'
            - '10xv4'
            - 'oak_v1'
            - 'sciRNAseq3'
            - 'marsseq_v1'
            - 'marsseq_v2'
            Or use 'seqspec' to specify a non-supported sequencing technique.
            """
        }

    emit:
        merged_samplesheet = merged_samplesheet
        data_output     = data_output_ch
        bc_whitelist    = bc_whitelist_ch
        vendor_results  = ch_vendor_results
        versions        = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
