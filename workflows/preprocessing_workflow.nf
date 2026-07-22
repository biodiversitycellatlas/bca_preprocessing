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
include { seqspec_workflow              } from '../subworkflows/local/pre-processing/preprocs_seqspec'

include { MERGE_FASTQS                  } from '../modules/local/custom/manipulate/merge_files/main'
include { DOWNLOAD_WHITELIST            } from '../modules/local/custom/manipulate/download_whitelist/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW TO RUN PRE-PROCESSING
        Sequencing-specific pre-processing of the data:
        - Parse Bioscience: Demultiplexing using groups of wells and mapping using split-pipe
        - BD Rhapsody: Removing variable bases and mapping using BD rhapsody pipeline
        - 10xv3, OAK seq & Ultima Genomics : Mapping using CellRanger
        - Sci-RNA-seq3: Pre-processing based on the sci-rocket pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow preprocessing_workflow {
    take:
        ch_samplesheet

    main:
        // Merge fastq files of duplicate sample IDs
        MERGE_FASTQS(ch_samplesheet)
        merged_samplesheet = MERGE_FASTQS.out.merged_files

        // Download it when provided as an URL(s), multiple URLs are separated by whitespace
        def whitelist_param = params.bc_whitelist?.toString()?.trim()
        def whitelist_is_url = whitelist_param &&
            whitelist_param.tokenize().every { it ==~ /(?i)^(https?|ftp):\/\/.+/ }

        def resolved_whitelist
        if (whitelist_is_url) {
            DOWNLOAD_WHITELIST(whitelist_param)

            // Keep the whitespace-separated format
            resolved_whitelist = DOWNLOAD_WHITELIST.out.whitelist
                .map { files -> [files].flatten()*.toString().join(' ') }
        } else {
            resolved_whitelist = params.bc_whitelist
        }

        // Check for protocol and run appropriate pre-processing workflow
        if (params.protocol.startsWith('parse_biosciences')) {
            parse_workflow(merged_samplesheet)
            data_output_ch = parse_workflow.out.data_output
            bc_whitelist_ch  = resolved_whitelist

        } else if (params.protocol == 'bd_rhapsody') {
            bd_rhapsody_workflow(merged_samplesheet)
            data_output_ch = bd_rhapsody_workflow.out.data_output
            bc_whitelist_ch  = resolved_whitelist

        } else if (params.protocol.startsWith('10x') || params.protocol == 'oak_seq' || params.protocol == 'ultima_genomics') {
            tenx_genomics_workflow(merged_samplesheet)
            data_output_ch = tenx_genomics_workflow.out.data_output
            bc_whitelist_ch  = resolved_whitelist

        } else if (params.protocol == 'sciRNAseq3') {
            sciRNAseq3_nogather_workflow(merged_samplesheet)
            data_output_ch   = sciRNAseq3_nogather_workflow.out.data_output
            bc_whitelist_ch  = sciRNAseq3_nogather_workflow.out.bc_whitelist.map { tup -> tup*.toString().join(' ') }

        } else if (params.seqspec_file && file(params.seqspec_file).exists() && params.protocol == 'seqspec') {
            seqspec_workflow(merged_samplesheet)
            data_output_ch = seqspec_workflow.out.data_output
            bc_whitelist_ch  = resolved_whitelist

        } else {
            error """
            Invalid sequencing technology specified. Use one of the following parameters for 'protocol':
            - 'parse_biosciences_WT_mini'
            - 'parse_biosciences_WT'
            - 'bd_rhapsody'
            - '10xv1'
            - '10xv2'
            - '10xv3'
            - '10xv4'
            - 'oak_seq'
            - 'ultima_genomics'
            - 'sciRNAseq3'
            Or use 'seqspec' to specify a non-supported sequencing technique.
            """
        }

    emit:
        merged_samplesheet = merged_samplesheet
        data_output     = data_output_ch
        bc_whitelist    = bc_whitelist_ch
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
