//
// Workflow with functionality specific to 'main.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { bd_rhapsody_workflow          } from '../subworkflows/preprocs_bd_rhapsody'
include { parse_workflow                } from '../subworkflows/preprocs_parse_biosciences'
include { tenx_genomics_workflow        } from '../subworkflows/preprocs_10xv3'
include { sciRNAseq3_workflow         } from '../subworkflows/preprocs_sciRNAseq3'
include { seqspec_workflow              } from '../subworkflows/preprocs_seqspec'

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
        if (params.protocol == 'parse_biosciences') {     
            data_output_ch = parse_workflow(ch_samplesheet)
            bc_whitelist_ch  = params.seqtech_parameters[params.protocol].bc_whitelist

        } else if (params.protocol == 'bd_rhapsody') {
            data_output_ch = bd_rhapsody_workflow(ch_samplesheet)
            bc_whitelist_ch  = Channel.value( params.seqtech_parameters[params.protocol].bc_whitelist )

        } else if (params.protocol == '10xv3' || params.protocol == 'oak_seq' || params.protocol == 'ultima_genomics') {
            tenx_genomics_workflow(ch_samplesheet)
            data_output_ch = tenx_genomics_workflow.out.data_output
            bc_whitelist_ch  = tenx_genomics_workflow.out.bc_whitelist

        } else if (params.protocol == 'sciRNAseq3') {
            sciRNAseq3_workflow(ch_samplesheet)
            data_output_ch   = sciRNAseq3_workflow.out.data_output
            bc_whitelist_ch  = sciRNAseq3_workflow.out.bc_whitelist.map { tup -> tup*.toString().join(' ') }

        } else if (params.seqspec_file && file(params.seqspec_file).exists() && params.protocol == 'seqspec') {     
            data_output_ch = seqspec_workflow(ch_samplesheet)
            bc_whitelist_ch  = Channel.value( params.seqtech_parameters[params.protocol].bc_whitelist )
            
        } else {
            error """
            Invalid sequencing technology specified. Use one of the following parameters for 'protocol': 
            - 'parse_biosciences' 
            - 'bd_rhapsody'
            - '10xv3' 
            - 'oak_seq'
            - 'ultima_genomics' 
            - 'sciRNAseq3'
            Or use 'seqspec' to specify a non-supported sequencing technique.
            """
        }

    emit:
        data_output     = data_output_ch
        bc_whitelist    = bc_whitelist_ch
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/