//
// Workflow with functionality specific to 'main.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { bd_rhapsody_workflow } from '../subworkflows/preprocs_bd_rhapsody.nf'
include { parse_workflow } from '../subworkflows/preprocs_parse.nf'
include { oak_seq_workflow } from '../subworkflows/preprocs_oak_seq.nf'
include { tenx_genomics_workflow } from '../subworkflows/preprocs_10x_genomics.nf'
include { seqspec_workflow } from '../subworkflows/preprocs_seqspec.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW TO RUN PRE-PROCESSING
        Sequencing-specific pre-processing of the data:
        - Parse Bioscience: Demultiplexing using groups of wells and mapping using split-pipe
        - BD Rhapsody: Removing variable bases and mapping using BD rhapsody pipeline
        - OAK seq & 10xv3: Mapping using CellRanger
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow preprocessing_workflow {
    take:
        sample_ids

    main:
        if (params.protocol == 'parse_biosciences') {     
            data_output = parse_workflow(sample_ids)

        } else if (params.protocol == 'bd_rhapsody') {
            data_output = bd_rhapsody_workflow(sample_ids)

        } else if (params.protocol == 'oak_seq') {
            data_output = oak_seq_workflow(sample_ids)

        } else if (params.protocol == '10xv3') {
            data_output = tenx_genomics_workflow(sample_ids)

        } else if (params.seqspec_file && file(params.seqspec_file).exists() && params.protocol == 'seqspec') {     
            data_output = seqspec_workflow(sample_ids)
            
        } else {
            error """
            Invalid sequencing technology specified. Use one of the following parameters for 'protocol': 
            - 'parse_biosciences' 
            - 'bd_rhapsody' 
            - 'oak_seq' 
            - '10xv3'
            Or use 'seqspec' to specify a custom workflow.
            """
        }

    emit:
        data_output
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/