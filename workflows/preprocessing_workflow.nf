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
include { oak_seq_workflow              } from '../subworkflows/preprocs_oak_seq'
include { tenx_genomics_workflow        } from '../subworkflows/preprocs_10xv3'
include { sci_rna_seq3_workflow         } from '../subworkflows/preprocs_sci_rna_seq3'
include { seqspec_workflow              } from '../subworkflows/preprocs_seqspec'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW TO RUN PRE-PROCESSING
        Sequencing-specific pre-processing of the data:
        - Parse Bioscience: Demultiplexing using groups of wells and mapping using split-pipe
        - BD Rhapsody: Removing variable bases and mapping using BD rhapsody pipeline
        - OAK seq & 10xv3: Mapping using CellRanger
        - Sci-RNA-seq3: Pre-processing based on the sci-rocket pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow preprocessing_workflow {
    take:
        ch_samplesheet

    main:
        if (params.protocol == 'parse_biosciences') {     
            data_output = parse_workflow(ch_samplesheet)

        } else if (params.protocol == 'bd_rhapsody') {
            data_output = bd_rhapsody_workflow(ch_samplesheet)

        } else if (params.protocol == 'oak_seq') {
            data_output = oak_seq_workflow(ch_samplesheet)

        } else if (params.protocol == '10xv3') {
            data_output = tenx_genomics_workflow(ch_samplesheet)

        } else if (params.protocol == 'sci_rna_seq3') {
            data_output = sci_rna_seq3_workflow(ch_samplesheet)

        } else if (params.seqspec_file && file(params.seqspec_file).exists() && params.protocol == 'seqspec') {     
            data_output = seqspec_workflow(ch_samplesheet)
            
        } else {
            error """
            Invalid sequencing technology specified. Use one of the following parameters for 'protocol': 
            - 'parse_biosciences' 
            - 'bd_rhapsody' 
            - 'oak_seq' 
            - '10xv3'
            - 'sci_rna_seq3'
            Or use 'seqspec' to specify a non-supported sequencing technique.
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