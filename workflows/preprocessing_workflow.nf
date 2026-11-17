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
include { scalebio_workflow             } from '../subworkflows/local/pre-processing/preprocs_scalebio'

include { MERGE_FASTQS                  } from '../modules/local/custom/manipulate/merge_files/main'
include { CRAM_TO_FASTQ                 } from '../modules/local/tools/samtools/cram_to_fastq/main'
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
        - ScaleBio: Correcting the barcodes over all reads and splitting the samples by RT well
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow preprocessing_workflow {
    take:
        ch_samplesheet

    main:
        // Merge fastq files of duplicate sample IDs. ScaleBio samples share their library's
        // FASTQs, which SCALEBIO_DEMUX reads once per lane, so they are not merged per sample.
        def scalebio = params.protocol in ['scalebio_quantum_rna', 'scalebio_rna_v1_1']
        def ch_versions = channel.empty()

        // CRAM input. Illumina-style CRAMs (both reads, index reads in the BC/QT tags) become
        // bcl-convert-named FASTQs here, so every step after reads FASTQs. Ultima CRAMs hold
        // ScaleBio barcodes in tags and go to the ScaleBio demultiplexing as they are.
        def ch_input = ch_samplesheet
        if (params.cram_type == 'illumina') {
            def ch_split = ch_samplesheet.branch { row ->
                cram:  row[0].input_format == 'cram'
                fastq: true
            }
            // Each CRAM is converted once, however many samples read it (a ScaleBio library
            // holds many samples), and its FASTQs go back to every one of them
            def reference = params.cram_reference ? file(params.cram_reference, checkIfExists: true) : []
            def ch_cram_files = ch_split.cram
                .flatMap { _meta, crams, _bc, _indices, _input_file -> crams }
                .unique { cram -> cram.toString() }
                .map { cram -> [[id: cram.name, path: cram.toString()], cram] }
            CRAM_TO_FASTQ(ch_cram_files, reference)
            ch_versions = ch_versions.mix(CRAM_TO_FASTQ.out.versions.first())

            // A sample's CRAMs back together in their samplesheet order, the cDNA read in the
            // slot the protocol reads it from
            def cdna_read = params.seqtech_parameters[params.protocol]?.cdna_read
            def ch_fastqs_by_cram = CRAM_TO_FASTQ.out.fastqs
                .map { cram_meta, r1, r2, indices -> [cram_meta.path, [r1, r2, [indices].flatten()]] }
            def ch_converted = ch_split.cram
                .flatMap { meta, crams, _bc, _indices, input_file ->
                    crams.withIndex().collect { cram, k -> [cram.toString(), meta, input_file, k, crams.size()] }
                }
                .combine(ch_fastqs_by_cram, by: 0)
                .map { _path, meta, input_file, k, n_parts, fastqs -> [groupKey(meta.id, n_parts), meta, input_file, k, fastqs] }
                .groupTuple()
                .map { _key, metas, input_files, order, fastqs ->
                    def parts = [order, fastqs].transpose().sort { part -> part[0] }.collect { part -> part[1] }
                    def r1s = parts.collect { p -> p[0] }
                    def r2s = parts.collect { p -> p[1] }
                    def indices = parts.collectMany { p -> p[2] }
                    cdna_read == 'R1' ? [metas[0], r1s, r2s, indices, input_files[0]] : [metas[0], r2s, r1s, indices, input_files[0]]
                }
            ch_input = ch_split.fastq.mix(ch_converted)
        }

        merged_samplesheet = ch_input
        if (!scalebio) {
            MERGE_FASTQS(ch_input)
            merged_samplesheet = MERGE_FASTQS.out.merged_files
            ch_versions = MERGE_FASTQS.out.versions
        }

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
            // With the round-1 barcodes merged per well, the reads carry only T barcodes in
            // bc1: their list replaces the first whitelist file, bc2 and bc3 stay as given
            bc_whitelist_ch  = params.parse_collapse_bc1
                ? resolved_whitelist.map { files -> [files] }.combine(parse_workflow.out.bc1_whitelist)
                    .map { files, t_list -> [t_list] + files.drop(1) }.first()
                : resolved_whitelist
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

        } else if (scalebio) {
            if (!whitelist_tokens) {
                error("ScaleBio corrects every barcode against its list: set 'bc_whitelist' to the barcode lists, one per segment in segment order (see docs/PROTOCOLS_AND_WHITELISTS.md).")
            }
            scalebio_workflow(ch_input, resolved_whitelist)
            data_output_ch   = scalebio_workflow.out.data_output
            // The canonical barcodes SCALEBIO_DEMUX writes into the synthetic barcode read
            bc_whitelist_ch  = scalebio_workflow.out.bc_whitelist
            ch_versions = ch_versions.mix(scalebio_workflow.out.versions)

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
            - 'scalebio_quantum_rna'
            - 'scalebio_rna_v1_1'
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
