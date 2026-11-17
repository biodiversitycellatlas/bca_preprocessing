//
// Subworkflow with functionality specific to the workflow 'preprocessing_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SCALEBIO_PREPARE_BARCODES } from '../../../modules/local/custom/demultiplex/scalebio_prepare_barcodes/main'
include { SCALEBIO_DEMUX            } from '../../../modules/local/custom/demultiplex/scalebio_demultiplex/main'
include { SCALEBIO_GATHER           } from '../../../modules/local/custom/demultiplex/scalebio_gather/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO RUN PRE-PROCESSING FOR SCALEBIO DATA
        A ScaleBio library holds many samples, told apart by their RT wells (and PCR
        wells). Samples read from the same FASTQs form one library, which is read once
        per lane: SCALEBIO_DEMUX corrects the barcodes spread over R1/R2/I1/I2, trims the
        cDNA and writes each sample's reads with a synthetic barcode read (the corrected
        barcodes, then the UMI). SCALEBIO_GATHER joins a sample's lanes for the mappers.
        Independent implementation (bin/scalebio_demux.py): no ScaleBio code is used.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow scalebio_workflow {
    take:
        ch_samplesheet      // [meta, [cDNA FASTQs], [barcode FASTQs], [index FASTQs], samplesheet], not merged;
                            // Ultima CRAMs: [meta, [CRAMs], [], [], samplesheet]
        ch_whitelist        // value channel: the barcode lists, one per segment in segment order

    main:
        def layout = params.seqtech_parameters[params.protocol].scalebio_layout
        def ch_versions = channel.empty()

        SCALEBIO_PREPARE_BARCODES(ch_whitelist, layout)
        ch_versions = ch_versions.mix(SCALEBIO_PREPARE_BARCODES.out.versions)

        // One task per lane of each library
        def ch_chunks = ch_samplesheet
            .toList()
            .flatMap { rows -> scalebioChunks(rows, layout) }

        def reference = params.cram_reference ? file(params.cram_reference, checkIfExists: true) : []
        SCALEBIO_DEMUX(ch_chunks, SCALEBIO_PREPARE_BARCODES.out.barcodes, reference)
        ch_versions = ch_versions.mix(SCALEBIO_DEMUX.out.versions.first())

        // Each sample's reads from every chunk it was read in, with those chunks' metrics
        def ch_parts = SCALEBIO_DEMUX.out.reads
            .flatMap { library, chunk, cdna, bc, metrics -> scalebioSampleParts(library, chunk, cdna, bc, metrics) }
            .groupTuple()

        def ch_gather_input = ch_samplesheet
            .map { meta, _cdna, _bc, _indices, input_file -> [meta.id, meta, input_file] }
            // No failOnMismatch: a sample without reads already fails in scalebioSampleParts,
            // and after any failed task it would add a misleading 'Join mismatch' error
            .join(ch_parts)
            .map { _id, meta, input_file, cdna_parts, bc_parts, metrics -> [meta, cdna_parts, bc_parts, metrics, input_file] }

        SCALEBIO_GATHER(ch_gather_input)
        ch_versions = ch_versions.mix(SCALEBIO_GATHER.out.versions.first())

        // The mappers' tuple: the synthetic barcode read takes the barcode slot, no index reads
        def ch_data_output = SCALEBIO_GATHER.out.reads
            .map { meta, cdna, bc, input_file -> [meta, cdna, bc, [], input_file] }

        // The canonical barcodes of every segment, in segment order: the mappers' whitelist.
        // Without one (QuantumScale: too many combinations for STARsolo), the corrected
        // barcodes are taken as they are ('--soloCBwhitelist None').
        def ch_bc_whitelist = params.seqtech_parameters[params.protocol].scalebio_mapping_whitelist
            ? SCALEBIO_PREPARE_BARCODES.out.whitelists.map { files -> [files].flatten().sort { f -> f.name } }
            : channel.value([])

    emit:
        data_output  = ch_data_output
        bc_whitelist = ch_bc_whitelist
        demux_stats  = SCALEBIO_GATHER.out.stats
        versions     = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Which index read a FASTQ holds, from its name (bcl-convert '_I1_'/'_I2_', or R3/R4, index1/2)
def scalebioIndexRead(String name) {
    if (name ==~ /.*(^|[_.])(I1|R3|index1)([_.].*|$)/) {
        return 1
    }
    if (name ==~ /.*(^|[_.])(I2|R4|index2)([_.].*|$)/) {
        return 2
    }
    return 0
}

/* 
    A FASTQ's lane key: its folder and name with the read token blanked, so R1/R2/I1/I2 of one
    lane share it. Null for names without a bcl-convert read token.
*/
def scalebioLaneKey(Path file) {
    def m = file.name =~ /^(.*)_(R1|R2|R3|R4|I1|I2)(_\d{3}\.f(ast)?q(\.gz)?)$/
    return m.matches() ? "${file.parent}/${m.group(1)}_X${m.group(3)}" : null
}

// Files sorted by full path (lanes of different runs may share a name), without sorting the input in place
def scalebioByPath(List files) {
    return files.sort(false) { f -> f.toString() }
}

// A file holding reads as SAM/BAM/CRAM (Ultima CRAMs reach the demultiplexing as they are)
def scalebioIsCram(Object file) {
    return file.toString() ==~ /(?i).*\.(cram|bam|sam)$/
}

/* 
    The SCALEBIO_DEMUX tasks: samples read from the same FASTQs (or Ultima CRAMs) form one
    library, read once per lane (or CRAM). Each task: [library, chunk, cDNA, barcode read,
    I1 (or []), I2 (or [])]; for an Ultima CRAM: [library, chunk, CRAM, [], [], []]
*/
def scalebioChunks(List rows, String layout) {
    def needs_i1 = layout == 'quantum'
    def libraries = rows.groupBy { row -> [row[1], row[2], row[3]].flatten().collect { f -> f.toString() }.sort() }
    def chunks = []
    libraries.values().eachWithIndex { lib_rows, n ->
        def crams = lib_rows[0][1].findAll { f -> scalebioIsCram(f) }
        def library = [
            id     : "lib${n + 1}",
            ultima : !crams.isEmpty(),
            samples: lib_rows.collect { row -> [id: row[0].id, rt: row[0].rt ?: '', pcr: row[0].p5 ?: ''] },
        ]
        def names = library.samples.collect { s -> s.id }.join(', ')
        if (library.ultima) {
            if (crams.size() != lib_rows[0][1].size() || lib_rows[0][2] || lib_rows[0][3]) {
                error("ScaleBio sample(s) ${names}: Ultima CRAMs hold every read; give them in 'cram' only, without FASTQs.")
            }
            scalebioByPath(crams).eachWithIndex { cram, k ->
                def lane_token = cram.name.find(/_L\d{3}[_.]/)
                chunks << [library, "${library.id}_${String.format('%02d', k + 1)}" + (lane_token ? lane_token[0..-2] : ''), cram, [], [], []]
            }
            return
        }
        def cdna = scalebioByPath(lib_rows[0][1])
        def bc   = scalebioByPath(lib_rows[0][2])
        def i1   = scalebioByPath(lib_rows[0][3].findAll { f -> scalebioIndexRead(f.name) == 1 })
        def i2   = scalebioByPath(lib_rows[0][3].findAll { f -> scalebioIndexRead(f.name) == 2 })
        // v1.1 reads only the i5 index: a lone index read (e.g. a single-index CRAM's, which
        // CRAM_TO_FASTQ writes as I1) is that one
        if (!needs_i1 && !i2 && i1) {
            i2 = i1
            i1 = []
        }
        def counts = [cDNA: cdna.size(), barcode: bc.size(), I2: i2.size()] + (needs_i1 ? [I1: i1.size()] : [:])
        if (counts.values().toSet().size() != 1) {
            error("ScaleBio sample(s) ${names}: every read needs one file per lane, found ${counts.collect { k, v -> "${v} ${k}" }.join(', ')}. " +
                "The index FASTQs go in 'fastq_indices' (I2${needs_i1 ? ' and I1' : ''}, named with _I1_/_I2_).")
        }
        cdna.indices.each { k ->
            def lane = [cdna[k], bc[k], i2[k]] + (needs_i1 ? [i1[k]] : [])
            def keys = lane.collect { f -> scalebioLaneKey(f) }
            if (!keys.contains(null) && keys.toSet().size() != 1) {
                error("ScaleBio sample(s) ${names}: the read files do not pair up by lane: ${lane.collect { f -> f.name }.join(', ')}.")
            }
            def lane_token = cdna[k].name.find(/_L\d{3}_/)
            def chunk = "${library.id}_${String.format('%02d', k + 1)}" + (lane_token ? lane_token[0..-2] : '')
            chunks << [library, chunk, cdna[k], bc[k], needs_i1 ? i1[k] : [], i2[k]]
        }
    }
    return chunks
}

// One SCALEBIO_DEMUX task's output split per sample: [sample id, cDNA part, barcode part, metrics]
def scalebioSampleParts(Map library, String chunk, Object cdna, Object bc, Path metrics) {
    def cdna_files = [cdna].flatten()
    def bc_files   = [bc].flatten()
    return library.samples.collect { s ->
        def c = cdna_files.find { f -> f.name == "${s.id}.${chunk}_cDNA.fastq.gz" }
        def b = bc_files.find { f -> f.name == "${s.id}.${chunk}_BC_UMI.fastq.gz" }
        if (!c || !b) {
            error("SCALEBIO_DEMUX wrote no reads file for sample '${s.id}' in chunk ${chunk}.")
        }
        [s.id, c, b, metrics]
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
