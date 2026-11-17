/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Barcode check
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BARCODE_CHECK runs bin/barcode_check.sh on every sample before mapping, and writes
    a verdict (pass, fail or skipped) to a one-row TSV. This class holds the Nextflow
    side: reading that row, and turning a failed sample into a stopped pipeline or a
    dropped sample, as params.mapping_rate_action says for the mapping-rate check.

    The verdict is a task output, not an exit code or a shared file: the task itself
    succeeded, and its numbers reach the dashboard also for a sample that failed.

    See "Barcode check" in docs/CONFIGURATION_PARAMETERS.md.
----------------------------------------------------------------------------------------
*/

class BarcodeCheck {

    /** Samples dropped in this run, for the summary at the end. */
    private static final List<Map> DROPPED = Collections.synchronizedList([])

    static boolean enabled(Object skip, Object minMatch) {
        return !skip && (minMatch ?: 0) as double > 0
    }

    /**
     * The verdict rows of a barcode-check TSV, each keyed on the header: one row in a
     * sample's own file, one per sample in pipeline_info/barcode_check.tsv.
     */
    static List<Map> readAll(Object tsv) {
        def lines = tsv.text.readLines().findAll { it.trim() }
        if (lines.size() < 2) {
            return []
        }
        def header = lines[0].split('\t', -1) as List
        return lines.drop(1).collect { line ->
            [header, line.split('\t', -1) as List].transpose().collectEntries { key, value -> [(key): value] }
        }
    }

    /** The verdict of one sample's file. */
    static Map read(Object tsv) {
        def rows = readAll(tsv)
        return rows ? rows[0] : [verdict: 'skipped', note: 'no verdict written']
    }

    static boolean failed(Map row) {
        return row.verdict == 'fail'
    }

    /** One sample of the message: its numbers, and what the reads show. */
    static String describe(Map row) {
        def line = "  - ${row.sample}: ${row.pct_match}% of the first ${row.reads} barcode reads are on the whitelist"
        if (row.pct_per_segment?.toString()?.contains(',')) {
            line += " (per segment: ${row.pct_per_segment}%)"
        }
        return row.note ? "${line}\n      ${row.note}" : line
    }

    /** The block shown when samples fail; 'aborted' says which action was taken. */
    static String failureMessage(List<Map> rows, Object minMatch, boolean aborted) {
        return [
            '',
            '====================== BARCODE CHECK FAILED =======================',
            "${rows.size()} sample(s) have barcodes that do not match the whitelist (threshold ${minMatch}%):",
            *rows.collect { row -> describe(row) },
            '',
            aborted
                ? 'The pipeline was stopped before mapping (mapping_rate_action = "abort").'
                : 'These samples were dropped; all other samples are processed (mapping_rate_action = "skip_sample").',
            '',
            "Likely causes: the whitelist of another kit version, a wrong 'protocol', whitelist files",
            'given in the wrong order, or the cDNA and barcode reads swapped in the samplesheet.',
            '',
            'To disable this check: --skip_barcode_check true (or lower --min_barcode_match).',
            '==================================================================='
        ].join('\n')
    }

    /** Remember a dropped sample for summarise(). */
    static void recordDropped(Map row) {
        DROPPED << row
    }

    /**
     * Repeat the dropped samples at the end of the run, where the warning given when
     * each was dropped has long scrolled away. Called from workflow.onComplete, so it
     * must never throw.
     */
    static void summarise(Map params, log) {
        try {
            if (DROPPED) {
                def rows = DROPPED.toList().sort { it.sample }
                log.warn(failureMessage(rows, params.min_barcode_match, false) +
                    "\nDetails: ${params.outdir}/pipeline_info/barcode_check.tsv")
            }
        } catch (Exception e) {
            log.warn("Barcode check summary could not be written: ${e.message}")
        }
    }
}
