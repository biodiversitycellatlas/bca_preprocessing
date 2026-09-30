/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Mapping-rate check
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    STARSOLO_ALIGN and ALEVIN_FRY run their mapper through bin/mapping_rate_guard.sh,
    which exits with EXIT_CODE once a sample is known to map below
    params.min_mapping_rate. This class holds the Nextflow side of that contract:

    - flagDir():       where the mappers of one sample share their verdict. It lives
                       under workDir, the one location every executor shares between
                       tasks, and is keyed on the session id, which -resume keeps, so
                       the task hash (and with it the cache) is stable across resumes.
    - reset():         clears that directory at startup, so a verdict from an earlier
                       attempt of this session cannot cancel a rerun with new settings.
    - errorStrategy(): 'terminate' (abort the pipeline) or 'ignore' (drop the sample)
                       for EXIT_CODE, the default rule from conf/base.config otherwise.
    - summarise():     reports the verdicts at the end of the run and writes them to
                       pipeline_info/mapping_rate_check.tsv.

    See "Mapping-rate check" in docs/CONFIGURATION_PARAMETERS.md.
----------------------------------------------------------------------------------------
*/

import java.nio.file.Files
import java.nio.file.Path

class MappingRateCheck {

    static final int EXIT_CODE = 42

    static final List<Integer> RETRY_CODES = ((130..145) + 104 + 175 + 255) as List<Integer>

    static final List<String> VALID_ACTIONS = ['abort', 'skip_sample']

    static boolean enabled(Map params) {
        return !params.skip_mapping_rate_check && (params.min_mapping_rate ?: 0) as double > 0
    }

    static Path flagDir(workflow) {
        return workflow.workDir.resolve('.mapping_rate_check').resolve(workflow.sessionId.toString())
    }

    static void reset(workflow) {
        def dir = flagDir(workflow)
        if (Files.isDirectory(dir)) {
            dir.toFile().deleteDir()
        }
    }

    static String errorStrategy(task, Map params) {
        if (task.exitStatus == EXIT_CODE) {
            return params.mapping_rate_action == 'skip_sample' ? 'ignore' : 'terminate'
        }
        return task.exitStatus in RETRY_CODES ? 'retry' : 'finish'
    }

    /**
     * Report every verdict at the end of the run. Returns true if any sample failed.
     * Called from workflow.onComplete, so it must never throw.
     */
    static boolean summarise(workflow, Map params, log) {
        try {
            def dir = flagDir(workflow)
            if (!Files.isDirectory(dir)) {
                return false
            }
            def rows = dir.toFile().listFiles()
                .findAll { it.name.endsWith('.verdict') }
                .sort { it.name }
                .collect { it.text.trim().split('\t') as List }
                .findAll { it.size() >= 8 }
            if (!rows) {
                return false
            }

            def out = new File("${params.outdir}/pipeline_info")
            out.mkdirs()
            new File(out, 'mapping_rate_check.tsv').text =
                (['sample', 'verdict', 'mapper', 'task', 'stage', 'reads', 'percent_mapped', 'threshold'].join('\t') + '\n') +
                rows.collect { it.join('\t') }.join('\n') + '\n'

            def failed = rows.findAll { it[1] == 'fail' }
            if (!failed) {
                return false
            }

            def aborted = params.mapping_rate_action != 'skip_sample'
            def lines = failed.collect { r ->
                def mapper = r[2] == 'star' ? 'STARsolo, uniquely mapped' : 'alevin-fry, mapped'
                "  - ${r[0]}: ${r[6]}% (${mapper}, ${r[4]} check after ${r[5]} reads)"
            }
            def msg = [
                '',
                '==================== MAPPING RATE CHECK FAILED ====================',
                "${failed.size()} sample(s) mapped below the ${params.min_mapping_rate}% threshold:",
                *lines,
                '',
                aborted
                    ? 'The pipeline was aborted (mapping_rate_action = "abort").'
                    : 'These samples were dropped; all other samples were processed (mapping_rate_action = "skip_sample").',
                '',
                'Recommendation: classify the unmapped reads with Kraken2 to find out what they are,',
                'e.g. rerun these samples with',
                '    --skip_mapping_rate_check true --perform_kraken true --star_generateBAM true',
                'or run kraken2 directly on a subsample of the cDNA FASTQ.',
                '',
                'To disable this check: --skip_mapping_rate_check true (or lower --min_mapping_rate).',
                "Details: ${params.outdir}/pipeline_info/mapping_rate_check.tsv",
                '==================================================================='
            ].join('\n')
            if (aborted) { log.error(msg) } else { log.warn(msg) }
            return true
        } catch (Exception e) {
            log.warn("Mapping-rate check summary could not be written: ${e.message}")
            return false
        }
    }
}
