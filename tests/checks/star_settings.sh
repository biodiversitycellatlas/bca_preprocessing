#!/usr/bin/env bash
# description: Check the strand default and the merging of star_extraargs (lib/StarSettings.groovy)
#
# Two STAR settings are resolved from more than one source:
#   - the strand: star_soloStrand when set; otherwise Reverse for a sample read from an
#     Ultima CRAM and Forward for every other sample. A wrong strand is silent: the run
#     finishes with few counts.
#   - the extra arguments: the protocol's, then star_extraargs. STAR refuses an option
#     given twice, so a custom option has to replace the protocol's option of the same
#     name, with all its values, and leave the protocol's other options in place.
#
# The config case reads conf/seqtech_parameters.config: the ScaleBio protocols must carry
# the settings docs/PROTOCOLS_AND_WHITELISTS.md says they share with ScaleRna.
#
# Needs nextflow, and is skipped without it.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

KEEP=0

usage() {
    cat <<EOF
Usage: tests/run_tests.sh star_settings [-- OPTIONS]
       tests/checks/star_settings.sh [OPTIONS]

Check lib/StarSettings.groovy and the ScaleBio STAR defaults.

Options:
  --keep          Keep the generated project and outputs for inspection.
  -h, --help      Show this message.

Cases:
  strand          star_soloStrand wins; unset, only Ultima CRAM samples are Reverse
  merge           custom options are added to the protocol's and replace the same option
  merge_invalid   star_extraargs that starts with a value is refused
  config          the ScaleBio protocols: PropUnique, 6 loci, no minimum score
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --keep)    KEEP=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) log_error "unknown argument: $1"; usage >&2; exit 2 ;;
    esac
done

log_header "star_settings"

if ! have_cmd nextflow; then
    record SKIP "nextflow" "nextflow not found"
    finish_check
    exit 0
fi

WORKDIR="$BCA_TEST_LOGDIR/star_settings"
rm -rf "$WORKDIR"
mkdir -p "$WORKDIR"
cleanup() { [[ "$KEEP" -eq 1 ]] || rm -rf "$WORKDIR"; }
trap cleanup EXIT

ln -s "$PROJECT_ROOT/lib" "$WORKDIR/lib"
ln -s "$PROJECT_ROOT/conf" "$WORKDIR/conf"
cat >"$WORKDIR/nextflow.config" <<'CFG'
includeConfig 'conf/seqtech_parameters.config'
CFG
cat >"$WORKDIR/main.nf" <<'NF'
def bd = '--outFilterScoreMinOverLread 0 --clip3pAdapterSeq AAAAAAAA --outFilterMatchNmin 25'

// [case, got, want]
def results = [
    ['strand', StarSettings.strand(null, 'illumina', 'fastq'),       'Forward'],
    ['strand', StarSettings.strand(null, 'illumina', 'cram'),        'Forward'],
    ['strand', StarSettings.strand(null, 'ultima', 'cram'),          'Reverse'],
    ['strand', StarSettings.strand(null, 'ultima', 'fastq'),         'Forward'],
    ['strand', StarSettings.strand(null, 'ultima', null),            'Forward'],
    ['strand', StarSettings.strand('', 'ultima', 'cram'),            'Reverse'],
    ['strand', StarSettings.strand('Forward', 'ultima', 'cram'),     'Forward'],
    ['strand', StarSettings.strand('Unstranded', 'ultima', 'cram'),  'Unstranded'],
    ['strand', StarSettings.strand('Reverse', 'illumina', 'fastq'),  'Reverse'],

    ['merge', StarSettings.mergeArgs(null, null),                    ''],
    ['merge', StarSettings.mergeArgs(bd, null),                      bd],
    ['merge', StarSettings.mergeArgs(bd, ''),                        bd],
    ['merge', StarSettings.mergeArgs(null, '--outSJtype None'),      '--outSJtype None'],
    ['merge', StarSettings.mergeArgs(bd, '--outSJtype None'),        bd + ' --outSJtype None'],
    ['merge', StarSettings.mergeArgs(bd, '--outFilterMatchNmin 30'),
        '--outFilterScoreMinOverLread 0 --clip3pAdapterSeq AAAAAAAA --outFilterMatchNmin 30'],
    ['merge', StarSettings.mergeArgs('--outFilterMultimapNmax 6', '--outSJtype None  --outFilterMultimapNmax 10'),
        '--outFilterMultimapNmax 10 --outSJtype None'],
    // an option with several values is replaced whole, and a negative number is a value
    ['merge', StarSettings.mergeArgs('--outSJfilterOverhangMin 30 12 12 12 --alignIntronMax 1000',
                                     '--outSJfilterOverhangMin -1 20 --clip3pNbases 3'),
        '--outSJfilterOverhangMin -1 20 --alignIntronMax 1000 --clip3pNbases 3'],
]
results.each { name, got, want ->
    println "${got == want ? 'OK' : 'BAD'}\t${name}\t'${got}'\twant '${want}'"
}

try {
    StarSettings.mergeArgs('--outFilterMultimapNmax 6', '10 --outSJtype None')
    println "BAD\tmerge_invalid\tno error"
} catch (IllegalArgumentException e) {
    println "${e.message.contains('star_extraargs') ? 'OK' : 'BAD'}\tmerge_invalid\t${e.message}"
}

['scalebio_quantum_rna', 'scalebio_rna_v1_1'].each { protocol ->
    def seqtech = params.seqtech_parameters[protocol]
    def got = [seqtech.star_soloMultiMappers, seqtech.star_outFilterScoreMin, StarSettings.mergeArgs(seqtech.star_extraargs, null)]
    def want = ['PropUnique', 0, '--outFilterMultimapNmax 6']
    println "${got == want ? 'OK' : 'BAD'}\tconfig\t${protocol} ${got}\twant ${want}"
}
NF

LOG="$BCA_TEST_LOGDIR/star_settings.log"
if (cd "$WORKDIR" && run_logged "$LOG" nextflow -q run main.nf); then
    for case_name in strand merge merge_invalid config; do
        expected=$(grep -c "^\s*\['$case_name'," "$WORKDIR/main.nf" || true)
        [[ "$case_name" == "merge_invalid" ]] && expected=1
        [[ "$case_name" == "config" ]] && expected=2
        ok=$(grep -c -P "^OK\t$case_name\t" "$LOG" || true)
        bad=$(grep -P "^BAD\t$case_name\t" "$LOG" | head -2 | tr '\n' ' ' || true)
        if [[ -n "$bad" ]]; then
            record FAIL "$case_name" "$bad"
        elif [[ "$ok" -eq "$expected" ]]; then
            record PASS "$case_name" "$ok of $expected as expected"
        else
            record FAIL "$case_name" "expected $expected results, got $ok; see star_settings.log"
        fi
    done
else
    record FAIL "run" "nextflow run failed, see star_settings.log"
fi

finish_check
