/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Analytical runs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Every mapping of a sample is an analytical run, named <sample><suffix>. Which runs
    exist follows from mapping_software, perform_geneext, geneext_downstream and
    run_method. mapping_workflow.nf maps them, restage_mapping.nf reads them back for
    post_mapping, and the parameter checks validate them, all from this one plan.

    - star_standard / alevin_standard: the runs against the standard annotation. With
      perform_geneext they only exist when geneext_downstream = 'both'.
    - geneext_input: the full-data STARsolo pass whose BAMs GeneExt reads, always named
      '_starsolo'. 'reuse_standard' when that is the standard STARsolo run itself,
      'own_run' when it is mapped only for GeneExt (stripped to a BAM with
      geneext_bam_only) and is neither reported nor carried downstream.
    - star_geneext / alevin_geneext: the re-mappings against the extended annotation.

    See "GeneExt" in docs/CONFIGURATION_PARAMETERS.md.
----------------------------------------------------------------------------------------
*/

class AnalyticalRuns {

    static final List<String> VALID_GENEEXT_DOWNSTREAM = ['geneext_only', 'both']

    static final String GENEEXT_INPUT_SUFFIX = '_starsolo'

    static Map plan(Map params) {
        def software   = params.mapping_software
        def star_full  = software in ['starsolo', 'both', 'alevin_starsolo']
        def star_sub   = software == 'alevin_subsampled_starsolo'
        def alevin     = software in ['alevin', 'both', 'alevin_starsolo', 'alevin_subsampled_starsolo']

        // GeneExt alone: only the pass it reads is mapped
        if (params.run_method == 'geneext_only') {
            return [
                star_standard   : null,
                alevin_standard : false,
                geneext_input   : 'own_run',
                star_geneext    : null,
                alevin_geneext  : false,
            ]
        }

        def geneext  = params.perform_geneext as boolean
        def standard = !geneext || params.geneext_downstream == 'both'

        def star_standard = !standard ? null
            : star_sub  ? '_subsampled_starsolo'
            : star_full ? '_starsolo'
            : null

        return [
            star_standard   : star_standard,
            alevin_standard : standard && alevin,
            geneext_input   : !geneext ? null
                : star_standard == GENEEXT_INPUT_SUFFIX ? 'reuse_standard'
                : 'own_run',
            star_geneext    : !geneext ? null
                : star_sub  ? '_geneext_subsampled_starsolo'
                : star_full ? '_geneext_starsolo'
                : null,
            alevin_geneext  : geneext && alevin,
        ]
    }

    /** The STARsolo runs that are reported and carried downstream, standard first. */
    static List<String> starSuffixes(Map plan) {
        return [plan.star_standard, plan.star_geneext].findAll { sfx -> sfx }
    }

    /** The alevin-fry runs that are reported and carried downstream, standard first. */
    static List<String> alevinSuffixes(Map plan) {
        return (plan.alevin_standard ? ['_alevinfry'] : []) + (plan.alevin_geneext ? ['_geneext_alevinfry'] : [])
    }
}
