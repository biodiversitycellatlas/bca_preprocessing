/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Vendor pipelines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Which vendor pipeline a protocol has, which settings enable it, and the startup
    checks on those settings. A vendor pipeline only runs when it is enabled, so a
    missing or misspelt setting used to skip it without a word:

    - 'external_pipeline_only' with no vendor pipeline enabled, or a protocol that has
      none, pre-processed the reads, ran nothing else and reported success;
    - a 'splitpipe_installation' that does not exist skipped split-pipe silently;
    - 'splitpipe_installation' without 'splitpipe_conda_env' failed only at the task;
    - a vendor setting given for another protocol was ignored;
    - a vendor pipeline run under a profile it cannot run under (Cell Ranger without a
      container engine, split-pipe or BD Rhapsody with one) failed only at the task,
      or ran whatever was on the host PATH.

    check() returns the errors and warnings; PIPELINE_INITIALISATION reports them.
    Vendor pipelines are not rerun with 'post_mapping', so there their settings only
    give a warning. An unknown protocol is left to the error in preprocessing_workflow.

    See "External pipelines only" in docs/usage.md.
----------------------------------------------------------------------------------------
*/

import java.nio.file.Files
import java.nio.file.Path

class VendorPipelines {

    // The protocols each vendor pipeline runs for, and the settings that belong to it
    static final Map<String, Map> PIPELINES = [
        cellranger: [
            name     : 'Cell Ranger',
            protocols: ['10xv2', '10xv3', '10xv4', 'oak_v1'],
            settings : ['perform_cellranger'],
        ],
        splitpipe: [
            name     : 'split-pipe',
            protocols: ['parse_biosciences_WT_mini', 'parse_biosciences_WT'],
            settings : ['splitpipe_installation', 'splitpipe_conda_env', 'splitpipe_demultiplex_script', 'splitpipe_kit'],
        ],
        rhapsody: [
            name     : 'the BD Rhapsody pipeline',
            protocols: ['bd_rhapsody_v1', 'bd_rhapsody_enhancedbeads'],
            settings : ['rhapsody_installation'],
        ],
    ]

    // Protocols the pipeline supports that have no vendor pipeline
    static final List<String> NO_VENDOR_PROTOCOLS = [
        'marsseq_v1', 'marsseq_v2', 'sciRNAseq3', 'scalebio_quantum_rna', 'scalebio_rna_v1_1', 'seqspec'
    ]

    // A flag is given when true, a path or name when it is not blank
    static boolean given(value) {
        if (value instanceof Boolean) {
            return value
        }
        return value != null && value.toString().trim() as boolean
    }

    static boolean pathExists(value) {
        def path = nextflow.Nextflow.file(value.toString().trim())
        return path instanceof Path && Files.exists(path)
    }

    // containerEngine is workflow.containerEngine: null when no container engine is enabled
    static Map check(Map params, String containerEngine = null) {
        def errors   = []
        def warnings = []
        def protocol = params.protocol?.toString()
        def method   = params.run_method?.toString()

        def givenSettings = { String key ->
            PIPELINES[key].settings.findAll { name -> given(params[name]) }
        }

        // Nothing is rerun before the mapping results are read back in
        if (method == 'post_mapping') {
            def ignored = PIPELINES.keySet().collectMany { key -> givenSettings.call(key) }
            if (ignored) {
                warnings << (
                    "${ignored.collect { name -> "'${name}'" }.join(', ')} ${ignored.size() == 1 ? 'is' : 'are'} ignored with 'run_method' = 'post_mapping':\n" +
                    "vendor pipelines are not rerun, only everything after mapping is."
                ).toString()
            }
            return [errors: errors, warnings: warnings]
        }

        def own = PIPELINES.find { _key, pipeline -> protocol in pipeline.protocols }?.key
        if (own == null && !(protocol in NO_VENDOR_PROTOCOLS)) {
            return [errors: errors, warnings: warnings]
        }

        // Settings of another protocol's vendor pipeline
        PIPELINES.each { key, pipeline ->
            if (key == own) {
                return
            }
            def stray = givenSettings.call(key)
            if (stray) {
                warnings << (
                    "${stray.collect { name -> "'${name}'" }.join(', ')} ${stray.size() == 1 ? 'is' : 'are'} ignored: ${stray.size() == 1 ? 'it belongs' : 'they belong'} to ${pipeline.name},\n" +
                    "which runs for ${pipeline.protocols.join(', ')}, not for protocol '${protocol}'."
                ).toString()
            }
        }

        // The protocol's own vendor pipeline: is it enabled, and can it run?
        def enabled = false
        if (own == 'cellranger') {
            enabled = given(params.perform_cellranger)
        }
        if (own == 'splitpipe' && given(params.splitpipe_installation)) {
            enabled = true
            if (!pathExists(params.splitpipe_installation)) {
                errors << "'splitpipe_installation' = '${params.splitpipe_installation}' does not exist. Correct the path, or unset it to skip split-pipe.".toString()
            }
            if (!given(params.splitpipe_conda_env)) {
                errors << "'splitpipe_installation' is set, but 'splitpipe_conda_env' is not: split-pipe runs from that conda environment (see docs/INSTALLATION_EXTERNAL_PIPELINES.md).".toString()
            }
        }
        if (own == 'rhapsody' && given(params.rhapsody_installation)) {
            enabled = true
            if (!pathExists(params.rhapsody_installation)) {
                errors << "'rhapsody_installation' = '${params.rhapsody_installation}' does not exist. Correct the path, or unset it to skip the BD Rhapsody pipeline.".toString()
            }
        }

        // Each vendor pipeline runs under one kind of profile only: Cell Ranger has a container
        // but no conda environment, split-pipe a conda environment but no container, and the BD
        // Rhapsody pipeline starts its own Singularity containers, which it cannot do from inside one
        if (enabled && own == 'cellranger' && !containerEngine) {
            warnings << (
                "Cell Ranger has no conda environment: with no container engine (e.g. -profile conda) it runs\n" +
                "the 'cellranger' on the PATH of the compute nodes. Use a container profile (docker, singularity, apptainer, podman) to run it from its image."
            ).toString()
        }
        if (enabled && own == 'splitpipe' && containerEngine) {
            warnings << (
                "split-pipe has no container: under ${containerEngine}, 'splitpipe_conda_env' is not activated and split-pipe\n" +
                "runs from the PATH of the compute nodes. Use -profile conda (or mamba) to run it from that environment."
            ).toString()
        }
        if (enabled && own == 'rhapsody' && containerEngine) {
            warnings << (
                "The BD Rhapsody pipeline starts its own containers with 'cwltool --singularity', but under ${containerEngine}\n" +
                "cwltool runs inside a container that has no singularity. Use -profile conda, with singularity or apptainer on the PATH of the compute nodes."
            ).toString()
        }

        if (method == 'external_pipeline_only' && !enabled) {
            def stop = "'run_method' = 'external_pipeline_only' runs only the pre-processing and the vendor pipeline"
            def fix  = own == 'cellranger' ? "Set 'perform_cellranger = true' to run Cell Ranger"
                : own == 'splitpipe'      ? "Set 'splitpipe_installation' and 'splitpipe_conda_env' to run split-pipe"
                : own == 'rhapsody'       ? "Set 'rhapsody_installation' to run the BD Rhapsody pipeline"
                : null
            errors << (fix
                ? "${stop},\nbut no vendor pipeline is enabled for protocol '${protocol}'. ${fix}, or use 'run_method' = 'standard'."
                : "${stop},\nbut protocol '${protocol}' has no vendor pipeline, so only the pre-processing would run. Use 'run_method' = 'standard'."
            ).toString()
        }

        return [errors: errors, warnings: warnings]
    }
}
