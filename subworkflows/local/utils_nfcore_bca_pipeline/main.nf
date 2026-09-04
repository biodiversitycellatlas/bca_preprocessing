//
// Subworkflow with functionality specific to the BCA pre-processing pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba', 'micromamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Fail fast on doublet options that would silently do nothing
    //
    // Doublet detection annotates the consensus calls, doublet filtering additionally removes
    // them from the matrix. Both are independent of 'ambient_rna_remover': the calls reach the
    // published AnnData through MTX_TO_H5AD whether or not CellSweep ran.
    //
    if (params.perform_doublet_filtering && !params.perform_doublet_detection) {
        error(
            "'perform_doublet_filtering' is set, but 'perform_doublet_detection' is false.\n" +
            "There would be no doublet calls to filter on: either enable 'perform_doublet_detection',\n" +
            "or unset 'perform_doublet_filtering'."
        )
    }

    //
    // Fail fast on an unknown run method or cell-calling method
    //
    def valid_run_methods = ["standard", "geneext_only", "external_pipeline_only", "post_mapping"]
    if (!(params.run_method in valid_run_methods)) {
        error("Unknown 'run_method' = '${params.run_method}'. Use one of: ${valid_run_methods.join(', ')}.")
    }

    def valid_cellfilter_methods = ["star_solocellfilter", "second_derivative", "manual_cutoff"]
    if (!(params.cellfilter_method in valid_cellfilter_methods)) {
        error("Unknown 'cellfilter_method' = '${params.cellfilter_method}'. Use one of: ${valid_cellfilter_methods.join(', ')}.")
    }

    //
    // Fail fast on a 'geneext_only' run that cannot produce a BAM
    //
    // GeneExt reads the alignments, so without a BAM there is nothing to extend from. Left
    // unchecked the run succeeds and writes no annotation at all: the BAM channel stays
    // empty, collect() on it emits nothing, and the merge and GeneExt never run.
    //
    if (params.run_method == "geneext_only" && !params.star_generateBAM && !params.geneext_bam_only) {
        error(
            "'run_method' = 'geneext_only' with both 'star_generateBAM' and 'geneext_bam_only'\n" +
            "unset produces no BAM for GeneExt to read, so no extended annotation would be\n" +
            "written. Leave 'geneext_bam_only' set, or enable 'star_generateBAM'."
        )
    }

    //
    // Create channel from input file provided through params.input
    //
    def samples_file = file(params.input)

    //
    // Fail fast on a 'post_mapping' run with nothing to pick up
    //
    // The mode reads a previous run's published mapping results back in rather than
    // mapping again, so those results have to be there before anything is scheduled.
    //
    if (params.run_method == "post_mapping") {
        def prev_dir = file(params.previous_outdir ?: params.outdir)
        if (!prev_dir.exists()) {
            error(
                "'run_method' = 'post_mapping' reads a previous run's mapping results from\n" +
                "${prev_dir}, which does not exist. Set 'previous_outdir' to that run's results directory."
            )
        }
        if (!file("${prev_dir}/mapping_STARsolo").exists() && !file("${prev_dir}/mapping_alevin").exists()) {
            error(
                "'run_method' = 'post_mapping' found neither 'mapping_STARsolo/' nor 'mapping_alevin/'\n" +
                "under ${prev_dir}, so there are no mapping results to pick up.\n" +
                "Set 'previous_outdir' to the results directory of the run to continue from."
            )
        }
    }

    //
    // Fail fast on manual cutoffs that are not actually given
    //
    // The cutoff is per sample and has no sensible default, so a missing one would
    // silently fall back to the second-derivative call it was meant to replace.
    //
    def samplesheet_rows = samples_file.splitCsv(header: true, sep: ',')

    // An empty column is not a cutoff, so it is not worth warning about
    def has_cutoffs = samplesheet_rows.any { row -> row.manual_cutoff?.toString()?.trim() }

    if (params.cellfilter_method == "manual_cutoff") {
        def without_cutoff = samplesheet_rows
            .findAll { row -> !(row.manual_cutoff?.toString()?.trim()) }
            .collect { row -> row.sample }
            .unique()

        if (without_cutoff) {
            error(
                "'cellfilter_method' = 'manual_cutoff' needs a 'manual_cutoff' column in the samplesheet\n" +
                "giving the UMI threshold for every sample, but it is missing for: ${without_cutoff.join(', ')}."
            )
        }
    } else if (has_cutoffs) {
        log.warn(
            "The samplesheet gives a 'manual_cutoff' for one or more samples, but 'cellfilter_method' is " +
            "'${params.cellfilter_method}'. Those cutoffs will be ignored; set " +
            "'cellfilter_method' to 'manual_cutoff' to apply them."
        )
    }

    Channel
    .fromPath(params.input)
    .splitCsv(header: true, sep: ',')
    .map { row ->
        // Concat metadata into a Map object
        def meta = [
            id             : row.sample,
            expected_cells : row.expected_cells ? row.expected_cells.toInteger() : null,
            manual_cutoff  : row.manual_cutoff ? row.manual_cutoff.toInteger() : null,
            p5             : row.p5 ? row.p5 : '',
            p7             : row.p7 ? row.p7 : '',
            rt             : row.rt ? row.rt : '',
        ]

        // Assign each FASTQ to its own variable (or null if missing)
        def fastq_cDNA   = row.fastq_cDNA   ? [file(row.fastq_cDNA)]   : []
        def fastq_BC_UMI = row.fastq_BC_UMI ? [file(row.fastq_BC_UMI)] : []
        def fastq_indices = row.fastq_indices ? [file(row.fastq_indices)] : []

        // Returng tuple
        [ meta.id, meta, fastq_cDNA, fastq_BC_UMI, fastq_indices, samples_file ]
    }
    // Group by sample ID
    .groupTuple(by: 0)

    .map { sample_id, metas, fastq_cDNAs, fastq_BC_UMIs, fastq_indices_list, samples_files ->

        def meta = metas[0]

        // Filter out empty lists, flatten valid files only
        def merged_fastq_cDNA = fastq_cDNAs.findResults { it }.flatten().findAll { it.exists() }
        def merged_fastq_BC_UMI = fastq_BC_UMIs.findResults { it }.flatten().findAll { it.exists() }
        def merged_fastq_indices = fastq_indices_list.findResults { it }.flatten()

        samples_file = samples_files[0]

        // Return tuple for merging FASTQs
        tuple(meta, merged_fastq_cDNA, merged_fastq_BC_UMI, merged_fastq_indices, samples_file)
    }
    .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}


//
// Generate methods description for MultiQC
//
def toolCitationText() {

    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics , 34(17), i884-i890. doi: /10.1093/bioinformatics/bty560</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047-3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
