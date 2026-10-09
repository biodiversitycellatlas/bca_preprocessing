//
// Workflow with functionality specific to 'main.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MC_GENE_TABLE             } from '../modules/local/tools/metacells/gene_table/main'
include { METACELL2                 } from '../modules/local/tools/metacells/run/main'
include { GENERATE_FILTERING_REPORT } from '../modules/local/custom/filtering_report/main'
include { APPLY_METACELL_FILTER     } from '../modules/local/tools/metacells/apply/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW TO GROUP CELLS INTO METACELLS AND FILTER THEM
        Runs on the h5ad files the filtering workflow publishes. Each sample's cell-called
        object is grouped into Metacell2 metacells, and every group is summarised for
        filtering_report.html, where the user chooses the metacells and genes to keep.

        That choice cannot be made inside a run, so applying it is a second pass:
        re-running with -resume and --metacell_selection reuses every task above and
        runs only APPLY_METACELL_FILTER, which writes the filtered UMI matrices.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow metacell_workflow {
    take:
        ch_h5ad
        ref_gtf
        mt_rrna_metrics
        antisense_metrics

    main:
        def ch_final_h5ad = channel.empty()

        // Same rule as filtering_workflow: unless this pipeline re-calls cells, alevin-fry's
        // full matrix is already the mapper's own cell call
        def alevin_full_is_cell_called = !(params.cellfilter_method in ["second_derivative", "manual_cutoff"])
        def sample_key = { meta -> meta.subMap(['id', 'mapping_method']) }

        def ch_candidates = ch_h5ad.filter { meta, _h5ad ->
            !(params.mc2_skip_subsampled && meta.id.endsWith('_subsampled_starsolo'))
        }

        // One cell-called object per sample and mapper; where a mapper offers both, the
        // re-called 'filtered' matrix is the one the doublet calls were made on
        def ch_cell_called = ch_candidates
            .filter { meta, _h5ad ->
                meta.datatype == 'filtered' || (meta.datatype == 'full' && alevin_full_is_cell_called)
            }
            .map { meta, h5ad -> [sample_key.call(meta), meta, h5ad] }
            .groupTuple(by: 0)
            .map { _key, metas, h5ads ->
                def i = metas.findIndexOf { m -> m.datatype == 'filtered' }
                i = i >= 0 ? i : 0
                [sample_key.call(metas[i]), metas[i], h5ads[i]]
            }

        /*
         * CellSweep only runs on the unfiltered matrix, so its annotations are lifted from
         * that sample's raw (STARsolo) or full (alevin-fry) object. Where the cell-called
         * object is that same full object, it already carries them.
         */
        def ch_mc_input
        if (params.ambient_rna_remover == "cellsweep") {
            def ch_ambient = ch_candidates
                .filter { meta, _h5ad -> meta.datatype in ['raw', 'full'] }
                .map { meta, h5ad -> [sample_key.call(meta), meta, h5ad] }

            ch_mc_input = ch_cell_called
                .join(ch_ambient, remainder: true)
                .filter { row -> row[1] instanceof Map }
                .map { row ->
                    def meta    = row[1]
                    def ameta   = row.size() > 3 ? row[3] : null
                    def ambient = row.size() > 4 ? row[4] : null
                    def usable  = ambient && ameta instanceof Map && ameta.datatype != meta.datatype
                    [meta, row[2], usable ? ambient : []]
                }
        } else {
            ch_mc_input = ch_cell_called.map { _key, meta, h5ad -> [meta, h5ad, []] }
        }

        // Built once from the reference annotation and shared by every sample
        MC_GENE_TABLE(
            ref_gtf.first(),
            params.gene_annotation ? file(params.gene_annotation, checkIfExists: true) : []
        )

        METACELL2(ch_mc_input, MC_GENE_TABLE.out.gene_table.first())

        // Sorted, so the collected inputs hash the same on every -resume
        def thresholds = [
            max_mito_pct    : params.mcfilter_max_mito_pct,
            max_doublet_pct : params.mcfilter_max_doublet_pct,
            max_alpha_hat   : params.mcfilter_max_alpha_hat,
            min_gene_umis   : params.mcfilter_min_gene_umis,
            pfam_presets    : params.mcfilter_pfam_presets,
        ]

        GENERATE_FILTERING_REPORT(
            METACELL2.out.summary.map { _meta, summary -> summary }.collect(sort: true).ifEmpty([]),
            mt_rrna_metrics.collect(sort: true).ifEmpty([]),
            antisense_metrics.collect(sort: true).ifEmpty([]),
            file("${projectDir}/bin/filtering_report.html", checkIfExists: true),
            file("${projectDir}/assets/bca_logo.svg", checkIfExists: true),
            METACELL2.out.summary
                .filter { meta, _summary -> meta.expected_cells != null }
                .map { meta, _summary -> "${meta.id}=${meta.expected_cells}".toString() }
                .collect(sort: true).ifEmpty([]),
            thresholds
        )

        def ch_versions = MC_GENE_TABLE.out.versions.mix(
            METACELL2.out.versions,
            GENERATE_FILTERING_REPORT.out.versions
        )

        /*
         * Second pass: apply the selection exported from the report. Only the samples it
         * covers are filtered; a sample left out of it is reported, not failed.
         */
        if (params.metacell_selection) {
            def selection_file = file(params.metacell_selection, checkIfExists: true)
            def selected_ids   = (new groovy.json.JsonSlurper().parse(selection_file).samples ?: [:]).keySet()

            def ch_to_apply = METACELL2.out.cells_h5ad.branch { meta, _h5ad ->
                selected: meta.id in selected_ids
                other:    true
            }

            ch_to_apply.other.view { meta, _h5ad ->
                "[WARNING] ${meta.id} is not in ${selection_file.name}, so no filtered matrix is written for it. " +
                "Export a selection covering it from filtering_report.html."
            }

            APPLY_METACELL_FILTER(ch_to_apply.selected, selection_file)
            ch_final_h5ad = APPLY_METACELL_FILTER.out.h5ad
            ch_versions   = ch_versions.mix(APPLY_METACELL_FILTER.out.versions)
        }

    emit:
        report          = GENERATE_FILTERING_REPORT.out.html
        summaries       = METACELL2.out.summary
        cells_h5ad      = METACELL2.out.cells_h5ad
        metacells_h5ad  = METACELL2.out.metacells_h5ad
        final_h5ad      = ch_final_h5ad
        versions        = ch_versions
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
