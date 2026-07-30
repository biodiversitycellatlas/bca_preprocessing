//
// Workflow with functionality specific to 'main.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MTX_TO_H5AD             } from '../modules/local/tools/scanpy/main'
include { CELLBENDER              } from '../modules/local/tools/cellbender/main'
include { CELLSWEEP               } from '../modules/local/tools/cellsweep/main'
include { MTX_TO_10X              } from '../modules/local/tools/mtx_to_10x/main'
include { SCRUBLET                } from '../modules/local/tools/scrublet/main'
include { SCDBLFINDER             } from '../modules/local/tools/scdblfinder/main'
include { COMBINE_DOUBLET_RESULTS } from '../modules/local/tools/doublet_combine/main'
include { DOUBLET_FILTER as DOUBLET_FILTER_RAW         } from '../modules/local/tools/doublet_filter/main'
include { DOUBLET_FILTER as DOUBLET_FILTER_CELL_CALLED } from '../modules/local/tools/doublet_filter/main'
include { DOWNLOAD_DEMUXAFY_SIF   } from '../modules/local/tools/demuxafy_sif/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW TO RUN FILTERING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow filtering_workflow {
    take:
        ch_starsolo_genefull50_raw
        ch_starsolo_genefull50_filtered
        af_mtx

    main:
        // Initialize reporting channels
        def ch_cb_html                  = Channel.empty()
        def ch_cs_ambient_hist_plot     = Channel.empty()
        def ch_cs_umap_comparison_plot  = Channel.empty()
        def ch_cs_top_genes             = Channel.empty()
        def ch_scrublet_histogram       = Channel.empty()
        def ch_doublet_filter_summary   = Channel.empty()
        def ch_doublet_filter_plot      = Channel.empty()

        // Resolve each mapper's output directory into an (mtx, barcodes, features) triplet
        def ch_starsolo_raw = ch_starsolo_genefull50_raw.map { meta, dir ->
            [ meta + [mapping_method: 'starsolo', datatype: 'raw'],
              file("${dir}/matrix.mtx"), file("${dir}/barcodes.tsv"), file("${dir}/features.tsv") ]
        }

        def ch_starsolo_filtered = ch_starsolo_genefull50_filtered.map { meta, dir ->
            [ meta + [mapping_method: 'starsolo', datatype: 'filtered'],
              file("${dir}/matrix.mtx"), file("${dir}/barcodes.tsv"), file("${dir}/features.tsv") ]
        }

        def ch_alevin_fry = af_mtx.map { meta, dir ->
            [ meta + [mapping_method: 'alevin-fry', datatype: 'full'],
              file("${dir}/quants_mat.mtx"), file("${dir}/quants_mat_cols.txt"), file("${dir}/quants_mat_rows.txt") ]
        }

        def ch_all_matrices = ch_starsolo_raw.mix(ch_starsolo_filtered, ch_alevin_fry)

        // Doublets are called on the cell-called matrices only
        def ch_cell_called_matrices = ch_all_matrices.filter { meta, _mtx, _barcodes, _features ->
            meta.datatype in ['filtered', 'full']
        }

        MTX_TO_H5AD(ch_all_matrices)

        // Ambient RNA removal runs on the unfiltered matrices
        def ch_raw_h5ad = MTX_TO_H5AD.out.h5ad.filter { meta, _h5ad ->
            meta.datatype in ['raw', 'full']
        }

        // Ambient RNA removal
        if (params.ambient_rna_remover == "cellbender" || params.perform_cellbender) {
            CELLBENDER(ch_raw_h5ad)
            ch_cb_html = CELLBENDER.out.cb_html
            // TODO: create script to convert pdf to imgs and table

        } else if (params.ambient_rna_remover == "cellsweep") {

            if (params.perform_doublet_detection) {
                // Resolve the Demuxafy .sif image once per run: reuse params.demuxafy_sif if provided, otherwise download it
                DOWNLOAD_DEMUXAFY_SIF()
                def ch_demuxafy_sif = DOWNLOAD_DEMUXAFY_SIF.out.sif_path_file
                    .map { it.text.trim() }
                    .first()

                MTX_TO_10X(ch_cell_called_matrices)

                SCRUBLET(MTX_TO_10X.out.tenx_dir, ch_demuxafy_sif)
                SCDBLFINDER(MTX_TO_10X.out.tenx_dir, ch_demuxafy_sif)

                // Both tools ran on the same matrix, so they share a meta and join on it directly
                COMBINE_DOUBLET_RESULTS(
                    SCRUBLET.out.scrublet_results.join(SCDBLFINDER.out.scdblfinder_results),
                    ch_demuxafy_sif
                )

                // Join on sample + mapping method instead.
                def ch_doublets_by_sample = COMBINE_DOUBLET_RESULTS.out.combined_results
                    .map { meta, doublets -> [meta.subMap(['id', 'mapping_method']), doublets] }

                def ch_raw_h5ad_with_doublets = ch_raw_h5ad
                    .map { meta, h5ad -> [meta.subMap(['id', 'mapping_method']), meta, h5ad] }
                    .join(ch_doublets_by_sample)
                    .map { _key, meta, h5ad, doublets -> [meta, h5ad, doublets] }

                // Filter both raw and cell-called matrices for doublets
                if (params.perform_doublet_filtering) {
                    DOUBLET_FILTER_RAW(ch_raw_h5ad_with_doublets)
                    DOUBLET_FILTER_CELL_CALLED(
                        MTX_TO_H5AD.out.h5ad
                            .filter { meta, _h5ad -> meta.datatype == 'filtered' }
                            .join(COMBINE_DOUBLET_RESULTS.out.combined_results)
                    )

                    // The flagged cells are gone from the matrix, so there is nothing left for
                    // CellSweep to annotate or project; each DOUBLET_FILTER run reports what it
                    // removed in its own summary plot instead.
                    ch_cellsweep_input = DOUBLET_FILTER_RAW.out.h5ad.map { meta, h5ad -> [meta, h5ad, []] }

                    ch_doublet_filter_summary = DOUBLET_FILTER_RAW.out.summary
                        .mix(DOUBLET_FILTER_CELL_CALLED.out.summary)
                    ch_doublet_filter_plot    = DOUBLET_FILTER_RAW.out.filtering_summary_plot
                        .mix(DOUBLET_FILTER_CELL_CALLED.out.filtering_summary_plot)

                } else {
                    // The calls travel alongside the matrix for CellSweep to annotate and project
                    ch_cellsweep_input = ch_raw_h5ad_with_doublets
                }

                ch_scrublet_histogram = SCRUBLET.out.scrublet_histogram

            } else {
                // No doublet calls to project: nothing to stage into CellSweep's optional input
                ch_cellsweep_input = ch_raw_h5ad.map { meta, h5ad -> [meta, h5ad, []] }
            }

            CELLSWEEP(ch_cellsweep_input)
            ch_cs_ambient_hist_plot     = CELLSWEEP.out.cs_ambient_hist_plot
            ch_cs_umap_comparison_plot  = CELLSWEEP.out.cs_umap_comparison_plot
            ch_cs_top_genes             = CELLSWEEP.out.cs_top_genes
        }

    emit:
        cb_html                     = ch_cb_html
        cs_ambient_hist_plot        = ch_cs_ambient_hist_plot
        cs_umap_comparison_plot     = ch_cs_umap_comparison_plot
        cs_top_genes                = ch_cs_top_genes
        scrublet_histogram          = ch_scrublet_histogram
        doublet_filter_summary      = ch_doublet_filter_summary
        doublet_filter_plot         = ch_doublet_filter_plot
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
