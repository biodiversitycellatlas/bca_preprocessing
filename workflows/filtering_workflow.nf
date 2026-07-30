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
include { DOUBLET_FILTER          } from '../modules/local/tools/doublet_filter/main'
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

        // Process STARsolo raw
        ch_starsolo_raw_ready = ch_starsolo_genefull50_raw.map { meta, dir ->
            def mtx      = file("${dir}/matrix.mtx")
            def barcodes = file("${dir}/barcodes.tsv")
            def features = file("${dir}/features.tsv")
            [meta, mtx, barcodes, features, 'starsolo', 'raw']
        }

        // Process STARsolo filtered
        ch_starsolo_filtered_ready = ch_starsolo_genefull50_filtered.map { meta, dir ->
            def mtx      = file("${dir}/matrix.mtx")
            def barcodes = file("${dir}/barcodes.tsv")
            def features = file("${dir}/features.tsv")
            [meta, mtx, barcodes, features, 'starsolo', 'filtered']
        }

        // Process alevin-fry
        ch_af_ready = af_mtx.map { meta, dir ->
            def mtx      = file("${dir}/quants_mat.mtx")
            def barcodes = file("${dir}/quants_mat_cols.txt")
            def features = file("${dir}/quants_mat_rows.txt")
            [meta, mtx, barcodes, features, 'alevin-fry', 'full']
        }

        // Combine all three channels
        ch_all_matrices = ch_starsolo_raw_ready.mix(ch_starsolo_filtered_ready, ch_af_ready)

        // Isolate the raw/full mtx triplets (pre-h5ad) for doublet detection
        ch_raw_mtx_triplet = ch_all_matrices
            .filter { meta, mtx, barcodes, features, mapping_method, datatype ->
                datatype == 'raw' || datatype == 'full'
            }

        // STARsolo's "raw" output includes every barcode with >=1 read (tens-to-hundreds of
        // thousands, mostly near-empty droplets), which breaks Scrublet/scDblFinder's automatic
        // doublet-rate estimate. STARsolo's own EmptyDrops-called ("filtered") cells are passed
        // through as an optional candidate-barcode list to restrict doublet calling to real cells.
        // alevin-fry's output is already knee-filtered at the permit-list stage, so it gets the
        // NO_FILE sentinel (no separate filtered channel flows into this workflow for it).
        ch_no_file = file("${projectDir}/assets/NO_FILE")
        ch_starsolo_filtered_barcodes = ch_starsolo_genefull50_filtered.map { meta, dir ->
            [[meta.id, 'starsolo', 'raw'], file("${dir}/barcodes.tsv")]
        }

        // Run the MTX to H5AD conversion for everything in the queue
        MTX_TO_H5AD(ch_all_matrices)

        // Isolate the raw/full data for ambient RNA removal
        ch_raw_h5ad = MTX_TO_H5AD.out.h5ad
            .filter { meta, h5ad, mapping_method, datatype ->
                // The 'filtered' datatype is ignored here
                datatype == 'raw' || datatype == 'full'
            }
            .map { meta, h5ad, mapping_method, datatype ->
                [meta, h5ad, mapping_method, datatype]
            }

        // Ambient RNA removal
        if (params.ambient_rna_remover == "cellbender" || params.perform_cellbender) {
            CELLBENDER(ch_raw_h5ad)
            ch_cb_html = CELLBENDER.out.cb_html
            // TODO: create script to convert pdf to imgs and table

        } else if (params.ambient_rna_remover == "cellsweep") {

            if (params.perform_doublet_filtering) {
                // Resolve the Demuxafy .sif image once per run: reuse params.demuxafy_sif if provided, otherwise download it
                DOWNLOAD_DEMUXAFY_SIF()
                ch_demuxafy_sif = DOWNLOAD_DEMUXAFY_SIF.out.sif_path_file
                    .map { it.text.trim() }
                    .first()

                MTX_TO_10X(ch_raw_mtx_triplet)

                // Attach the optional filtered-barcode list by the same composite key.
                ch_tenx_with_filter = MTX_TO_10X.out.tenx_dir
                    .map { meta, tenx_dir, mapping_method, datatype -> [[meta.id, mapping_method, datatype], meta, tenx_dir] }
                    .join(ch_starsolo_filtered_barcodes, remainder: true)
                    .map { key, meta, tenx_dir, filtered_barcodes -> [meta, tenx_dir, key[1], key[2], filtered_barcodes ?: ch_no_file] }

                SCRUBLET(ch_tenx_with_filter, ch_demuxafy_sif)
                SCDBLFINDER(ch_tenx_with_filter, ch_demuxafy_sif)

                // Composite key (meta.id, mapping_method, datatype): ch_raw_h5ad can carry
                // multiple entries per meta.id (e.g. starsolo 'raw' and alevin-fry 'full'),
                // so joining on meta alone would cross-join those branches.
                ch_scrublet_keyed = SCRUBLET.out.scrublet_results
                    .map { meta, tsv, mapping_method, datatype -> [[meta.id, mapping_method, datatype], meta, tsv] }
                ch_scdblfinder_keyed = SCDBLFINDER.out.scdblfinder_results
                    .map { meta, tsv, mapping_method, datatype -> [[meta.id, mapping_method, datatype], tsv] }

                ch_doublet_pairs = ch_scrublet_keyed
                    .join(ch_scdblfinder_keyed)
                    .map { key, meta, scrublet_tsv, scdblfinder_tsv -> [meta, scrublet_tsv, scdblfinder_tsv, key[1], key[2]] }

                COMBINE_DOUBLET_RESULTS(ch_doublet_pairs, ch_demuxafy_sif)

                ch_raw_h5ad_keyed = ch_raw_h5ad
                    .map { meta, h5ad, mapping_method, datatype -> [[meta.id, mapping_method, datatype], h5ad] }
                ch_combined_keyed = COMBINE_DOUBLET_RESULTS.out.combined_results
                    .map { meta, tsv, mapping_method, datatype -> [[meta.id, mapping_method, datatype], meta, tsv, mapping_method, datatype] }

                ch_doublet_filter_input = ch_combined_keyed
                    .join(ch_raw_h5ad_keyed)
                    .map { key, meta, combined_tsv, mapping_method, datatype, h5ad -> [meta, h5ad, combined_tsv, mapping_method, datatype] }

                DOUBLET_FILTER(ch_doublet_filter_input)

                ch_cellsweep_input          = DOUBLET_FILTER.out.h5ad
                ch_doublet_filter_summary   = DOUBLET_FILTER.out.summary
                ch_doublet_filter_plot      = DOUBLET_FILTER.out.filtering_summary_plot
                ch_scrublet_histogram       = SCRUBLET.out.scrublet_histogram

            } else {
                ch_cellsweep_input = ch_raw_h5ad
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
