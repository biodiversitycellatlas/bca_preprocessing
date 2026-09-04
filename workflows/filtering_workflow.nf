//
// Workflow with functionality specific to 'main.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MTX_TO_H5AD             } from '../modules/local/tools/scanpy/main'
include { CELLSWEEP               } from '../modules/local/tools/cellsweep/main'
include { MTX_TO_10X              } from '../modules/local/tools/mtx_to_10x/main'
include { SCRUBLET                } from '../modules/local/tools/scrublet/main'
include { SCDBLFINDER             } from '../modules/local/tools/scdblfinder/main'
include { COMBINE_DOUBLET_RESULTS } from '../modules/local/tools/doublet_combine/main'
include { DOUBLET_FILTER as DOUBLET_FILTER_RAW         } from '../modules/local/tools/doublet_filter/main'
include { DOUBLET_FILTER as DOUBLET_FILTER_CELL_CALLED } from '../modules/local/tools/doublet_filter/main'
include { DOWNLOAD_DEMUXAFY_SIF   } from '../modules/local/tools/demuxafy_sif/main'
include { COLLAPSE_ALEVIN_USA     } from '../modules/local/custom/manipulate/collapse_alevin_usa/main'
include { COLLAPSE_ALEVIN_USA as COLLAPSE_ALEVIN_UNSPLICED } from '../modules/local/custom/manipulate/collapse_alevin_usa/main'
include { VELOCITY_H5AD as VELOCITY_H5AD_STARSOLO } from '../modules/local/tools/velocity_h5ad/main'
include { VELOCITY_H5AD as VELOCITY_H5AD_ALEVIN   } from '../modules/local/tools/velocity_h5ad/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELPERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
 * Attach an optional per-sample annotation to a channel of matrices, appending it as the
 * last element and an empty list where it is absent.
 */
def attach_annotation(ch_left, ch_annotation, Closure key_of) {
    return ch_left
        .map { row -> [key_of.call(row[0])] + (row as List) }
        .join(ch_annotation, remainder: true)
        .filter { row -> row[1] != null }
        .map { row -> row[1..-2] + [row[-1] ?: []] }
}

/* A matrix directory written by mtx_io.write_triplet, resolved back into a triplet. */
def resolve_triplet(ch_dirs) {
    return ch_dirs.map { meta, dir ->
        [ meta, file("${dir}/matrix.mtx"), file("${dir}/barcodes.tsv"), file("${dir}/features.tsv") ]
    }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW TO RUN FILTERING
        The matrices travel as an (mtx, barcodes, features) triplet all the way through:
        doublet detection, optional doublet filtering and ambient-RNA denoising each read
        and write that triplet, and MTX_TO_H5AD runs last, packing one .h5ad per matrix
        that carries every annotation the stages before it produced.

        Each stage is optional and passes its input straight through when switched off, so
        with everything off this is the plain triplet -> h5ad conversion it has always been.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow filtering_workflow {
    take:
        ch_starsolo_genefull50_raw
        ch_starsolo_genefull50_filtered
        af_mtx
        af_filtered_mtx
        ch_starsolo_velocyto_filtered

    main:
        // Initialize reporting channels
        def ch_velocity_h5ad            = Channel.empty()
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

        // Alevin-fry's matrices are summed to one column per gene
        def ch_alevin_dirs = af_mtx
            .map { meta, dir -> [ meta + [mapping_method: 'alevin-fry', datatype: 'full'], dir ] }
            .mix(
                af_filtered_mtx.map { meta, dir ->
                    [ meta + [mapping_method: 'alevin-fry', datatype: 'filtered'], dir ]
                }
            )

        COLLAPSE_ALEVIN_USA(ch_alevin_dirs, params.alevin_usa_counts ?: 'SUA')

        // Fix alevin-fry's output into the same triplet format as STARsolo
        def ch_alevin_fry = COLLAPSE_ALEVIN_USA.out.matrix.map { meta, dir ->
            [ meta,
              file("${dir}/quants_mat.mtx"), file("${dir}/quants_mat_rows.txt"), file("${dir}/quants_mat_cols.txt") ]
        }

        def ch_all_matrices = ch_starsolo_raw.mix(ch_starsolo_filtered, ch_alevin_fry)

        /*
         * RNA velocity outputs. These are a terminal deliverable, deliberately kept out of
         * ch_all_matrices: everything after that funnel branches on meta.datatype, so mixing
         * them in would fan the h5ad, 10x-export, doublet and ambient-RNA steps out over the
         * extra datatypes -- and calling doublets or ambient RNA on unspliced counts is
         * meaningless.
         */
        if (params.perform_velocity) {

            // The unspliced block on its own, alevin-fry's counterpart of STARsolo's Velocyto
            // unspliced matrix. Retagged so it does not land in the same publishDir as the
            // gene-level matrix collapsed from the params.alevin_usa_counts blocks.
            def ch_alevin_unspliced_in = ch_alevin_dirs.map { meta, dir ->
                [ meta + [datatype: "${meta.datatype}_unspliced"], dir ]
            }

            COLLAPSE_ALEVIN_UNSPLICED(ch_alevin_unspliced_in, 'U')

            VELOCITY_H5AD_STARSOLO(ch_starsolo_velocyto_filtered.map { meta, dir ->
                [ meta + [mapping_method: 'starsolo', datatype: 'filtered'], dir ]
            })

            // The USA matrix carries all three blocks, so the velocity object is built from it
            // directly rather than from the collapsed output
            VELOCITY_H5AD_ALEVIN(ch_alevin_dirs)

            ch_velocity_h5ad = VELOCITY_H5AD_STARSOLO.out.h5ad.mix(VELOCITY_H5AD_ALEVIN.out.h5ad)
        }

        // Check the cell filtering method: unless this pipeline re-calls cells, alevin-fry's
        // full matrix is already the mapper's own cell call
        def alevin_full_is_cell_called = !(params.cellfilter_method in ["second_derivative", "manual_cutoff"])

        /*
         * Three views of the same matrices, each holding at most one entry per
         * [id, mapping_method] so that the annotation joins below stay unambiguous.
         */
        def ch_ambient_matrices = ch_all_matrices.filter { meta, _mtx, _barcodes, _features ->
            meta.datatype in ['raw', 'full']
        }

        def ch_filtered_matrices = ch_all_matrices.filter { meta, _mtx, _barcodes, _features ->
            meta.datatype == 'filtered'
        }

        def ch_cell_called_matrices = ch_all_matrices.filter { meta, _mtx, _barcodes, _features ->
            meta.datatype == 'filtered' || (meta.datatype == 'full' && alevin_full_is_cell_called)
        }

        // Keyed on sample + mapping method, so calls made on the cell-called matrix can be carried over to that sample's raw matrix
        def sample_key = { meta -> meta.subMap(['id', 'mapping_method']) }

        /*
         * Doublet detection, on the cell-called matrices of both mappers.
         *
         * Independent of params.ambient_rna_remover: the calls are an annotation in their own
         * right, and MTX_TO_H5AD carries them whether or not CellSweep ran.
         */
        def ch_calls_by_meta   = Channel.empty()
        def ch_calls_by_sample = Channel.empty()

        if (params.perform_doublet_detection) {

            // Resolve the Demuxafy .sif image once per run: reuse params.demuxafy_sif if provided, otherwise download it
            DOWNLOAD_DEMUXAFY_SIF()

            // Mount the .sif straight from the task work directory
            def ch_demuxafy_sif = DOWNLOAD_DEMUXAFY_SIF.out.sif_file
                .map { it.toRealPath().toString() }

            MTX_TO_10X(ch_cell_called_matrices)

            SCRUBLET(MTX_TO_10X.out.tenx_dir, ch_demuxafy_sif)
            SCDBLFINDER(MTX_TO_10X.out.tenx_dir, ch_demuxafy_sif)

            // Both tools ran on the same matrix, so they share a meta and join on it directly. 
            COMBINE_DOUBLET_RESULTS(
                SCRUBLET.out.scrublet_results.join(SCDBLFINDER.out.scdblfinder_results),
                ch_demuxafy_sif
            )

            ch_calls_by_meta      = COMBINE_DOUBLET_RESULTS.out.combined_results
            ch_calls_by_sample    = ch_calls_by_meta.map { meta, calls -> [sample_key.call(meta), calls] }
            ch_scrublet_histogram = SCRUBLET.out.scrublet_histogram

            // Warn once for a sample that lost a caller: it continues through every stage
            // below unannotated and unfiltered, rather than ending the run.
            ch_cell_called_matrices
                .map { meta, _mtx, _barcodes, _features -> [sample_key.call(meta), meta] }
                .join(ch_calls_by_sample, remainder: true)
                .filter { row -> row[1] != null && row[2] == null }
                .view { _key, meta, _calls ->
                    "[WARNING] No doublet calls for ${meta.id} (${meta.mapping_method}): a caller " +
                    "could not run on this matrix, which usually means too few cells were called. " +
                    "The sample continues without doublet annotation."
                }
        }

        /*
         * Optional doublet removal.
         */
        if (params.perform_doublet_filtering) {

            def ch_ambient_to_filter  = attach_annotation(ch_ambient_matrices, ch_calls_by_sample, sample_key)
            def ch_filtered_to_filter = attach_annotation(ch_filtered_matrices, ch_calls_by_meta, { meta -> meta })

            DOUBLET_FILTER_RAW(ch_ambient_to_filter.filter { row -> row[-1] })
            DOUBLET_FILTER_CELL_CALLED(ch_filtered_to_filter.filter { row -> row[-1] })

            ch_ambient_matrices = resolve_triplet(DOUBLET_FILTER_RAW.out.matrix)
                .mix(ch_ambient_to_filter.filter { row -> !row[-1] }.map { row -> row[0..-2] })

            ch_filtered_matrices = resolve_triplet(DOUBLET_FILTER_CELL_CALLED.out.matrix)
                .mix(ch_filtered_to_filter.filter { row -> !row[-1] }.map { row -> row[0..-2] })

            ch_doublet_filter_summary = DOUBLET_FILTER_RAW.out.summary
                .mix(DOUBLET_FILTER_CELL_CALLED.out.summary)
            ch_doublet_filter_plot    = DOUBLET_FILTER_RAW.out.filtering_summary_plot
                .mix(DOUBLET_FILTER_CELL_CALLED.out.filtering_summary_plot)
        }

        // The matrices as they now stand, with the consensus calls alongside them. With
        // detection switched off there is nothing to join against, so the empty slot is
        // filled in directly rather than by a join against an empty channel.
        def ch_ambient_with_calls = params.perform_doublet_detection
            ? attach_annotation(ch_ambient_matrices, ch_calls_by_sample, sample_key)
            : ch_ambient_matrices.map { row -> (row as List) + [[]] }

        def ch_filtered_with_calls = params.perform_doublet_detection
            ? attach_annotation(ch_filtered_matrices, ch_calls_by_meta, { meta -> meta })
            : ch_filtered_matrices.map { row -> (row as List) + [[]] }

        /*
         * Ambient RNA removal, on the unfiltered matrices.
         *
         * CellSweep gets exactly one of two things: the untouched matrix plus the calls, which
         * it annotates and projects onto its UMAP, or the matrix STAGE 2 already removed them
         * from and no calls at all. CellSweep's own guidance is to remove doublets first;
         * annotating them is the conservative default here, and removal stays opt-in.
         */
        def ch_cellsweep_by_sample = Channel.empty()

        if (params.ambient_rna_remover == "cellsweep") {

            CELLSWEEP(
                params.perform_doublet_filtering
                    ? ch_ambient_matrices.map { row -> (row as List) + [[]] }
                    : ch_ambient_with_calls
            )

            ch_cellsweep_by_sample      = CELLSWEEP.out.cs_full_h5ad.map { meta, h5ad -> [sample_key.call(meta), h5ad] }
            ch_cs_ambient_hist_plot     = CELLSWEEP.out.cs_ambient_hist_plot
            ch_cs_umap_comparison_plot  = CELLSWEEP.out.cs_umap_comparison_plot
            ch_cs_top_genes             = CELLSWEEP.out.cs_top_genes
        }

        /*
         * Published AnnData objects, one per matrix, assembled last so that
         * each carries whatever the stages above produced for it. Only the ambient matrices
         * have denoised counts to merge; the cell-called ones get the doublet calls alone.
         */
        def ch_h5ad_ambient = params.ambient_rna_remover == "cellsweep"
            ? attach_annotation(ch_ambient_with_calls, ch_cellsweep_by_sample, sample_key)
            : ch_ambient_with_calls.map { row -> (row as List) + [[]] }

        def ch_h5ad_filtered = ch_filtered_with_calls.map { row -> (row as List) + [[]] }

        MTX_TO_H5AD(ch_h5ad_ambient.mix(ch_h5ad_filtered))

    emit:
        h5ad                        = MTX_TO_H5AD.out.h5ad
        velocity_h5ad               = ch_velocity_h5ad
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
