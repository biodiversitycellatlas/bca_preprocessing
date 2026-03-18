//
// Workflow with functionality specific to 'main.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MTX_TO_H5AD   } from '../modules/local/tools/scanpy/main'
include { CELLBENDER    } from '../modules/local/tools/cellbender/main'
include { CELLSWEEP     } from '../modules/local/tools/cellsweep/main'


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
        } else if (params.ambient_rna_remover == "cellsweep") {
            CELLSWEEP(ch_raw_h5ad)
        }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
