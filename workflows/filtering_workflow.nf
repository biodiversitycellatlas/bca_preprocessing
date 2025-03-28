include { CELLBENDER } from '../modules/cellbender'

workflow filtering_workflow {
    take:
        raw_matrix
    main:
        // Ambient RNA removal using CellBender
        if (params.perform_cellbender) {
            CELLBENDER(raw_matrix)
        } else {
            log.info "Skipping Cellbender steps as 'perform_cellbender' is false."
        }
}
