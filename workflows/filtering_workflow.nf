include { CELLBENDER } from '../modules/cellbender'

workflow filtering_workflow {
    take:
        raw_matrix
    main:
        // Ambient RNA removal using CellBender
        CELLBENDER(raw_matrix)
    emit:
        CELLBENDER.out
        
}
