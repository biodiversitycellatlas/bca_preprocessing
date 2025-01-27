// ===================  Filtering Workflow ================  \\ 
// 

include { CELLBENDER } from '../modules/cellbender'

workflow filtering_workflow {
    take:
        raw_matrix
    
    main:
        // Ambient RNA removal
        CELLBENDER(raw_matrix)
           
        filtered_h5 = CELLBENDER.out

    emit:
        filtered_h5
        
}
