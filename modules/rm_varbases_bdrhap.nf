// =============  REMOVE VARIABLE BASES FROM FASTQ ============= \\ 

process RM_VARBASES {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), path(fastq_files)

    output:
    tuple val(sample_id), path("noVB_${sample_id}_R{1,2}_001.fastq.gz"), emit: fastq_noVB_files

    script:
    """
    # Removing the Variable Bases and writing output to a new fastq
    zcat ${sample_id}_R1_001.fastq.gz | awk '{
            if (NR % 4 == 2) {  # Sequence
                l = length(\$0)
                sub(/^A|^GT|^TCA/, "", \$0)  # Remove variable bases (A, GT, TCA)
                d = l - length(\$0)  # Calculate how many bases were removed
                seq = \$0  # Store modified sequence
            } 
            else if (NR % 4 == 0) { 
            # Remove same number of characters from quality string
                \$0 = substr(\$0, d + 1)  
            }
        print
        }' > noVB_${sample_id}_R1_001.fastq
        
    gzip noVB_${sample_id}_R1_001.fastq
    # Renaming R2 by creating a new symbolic link
    ln -s "${sample_id}_R2_001.fastq.gz" "noVB_${sample_id}_R2_001.fastq.gz"
    """
}