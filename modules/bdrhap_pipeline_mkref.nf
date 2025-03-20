process BDRHAP_PIPELINE_MKREF {
    output:
    file("BD_Rhapsody_Reference_Files.tar.gz")

    script:
    """
    echo "\n\n===============  BD Rhapsody pipeline - mkref  ==============="

    ref_dir="${params.resDir}/genome/BDrhap_reference/star_index"
    
    mkdir -p \$ref_dir
    cd \$ref_dir

    # Retrieve the first accession number
    first_fastq=\$(ls "${params.resDir}/fastq/" | head -n1)  
    echo "\${first_fastq}"

    # Calculate SJDB overhang using the first read from the first fastq file
    sjdb_overhang=\$(zcat ${params.resDir}/fastq/\${first_fastq} 2>/dev/null | awk 'NR==2 {print length(\$0)-1; exit}' || echo "") 

    STAR --runMode genomeGenerate \\
        --genomeFastaFiles ${params.ref_fasta} \\
        --sjdbGTFfile ${params.ref_star_gtf} \\
        --sjdbOverhang 150 \\
        --genomeSAindexNbases 12 \\
        --genomeSAsparseD 1

    # Copy the reference file to the BD Rhapsody reference directory
    cp ${params.ref_star_gtf} \$ref_dir

    # Navigate two directories down to compress the reference files
    cd ../../

    # Compressing the directory into tar gz format
    tar -czvf BD_Rhapsody_Reference_Files.tar.gz BDrhap_reference/

    # Final PATH
    # /users/asebe/bvanwaardenburg/git/data/241106_BD_Rhapsody_Nvec/Nvec_BDpipeline/BD_Rhapsody_Reference_Files.tar.gz
    """
}