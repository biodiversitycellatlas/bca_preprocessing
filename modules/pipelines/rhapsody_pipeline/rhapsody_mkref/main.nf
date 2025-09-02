process BDRHAP_PIPELINE_MKREF {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    output:
    path("*")

    script:
    """
    echo "\n\n===============  BD Rhapsody pipeline - mkref  ==============="

    ref_dir="./BDrhap_reference"
    
    mkdir -p \${ref_dir}
    cd \${ref_dir}

    # Retrieve the first accession number
    first_fastq=\$(ls "${params.outdir}/fastq/" | head -n1)  
    echo "\${first_fastq}"

    # Calculate SJDB overhang using the first read from the first fastq file
    sjdb_overhang=\$(zcat ${params.outdir}/fastq/\${first_fastq} 2>/dev/null | awk 'NR==2 {print length(\$0)-1; exit}' || echo "") 

    STAR --runMode genomeGenerate \\
        --genomeFastaFiles ${params.ref_fasta} \\
        --sjdbGTFfile ${params.ref_gtf} \\
        --sjdbOverhang 150 \\
        --genomeSAindexNbases 12 \\
        --genomeSAsparseD 1
    
    # Navigate a directory down to copy the GTF inside the folder
    mv GenomeDir/ star_index/

    # Copy the reference file to the BD Rhapsody reference directory
    cp ${params.ref_gtf} .

    # Navigate a directory down to compress the reference files
    cd ../

    # Compressing the directory into tar gz format
    tar -czvf BD_Rhapsody_Reference_Files.tar.gz BDrhap_reference/
    """
}