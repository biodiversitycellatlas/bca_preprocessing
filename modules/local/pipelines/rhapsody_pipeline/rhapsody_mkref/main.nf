process BDRHAP_PIPELINE_MKREF {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    output:
    path("BD_Rhapsody_Reference_Files.tar.gz")

    script:
    """
    echo "\n\n===============  BD Rhapsody pipeline - mkref  ==============="

    # Directory matches what the BD Rhapsody pipeline expects
    ref_dir="BD_Rhapsody_Reference_Files"
    mkdir -p "\${ref_dir}"

    # Retrieve the first accession number
    first_fastq=\$(ls "${params.outdir}/fastq/" | head -n1)
    echo "First FASTQ: \${first_fastq}"

    # Calculate SJDB overhang using the first read from the first fastq file
    sjdb_overhang=\$(zcat ${params.outdir}/fastq/\${first_fastq} 2>/dev/null | awk 'NR==2 {print length(\$0)-1; exit}' || echo "")
    echo "SJDB overhang: \${sjdb_overhang}"

    # Build STAR genome index directly into <ref_dir>/star_index
    STAR --runMode genomeGenerate \\
        --genomeDir "\${ref_dir}/star_index" \\
        --genomeFastaFiles ${params.ref_fasta} \\
        --sjdbGTFfile ${params.ref_gtf} \\
        --sjdbOverhang 150 \\
        --genomeSAindexNbases 12 \\
        --genomeSAsparseD 1

    # Copy the GTF into the reference directory
    cp ${params.ref_gtf} "\${ref_dir}/"

    # Create the tar.gz archive with the expected name and structure
    tar -czvf BD_Rhapsody_Reference_Files.tar.gz "\${ref_dir}"
    """
}
