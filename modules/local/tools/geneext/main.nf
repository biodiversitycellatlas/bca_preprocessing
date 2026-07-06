process GENE_EXT {
    publishDir "${params.outdir}/gene_ext", mode: 'copy'
    label 'process_high'

    // conda "${projectDir}/submodules/GeneExt/environment.yaml"
    conda "${moduleDir}/environment.yml"

    input:
    path(bam_file)
    path(bam_index)

    output:
    path("geneext.gtf"), emit: gtf
    path "versions.yml", emit: versions

    script:
    def subsamplebam = params.geneext_subsamplebam ? "--subsamplebam ${params.geneext_subsamplebam}" : ""
    """
    echo "\n\n==================  GENE EXTENSION =================="
    echo "BAM file: ${bam_file}"
    echo "BAM index: ${bam_index}"
    echo "Original GTF: ${params.ref_gtf}"

    # Remove temporary directory if it exists
    if [ -d "tmp" ]; then rm -r tmp; fi

    # Run GeneExt
    python ${projectDir}/submodules/GeneExt/geneext.py \\
        -g ${params.ref_gtf} \\
        -b ${bam_file} \\
        -o geneext.gtf \\
        -t tmp_geneext \\
        -j 4 \\
        -force ${subsamplebam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
