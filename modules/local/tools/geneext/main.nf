process GENE_EXT {
    publishDir path: { "${params.outdir}/gene_ext" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    label 'process_long'

    // conda "${projectDir}/submodules/GeneExt/environment.yaml"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_bedtools_samtools_pysam_pruned:1deba2d2d8b44e3b' :
        'community.wave.seqera.io/library/python_bedtools_samtools_pysam_pruned:baa04298159443b3' }"
    // The script lives in a submodule. A container only sees the task directory and bin/,
    // unless the pipeline happens to share a parent directory with the work directory
    containerOptions { "${workflow.containerEngine in ['docker', 'podman'] ? '-v' : '-B'} ${projectDir}/submodules/GeneExt:${projectDir}/submodules/GeneExt" }

    input:
    path(bam_file)
    path(bam_index)
    path(ref_gtf, stageAs: 'ref_input/*')   // own folder, so it cannot clash with geneext.gtf

    output:
    path("geneext.gtf")                 , emit: gtf
    path("geneext.gtf.Report.html")     , emit: report , optional: true
    path("geneext.gtf.GeneExt.log")     , emit: log    , optional: true
    path "versions.yml"                 , emit: versions

    script:
    def subsamplebam = params.geneext_subsamplebam ? "--subsamplebam ${params.geneext_subsamplebam}" : ""
    """
    echo -e "\\n\\n==================  GENE EXTENSION =================="
    echo "BAM file: ${bam_file}"
    echo "BAM index: ${bam_index}"
    echo "Original GTF: ${ref_gtf}"

    # Remove temporary directory if it exists
    if [ -d "tmp" ]; then rm -r tmp; fi

    # Run GeneExt
    python ${projectDir}/submodules/GeneExt/geneext.py \\
        -g ${ref_gtf} \\
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
