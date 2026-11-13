process GENE_EXT {
    publishDir path: { "${params.outdir}/gene_ext" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    
    // Peaks of 1.4 GB (C. elegans) and 2.8 GB (fly) do not follow the BAM size, so memory is
    // fixed; time follows the BAM (params.dynamic_time), never below 2 h.
    cpus   4
    memory { 6.GB * task.attempt }
    time   { BcaResources.scaledTime(
        params.dynamic_time?.GENE_EXT, [bam_file], task.attempt, 2) }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_bedtools_samtools_pysam_pruned:f0dc809bbf3297e3' :
        'community.wave.seqera.io/library/python_bedtools_samtools_pysam_pruned:924a6b1d51a277f5' }"
    
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
        -j ${task.cpus} \\
        -force ${subsamplebam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
