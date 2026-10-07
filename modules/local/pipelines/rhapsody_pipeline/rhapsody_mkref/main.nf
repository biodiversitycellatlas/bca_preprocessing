process BDRHAP_PIPELINE_MKREF {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    // The STAR index is built by STARSOLO_INDEX; it is staged straight into the layout the
    // BD Rhapsody pipeline expects inside its reference archive
    input:
    path(star_index, stageAs: 'BD_Rhapsody_Reference_Files/star_index/*')
    path ref_gtf

    output:
    path("BD_Rhapsody_Reference_Files.tar.gz"), emit: reference
    path "versions.yml",                        emit: versions

    script:
    """
    echo -e "\\n\\n===============  BD Rhapsody pipeline - mkref  ==============="

    # Copy the GTF into the reference directory
    cp -L ${ref_gtf} BD_Rhapsody_Reference_Files/

    # Create the tar.gz archive with the expected name and structure. -h stores the staged
    # index files rather than the links to them; STARSOLO_INDEX's versions.yml is left out.
    tar -chzvf BD_Rhapsody_Reference_Files.tar.gz --exclude='versions.yml' BD_Rhapsody_Reference_Files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: \$(tar --version | head -n1 | sed 's/^.* //')
    END_VERSIONS
    """
}
