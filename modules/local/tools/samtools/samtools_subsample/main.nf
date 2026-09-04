process SAMTOOLS_SUBSAMPLE {
    tag "${meta.id}"
    label 'process_low2'

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/samtools:1.22.1--9a10f06c24cdf05f"

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("${meta.id}_subsampled.bam"), emit: subsampled_bam
    path "versions.yml",                               emit: versions

    script:
    // Matches the seed hard-coded in SUBSAMPLE_FASTQS, so both subsamplers are reproducible the same way
    def seed = 100
    """
    echo "\n\n==================  SUBSAMPLE BAM  =================="
    echo "Sample ID: ${meta.id}"
    echo "BAM file: ${bam_file}"
    echo "Target reads: ${params.geneext_subsample_nreads}"

    TARGET=${params.geneext_subsample_nreads}

    # -F 0x900 drops secondary and supplementary records, so this counts reads rather
    # than alignments and the target means what it says
    TOTAL=\$(samtools view -c -@ ${task.cpus} -F 0x900 ${bam_file})

    # "none" means there is nothing to remove: either the sample is already at or under the
    # target, or the fraction rounds up to 1 and rewriting the BAM would drop no reads
    FRAC=\$(awk -v t="\$TOTAL" -v n="\$TARGET" 'BEGIN {
        if (t <= n) { print "none"; exit }
        f = n / t
        if (f >= 0.999999) { print "none"; exit }
        printf "%.6f", f
    }')

    if [[ "\$FRAC" == "none" ]]; then
        echo "\$TOTAL reads is at or under the target of \$TARGET: passing the BAM through unchanged"
        ln -s ${bam_file} ${meta.id}_subsampled.bam
    else
        echo "Subsampling \$TOTAL reads down to ~\$TARGET (fraction \$FRAC, seed ${seed})"
        # --subsample decides per read name, so every record of a read is kept or dropped together
        samtools view -@ ${task.cpus} -b \\
            --subsample "\$FRAC" \\
            --subsample-seed ${seed} \\
            -o ${meta.id}_subsampled.bam \\
            ${bam_file}

        echo "Reads after subsampling: \$(samtools view -c -@ ${task.cpus} -F 0x900 ${meta.id}_subsampled.bam)"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed '1!d; s/samtools //')
    END_VERSIONS
    """
}
