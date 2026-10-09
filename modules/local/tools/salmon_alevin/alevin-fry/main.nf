process ALEVIN_FRY {
    publishDir path: { "${params.outdir}/mapping_alevin/${meta.id}" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${meta.id}"
    label 'process_high'

    // Exit 42 from bin/mapping_rate_guard.sh (mapping rate below params.min_mapping_rate)
    // aborts the pipeline or drops the sample; other exit codes follow conf/base.config.
    errorStrategy { MappingRateCheck.errorStrategy(task, params) }

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/alevin-fry_piscem_salmon_simpleaf_pruned:c71cfb476b003414"

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file)
    path(bc_whitelist, stageAs: 'whitelist_?/*')  // in barcode-segment order; empty without a whitelist
    path(splici_index_reference)
    path(salmon_index)

    output:
    tuple val(meta), path("${meta.id}_*"), emit: mapping_files
    tuple val(meta), path("${meta.id}_run/aux_info/meta_info.json"),    emit: af_meta_info
    tuple val(meta), path("${meta.id}_counts/quant.json"),              emit: af_quant_json
    tuple val(meta), path("${meta.id}_counts/alevin"),                  emit: af_mtx
    tuple val(meta), path("${meta.id}_counts/cell_meta.tsv"),           emit: af_cell_meta, optional: true
    path "versions.yml",                                                emit: versions

    script:
    // Retrieve settings from custom parameters if set, otherwise from conf/seqtech_parameters.config
    def bc_geom = params.alevin_bc_geometry ?: params.seqtech_parameters[params.protocol].alevin_bc_geometry
    def umi_geom = params.alevin_umi_geometry ?: params.seqtech_parameters[params.protocol].alevin_umi_geometry
    def read_geom = params.alevin_read_geometry ?: params.seqtech_parameters[params.protocol].alevin_read_geometry

    if (!bc_geom || !umi_geom || !read_geom) {
        error "No alevin geometry defined for protocol '${params.protocol}'. Set 'alevin_bc_geometry', 'alevin_umi_geometry' and 'alevin_read_geometry' in the configuration file."
    }

    // Mapping-rate check: salmon runs under bin/mapping_rate_guard.sh, which reads
    // percent_mapped from meta_info.json before the alevin-fry steps, and cancels salmon
    // early if STARsolo already found this sample below the threshold.
    def mapping_guard = MappingRateCheck.enabled(params)
        ? "mapping_rate_guard.sh run --mapper salmon" +
          " --flag-dir ${MappingRateCheck.flagDir(workflow)} --sample ${meta.base_id ?: meta.id} --label ${meta.id}" +
          " --min ${params.min_mapping_rate} --poll ${params.mapping_rate_poll_secs ?: 60}" +
          " --report ./${meta.id}_run/aux_info/meta_info.json --"
        : ''

    // Count reads on the strand STARsolo counts (--soloStrand), so both mappers see the same molecules
    def expected_ori = [Forward: 'fw', Reverse: 'rc', Unstranded: 'both'][params.star_soloStrand]
    if (!expected_ori) {
        error "star_soloStrand '${params.star_soloStrand}' has no alevin-fry orientation; use Forward, Reverse or Unstranded."
    }

    // Permit list. When the pipeline re-calls cells itself, alevin-fry keeps every barcode on
    // the whitelist (--unfiltered-pl), so its 'full' matrix is a raw matrix like STARsolo's and
    // the cell call is made on the whole curve. With star_solocellfilter, alevin-fry's own knee
    // is the cell call, since the pipeline uses its 'full' matrix as called cells.
    def whitelists = (bc_whitelist instanceof List ? bc_whitelist : [bc_whitelist]).findAll { it }
    def recalls_cells = params.cellfilter_method in ['second_derivative', 'manual_cutoff']
    def knee_reason = !recalls_cells ? "cellfilter_method '${params.cellfilter_method}' uses alevin-fry's own knee as the cell call"
        : !whitelists ? "no barcode whitelist for protocol '${params.protocol}'"
        : ''

    // A whitelist given in parts (one per barcode segment) is expanded into every combination
    // of the parts, in segment order, up to this many barcodes; beyond it the knee is used
    def max_permit_barcodes = 10000000
    """
    echo -e "\\n\\n==================  ALEVIN-FRY =================="
    echo "Sample ID: ${meta}"
    echo "Salmon Index: ${salmon_index}"
    echo "Splici reference: ${splici_index_reference}"
    echo "cDNA read: ${fastq_cDNA}"
    echo "CB/UMI read: ${fastq_BC_UMI}"
    echo "Geometry (bc / umi / read): ${bc_geom} / ${umi_geom} / ${read_geom}"
    echo "Expected orientation: ${expected_ori} (star_soloStrand ${params.star_soloStrand})"
    echo "Resolution: ${params.alevin_resolution}"


    echo -e "\\n\\n-------------  Salmon Alevin -------------------"
    ${mapping_guard} salmon alevin \\
        -i ${salmon_index} \\
        -l A \\
        -1 ${fastq_BC_UMI} \\
        -2 ${fastq_cDNA} \\
        -p ${task.cpus} \\
        --bc-geometry "${bc_geom}" \\
        --umi-geo "${umi_geom}" \\
        --read-geo "${read_geom}" \\
        -o ./${meta.id}_run \\
        --justAlign

    echo -e "\\n\\n-------------  generate permit -------------------"
    permit_args="-k"
    knee_reason="${knee_reason}"
    if [ -z "\$knee_reason" ]; then
        whitelists=(${whitelists.join(' ')})
        n_barcodes=1
        for wl in "\${whitelists[@]}"; do
            n_barcodes=\$(( n_barcodes * \$(zcat -f "\$wl" | tr -d '\\r' | awk 'NF' | wc -l) ))
        done

        if [ "\$n_barcodes" -le ${max_permit_barcodes} ]; then
            # One barcode per line, the segments' combinations in segment order
            zcat -f "\${whitelists[0]}" | tr -d '\\r' | awk 'NF { print \$1 }' > permit_whitelist.txt
            for wl in "\${whitelists[@]:1}"; do
                awk 'NR == FNR { if (NF) seg[++n] = \$1; next }
                     { for (i = 1; i <= n; i++) print \$0 seg[i] }' <(zcat -f "\$wl" | tr -d '\\r') permit_whitelist.txt > permit_whitelist.tmp
                mv permit_whitelist.tmp permit_whitelist.txt
            done
            echo "Unfiltered permit list: \$(wc -l < permit_whitelist.txt) barcodes from \${#whitelists[@]} whitelist file(s)"
            # Every barcode with a read, as in STARsolo's raw matrix: the low end is where the
            # ambient-RNA step finds its empty droplets (alevin-fry's own default is 10 reads)
            permit_args="--unfiltered-pl permit_whitelist.txt --min-reads 1"
        else
            knee_reason="the whitelist expands to \$n_barcodes barcodes, more than ${max_permit_barcodes}"
        fi
    fi
    if [ -n "\$knee_reason" ]; then
        echo "Knee permit list: \$knee_reason"
    fi

    alevin-fry generate-permit-list \\
        -i ./${meta.id}_run \\
        -d ${expected_ori} \\
        --output-dir ./${meta.id}_out_permit \\
        \$permit_args

    echo -e "\\n\\n-------------  collate -------------------"
    alevin-fry collate \\
        -i ./${meta.id}_out_permit \\
        -t ${task.cpus} \\
        -r ./${meta.id}_run

    echo -e "\\n\\n-------------  quant -------------------"
    alevin-fry quant \\
        -m ${splici_index_reference}/*t2g_3col.tsv \\
        -i ./${meta.id}_out_permit \\
        -o ./${meta.id}_counts \\
        -t ${task.cpus} \\
        -r ${params.alevin_resolution}

    rm -f permit_whitelist.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(salmon --version | sed 's/salmon //')
        alevin-fry: \$(alevin-fry --version | sed 's/alevin-fry //')
    END_VERSIONS
    """
}
