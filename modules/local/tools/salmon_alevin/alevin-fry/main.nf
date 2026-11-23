process ALEVIN_FRY {
    publishDir path: { "${params.outdir}/mapping_alevin/${meta.id}" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${meta.id}"


    // Time follows the FASTQs and memory follows the salmon index size.
    // Salmon's mapping stops gaining at about 8 threads (66 s for 30 M pairs at 8, 12 and 16), so more only adds cost.
    cpus   8
    time   { BcaResources.scaledTime(
        params.dynamic_time?.ALEVIN_FRY, [fastq_cDNA, fastq_BC_UMI], task.attempt, 4) }
    memory { BcaResources.scaledMemory(
        params.dynamic_memory?.ALEVIN_FRY, [fastq_cDNA, fastq_BC_UMI], task.attempt, 64, salmon_index) }

    // Exit 42 from bin/mapping_rate_guard.sh (mapping rate below params.min_mapping_rate)
    // aborts the pipeline or drops the sample; other exit codes follow conf/base.config.
    errorStrategy { MappingRateCheck.errorStrategy(task, params) }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/alevin-fry_piscem_salmon_simpleaf_pruned:c71cfb476b003414' :
        'community.wave.seqera.io/library/alevin-fry_piscem_salmon_simpleaf_pruned:0c7b7b78674d62a2' }"

    input:
    tuple val(meta), path(fastq_cDNA), path(fastq_BC_UMI), path(fastq_indices), path(input_file)
    path(bc_whitelist, stageAs: 'whitelist_?/*')        // in barcode-segment order; empty without a whitelist
    path(seqspec_layout, stageAs: 'seqspec_layout/*')   // SEQSPEC_LAYOUT's layout for protocol 'seqspec', empty otherwise
    path(splici_index_reference)
    path(salmon_index)
    val mapping_flag_dir                                // MappingRateCheck.flagDir(workflow), as a string

    output:
    tuple val(meta), path("${meta.id}_*"), emit: mapping_files
    tuple val(meta), path("${meta.id}_run/aux_info/meta_info.json"),    emit: af_meta_info
    tuple val(meta), path("${meta.id}_counts/quant.json"),              emit: af_quant_json
    tuple val(meta), path("${meta.id}_counts/alevin"),                  emit: af_mtx
    tuple val(meta), path("${meta.id}_counts/cell_meta.tsv"),           emit: af_cell_meta, optional: true
    path "versions.yml",                                                emit: versions

    script:
    // Retrieve settings from custom parameters if set, otherwise from conf/seqtech_parameters.config.
    // Protocol 'seqspec' fills the ones not set from the layout derived from the seqspec file
    // (SEQSPEC_LAYOUT), in the script below.
    def seqtech = params.seqtech_parameters[params.protocol]
    def bc_geom = params.alevin_bc_geometry ?: seqtech.alevin_bc_geometry
    def umi_geom = params.alevin_umi_geometry ?: seqtech.alevin_umi_geometry
    def read_geom = params.alevin_read_geometry ?: seqtech.alevin_read_geometry

    // The seqspec layout replaces the protocol's geometries, but not geometries set in the run
    def geom_from_layout = seqspec_layout ? [
        bc:   params.alevin_bc_geometry ? '' : 'yes',
        umi:  params.alevin_umi_geometry ? '' : 'yes',
        read: params.alevin_read_geometry ? '' : 'yes'
    ] : [bc: '', umi: '', read: '']

    if (!seqspec_layout && (!bc_geom || !umi_geom || !read_geom)) {
        error "No alevin geometry defined for protocol '${params.protocol}'. Set 'alevin_bc_geometry', 'alevin_umi_geometry' and 'alevin_read_geometry' in the configuration file."
    }

    // Mapping-rate check: salmon runs under bin/mapping_rate_guard.sh, which reads
    // percent_mapped from meta_info.json before the alevin-fry steps, and cancels salmon
    // early if STARsolo already found this sample below the threshold.
    def mapping_guard = MappingRateCheck.enabled(params.skip_mapping_rate_check, params.min_mapping_rate)
        ? "mapping_rate_guard.sh run --mapper salmon" +
          " --flag-dir ${mapping_flag_dir} --sample ${meta.base_id ?: meta.id} --label ${meta.id}" +
          " --min ${params.min_mapping_rate} --poll ${params.mapping_rate_poll_secs ?: 60}" +
          " --report ./${meta.id}_run/aux_info/meta_info.json --"
        : ''

    // Count reads on the strand STARsolo counts (--soloStrand), so both mappers see the same molecules
    def star_soloStrand = StarSettings.strand(params.star_soloStrand, params.cram_type, meta.input_format)
    def expected_ori = [Forward: 'fw', Reverse: 'rc', Unstranded: 'both'][star_soloStrand]
    if (!expected_ori) {
        error "star_soloStrand '${star_soloStrand}' has no alevin-fry orientation; use Forward, Reverse or Unstranded."
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

    // sciRNAseq3 (alevin_cb_segments in conf/seqtech_parameters.config): its 40 nt cell barcode
    // is longer than salmon alevin takes, so salmon reads read 1 rewritten as a code of the
    // cell's place in the sample's barcode lists, then the UMI (bin/sci_cb_code.sh, through a
    // process substitution). The barcodes are corrected already, so the permit list is the
    // codes that occur, and the matrix's barcodes are decoded to STARsolo's form at the end.
    // Not with geometries set in the run, which describe the reads as they are.
    def cb_code = seqtech.alevin_cb_segments && !params.alevin_bc_geometry && !params.alevin_umi_geometry
    def cb_args = "--lists ${whitelists.join(' ')} --digits ${seqtech.alevin_cb_digits}"
    def read1 = !cb_code ? fastq_BC_UMI
        : "<(sci_cb_code.sh encode ${cb_args} --segments ${seqtech.alevin_cb_segments} --umi ${seqtech.alevin_cb_umi} --r1 ${fastq_BC_UMI} --out - --permit sci_cb_permit.txt)"

    // The encoder has finished, so the permit list is complete; a failed encoder stops the task
    def cb_wait = cb_code ? 'wait $!' : ''
    def permit_lists = cb_code ? ['sci_cb_permit.txt'] : whitelists
    def rows = "./${meta.id}_counts/alevin/quants_mat_rows.txt"

    def cb_decode = !cb_code ? '' : """
    echo -e "\\n\\n-------------  Decode cell barcodes -------------------"
    # The matrix's barcodes back to the cell barcode, its segments joined with '_' as STARsolo
    # writes them. featureDump.txt and the permit list files keep the codes, as alevinQC joins
    # them with each other; cell_barcode_codes.tsv pairs each code with its barcode.
    sci_cb_code.sh decode ${cb_args} --with-code < ${rows} > ./${meta.id}_counts/cell_barcode_codes.tsv
    cut -f2 ./${meta.id}_counts/cell_barcode_codes.tsv > ${rows}
    rm -f sci_cb_permit.txt
    """
    """
    echo -e "\\n\\n==================  ALEVIN-FRY =================="
    echo "Sample ID: ${meta}"
    echo "Salmon Index: ${salmon_index}"
    echo "Splici reference: ${splici_index_reference}"
    echo "cDNA read: ${fastq_cDNA}"
    echo "CB/UMI read: ${fastq_BC_UMI}"
    BC_GEOM="${bc_geom ?: ''}"
    UMI_GEOM="${umi_geom ?: ''}"
    READ_GEOM="${read_geom ?: ''}"
    if [[ -n "${seqspec_layout}" ]]; then
        # Geometry derived from the seqspec file (SEQSPEC_LAYOUT, bin/seqspec_layout.py)
        source ${seqspec_layout}
        if [[ "${geom_from_layout.bc}" == "yes" ]]; then BC_GEOM="\${SEQSPEC_ALEVIN_BC}"; fi
        if [[ "${geom_from_layout.umi}" == "yes" ]]; then UMI_GEOM="\${SEQSPEC_ALEVIN_UMI}"; fi
        if [[ "${geom_from_layout.read}" == "yes" ]]; then READ_GEOM="\${SEQSPEC_ALEVIN_READ}"; fi
    fi
    if [[ -z "\${BC_GEOM}" || -z "\${UMI_GEOM}" || -z "\${READ_GEOM}" ]]; then
        echo "ERROR: no alevin geometry for protocol '${params.protocol}' (bc '\${BC_GEOM}', umi '\${UMI_GEOM}', read '\${READ_GEOM}'). Set 'alevin_bc_geometry', 'alevin_umi_geometry' and 'alevin_read_geometry'." >&2
        exit 1
    fi
    echo "Geometry (bc / umi / read): \${BC_GEOM} / \${UMI_GEOM} / \${READ_GEOM}"
    echo "Expected orientation: ${expected_ori} (star_soloStrand ${star_soloStrand})"
    echo "Resolution: ${params.alevin_resolution}"


    echo -e "\\n\\n-------------  Salmon Alevin -------------------"
    ${mapping_guard} salmon alevin \\
        -i ${salmon_index} \\
        -l A \\
        -1 ${read1} \\
        -2 ${fastq_cDNA} \\
        -p ${task.cpus} \\
        --bc-geometry "\${BC_GEOM}" \\
        --umi-geo "\${UMI_GEOM}" \\
        --read-geo "\${READ_GEOM}" \\
        -o ./${meta.id}_run \\
        --justAlign
${cb_wait}
    echo -e "\\n\\n-------------  generate permit -------------------"
    permit_args="-k"
    knee_reason="${knee_reason}"
    if [ -z "\$knee_reason" ]; then
        whitelists=(${permit_lists.join(' ')})
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
${cb_decode}
    rm -f permit_whitelist.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(salmon --version | sed 's/salmon //')
        alevin-fry: \$(alevin-fry --version | sed 's/alevin-fry //')
    END_VERSIONS
    """
}
