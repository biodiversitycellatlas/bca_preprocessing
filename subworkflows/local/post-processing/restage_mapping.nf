//
// Subworkflow with functionality specific to the workflow 'post_mapping_workflow.nf'
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO RE-STAGE PUBLISHED MAPPING RESULTS
        Rebuilds the channels QC_mapping_workflow emits, reading them back from the
        results a previous run published instead of re-running the mappers.

        The analytical run IDs are rebuilt from the samplesheet with the same suffixes
        mapping_workflow.nf appends, so 'mapping_software' and 'perform_geneext' must
        still describe the run being picked up (they are recorded in the previous run's
        pipeline_info/run_config_*.txt). Any analytical run whose published outputs are
        missing or incomplete is reported and skipped, rather than failing the run.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow restage_mapping_workflow {
    take:
        ch_samplesheet      // channel: [meta, fastq_cDNA, fastq_BC_UMI, fastq_indices, samplesheet]

    main:
        def prev_dir = (params.previous_outdir ?: params.outdir).toString()

        // Suffixes mirroring apply_suffix in mapping_workflow.nf, so the directory names
        // resolved here are the ones the previous run published under
        def star_suffixes = []
        if (params.mapping_software == "alevin_subsampled_starsolo") {
            star_suffixes << "_subsampled_starsolo"
        } else if (params.mapping_software in ["starsolo", "both", "alevin_starsolo"]) {
            star_suffixes << "_starsolo"
        }
        if (star_suffixes && params.perform_geneext) {
            star_suffixes << "_geneext_starsolo"
        }

        def alevin_suffixes = params.mapping_software in ["alevin", "both", "alevin_starsolo", "alevin_subsampled_starsolo"]
            ? ["_alevinfry"]
            : []

        // The samplesheet is re-read every run, so expected_cells and manual_cutoff are
        // whatever it says now -- that is what makes the cells re-callable
        def ch_meta = ch_samplesheet.map { row -> row[0] }

        def with_suffix = { meta, suffix ->
            def new_meta = meta.clone()
            new_meta.base_id = meta.id
            new_meta.id = meta.id + suffix
            new_meta.geneext = (suffix == "_geneext_starsolo")
            new_meta
        }

        // Report everything that is missing once per analytical run, rather than failing on
        // the first gap. 'required' is a list of [label, path] pairs.
        def keep_complete = { meta, dir, required ->
            if (!dir.exists()) {
                log.warn "post_mapping: no published results for '${meta.id}' at ${dir}; skipping it"
                return false
            }
            def missing = required.findAll { pair -> !pair[1].exists() }.collect { pair -> pair[0] }
            if (missing) {
                log.warn "post_mapping: '${meta.id}' is missing ${missing.join(', ')} under ${dir}; skipping it"
                return false
            }
            return true
        }

        /*
         * STARsolo
         */
        def ch_star_runs = ch_meta
            .flatMap { meta ->
                star_suffixes.collect { sfx ->
                    [ with_suffix(meta, sfx), file("${prev_dir}/mapping_STARsolo/${meta.id}${sfx}") ]
                }
            }
            .map { meta, dir ->
                def solo = file("${dir}/${meta.id}_Solo.out")
                [ meta, dir, solo, file("${solo}/GeneFull_Ex50pAS") ]
            }
            // Velocyto is deliberately absent from the required list below: a previous run made
            // without 'perform_velocity' has to restage as cleanly as one made with it
            .filter { meta, dir, solo, gene ->
                keep_complete(meta, dir, [
                    [ "${meta.id}_Solo.out",                  solo ],
                    [ "GeneFull_Ex50pAS",                     gene ],
                    [ "GeneFull_Ex50pAS/raw",                 file("${gene}/raw") ],
                    [ "GeneFull_Ex50pAS/Summary.csv",         file("${gene}/Summary.csv") ],
                    [ "GeneFull_Ex50pAS/CellReads.stats",     file("${gene}/CellReads.stats") ],
                    [ "GeneFull_Ex50pAS/UMIperCellSorted.txt", file("${gene}/UMIperCellSorted.txt") ],
                    [ "${meta.id}_Log.out",                   file("${dir}/${meta.id}_Log.out") ],
                    [ "${meta.id}_Log.final.out",             file("${dir}/${meta.id}_Log.final.out") ],
                ])
            }

        ch_star_runs
            .multiMap { meta, dir, solo, gene ->
                solodir:    [meta, solo]
                raw:        [meta, file("${gene}/raw")]
                filtered:   [meta, file("${gene}/filtered")]
                summary:    [meta, file("${gene}/Summary.csv")]
                cellreads:  [meta, file("${gene}/CellReads.stats")]
                umipercell: [meta, file("${gene}/UMIperCellSorted.txt")]
                log_out:    [meta, file("${dir}/${meta.id}_Log.out")]
                log_final:  [meta, file("${dir}/${meta.id}_Log.final.out")]
                bam:        [meta, file("${dir}/${meta.id}_Aligned.sortedByCoord.out.bam")]
                velo_raw:   [meta, file("${solo}/Velocyto/raw")]
            }
            .set { ch_star }

        // All optional upstream: the BAM only exists with star_generateBAM, STARsolo's own
        // filtered matrix is dropped whenever this pipeline re-calls cells, and Velocyto/raw
        // only exists when the previous run set 'perform_velocity'
        def ch_star_bam = ch_star.bam.filter { _meta, bam -> bam.exists() }
        def ch_star_filtered = ch_star.filtered.filter { _meta, dir -> dir.exists() }
        def ch_star_velo_raw = ch_star.velo_raw.filter { _meta, dir -> dir.exists() }

        /*
         * alevin-fry
         */
        def ch_alevin_runs = ch_meta
            .flatMap { meta ->
                alevin_suffixes.collect { sfx ->
                    [ with_suffix(meta, sfx), file("${prev_dir}/mapping_alevin/${meta.id}${sfx}") ]
                }
            }
            .filter { meta, dir ->
                keep_complete(meta, dir, [
                    [ "${meta.id}_run/aux_info/meta_info.json", file("${dir}/${meta.id}_run/aux_info/meta_info.json") ],
                    [ "${meta.id}_counts/quant.json",           file("${dir}/${meta.id}_counts/quant.json") ],
                    [ "${meta.id}_counts/alevin",               file("${dir}/${meta.id}_counts/alevin") ],
                ])
            }

        ch_alevin_runs
            .multiMap { meta, dir ->
                meta_info:  [meta, file("${dir}/${meta.id}_run/aux_info/meta_info.json")]
                quant_json: [meta, file("${dir}/${meta.id}_counts/quant.json")]
                mtx:        [meta, file("${dir}/${meta.id}_counts/alevin")]
                cell_meta:  [meta, file("${dir}/${meta.id}_counts/cell_meta.tsv")]
            }
            .set { ch_alevin }

        def ch_alevin_cell_meta = ch_alevin.cell_meta.filter { _meta, tsv -> tsv.exists() }

        // The analytical runs that were actually found, for the record
        def ch_mapped_metas = ch_star_runs
            .map { meta, _dir, _solo, _gene -> meta }
            .mix(ch_alevin_runs.map { meta, _dir -> meta })

    emit:
        mapped_metas                 = ch_mapped_metas
        star_solodir                 = ch_star.solodir
        starsolo_genefull50_raw      = ch_star.raw
        starsolo_genefull50_filtered = ch_star_filtered
        starsolo_velocyto_raw        = ch_star_velo_raw
        star_summaries               = ch_star.summary
        star_cellreads               = ch_star.cellreads
        star_umipercell              = ch_star.umipercell
        star_log                     = ch_star.log_out
        star_final_log               = ch_star.log_final
        starsolo_bam                 = ch_star_bam
        af_meta_info                 = ch_alevin.meta_info
        af_quant_json                = ch_alevin.quant_json
        af_mtx                       = ch_alevin.mtx
        af_cell_meta                 = ch_alevin_cell_meta
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
