# biodiversitycellatlas/bca_preprocessing: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.4] - 2026-11-19

### Breaking changes

- Nextflow 25.10 or later is required (was 25.04), for nf-schema 2.7.3. Nextflow 26.04 is supported.
- New second-derivative cell caller (`cellfilter_method = "second_derivative"`, the default). The cell boundary is the significant bend of the barcode-rank curve that separates barcodes unlike ambient RNA from those like it, so it also works without a plateau of empty droplets (split-pool protocols). `expected_cells` only breaks ties. Cell counts change; on libraries with shallow cell types they usually rise. `<id>_knee_data.json` records how the boundary was chosen.
- STARsolo's analysed matrix follows `star_soloMultiMappers`, which defaults to `"EM"`, so `raw/UniqueAndMult-EM.mtx` (with fractional counts) is used downstream instead of `matrix.mtx`. Set `"Unique"` for the old behaviour.
- Parse Biosciences: a cell's oligo-dT and random-hexamer (and targeted) round-1 barcodes are counted as one cell, as split-pipe does (`PARSEBIO_COLLAPSE_BC1`). This needs `bc_whitelist_parse_splitwells`; set `parse_collapse_bc1 = false` for the old behaviour.
- Parse: round-1 barcodes are matched by edit distance (up to 2), as split-pipe does, instead of Hamming distance. Barcodes with an `N` no longer match.
- `geneext_downstream` is new and defaults to `"geneext_only"`, so only GeneExt-annotated runs continue downstream. Use `"both"` for the old behaviour.
- `star_extraargs` is added to the protocol's own extra STAR options instead of replacing them. An option given there replaces the protocol's option of the same name.
- `star_soloStrand` defaults to `null`: `"Forward"`, or `"Reverse"` for Ultima CRAMs.
- `perform_trimming` defaults to `null`, the protocol's own default, so `protocol = "seqspec"` no longer trims with fastp. Set `true` to trim.
- `validate_params = true` is now the default.
- `protocol = "ultima_genomics"` is removed: use the `10xv*` protocol of the library chemistry.
- `ref_gtf_addfeature` is no longer added to the mapping annotation. It only feeds the rRNA read metrics.
- `ref_gtf_alt`, `perform_featurecounts` and `demuxafy_sif` are deprecated and have no effect. `code_dir` and `publish_dir_mode` are removed.
- `ambient_rna_remover` is replaced by the boolean `perform_cellsweep` (default `true`), and `gene_annotation*` is renamed `pfam_annotation*`, so the PFAM domain file is not mistaken for the GTF. The old names are deprecated and have no effect.
- The execution trace, report and timeline are written to `execution_reports/` in the launch directory instead of `<outdir>/pipeline_info/`.
- Relative FASTQ paths in the samplesheet are read from the samplesheet's folder.
- A missing FASTQ, or a `star_index` that is set but missing, now stops the run.

### New features

- Metacell filtering (`perform_metacells = true`): Metacell2 metacells, an interactive `filtering_report.html` to filter cells, metacells and genes, and `APPLY_METACELL_FILTER` to rerun with a saved selection. `pfam_annotation` adds protein domains from pfam_scan, eggNOG-mapper or TSV files.
- ScaleBio protocols `scalebio_quantum_rna` and `scalebio_rna_v1_1`, demultiplexed by `SCALEBIO_DEMUX`.
- CRAM input: unaligned Illumina or Ultima Genomics CRAM/BAM files in a new `cram` samplesheet column.
- seqspec layouts for STARsolo and alevin-fry, including split barcodes and whitelists from the spec's onlists.
- Barcode check: a sample whose barcodes do not match the whitelist is stopped before mapping.
- Mapping-rate check: a sample below `min_mapping_rate` is skipped or aborts the pipeline.
- Sample groups: the optional `qc_report` and `filtering_report` samplesheet columns combine samples after mapping, in the dashboard or as one Metacell2 input.
- `oak_v1` protocol with i7/i5 demultiplexing, replacing `oak_seq`.
- RNA velocity (`perform_velocity`) for STARsolo and alevin-fry, with velocity-ready `.h5ad` files.
- Docker, Podman and containers for every module, so every profile runs without conda.
- `test` profile for tiny test data.
- Combined `versions.yml` in the output directory.
- alevin-fry can now be re-mapped against the GeneExt annotation.
- New parameters: `h5ad_compression`, `geneext_bam_only`, `geneext_subsample_nreads`, `alevin_resolution` (default `"cr-like-em"`), `doublet_rate_per_1000` and `star_genomeChrBinNbits`.

### Improvements

- Doublet detection no longer needs the Demuxafy image: its Scrublet, scDblFinder and consensus scripts are ported into `bin/` with their own environments. scDblFinder runs with a fixed seed.
- Read metrics: `CALC_READ_METRICS` replaces `CALC_MT_RRNA` and `CALC_ANTISENSE`. It makes a single BAM pass and adds called-cell and strand-split rows. Antisense is taken from STARsolo's `CellReads.stats`.
- alevin-fry uses the same reference as STARsolo, with the same counting rules.
- Cell calling is one workflow for both mappers, and reads the raw count matrix for both.
- Input paths are resolved and validated the same way everywhere. Input files are staged into tasks, so relative and cloud paths work.
- GeneExt: each sample is subsampled before pooling by default.
- Filtering: the published `.h5ad` carries the doublet and CellSweep annotations.
- Resource requests follow measured use, from local and cluster traces. Memory, time and CPUs scale with each task's input (`dynamic_memory`, `dynamic_time`, `dynamic_cpus`), first attempts are sized for the typical case and retries reach the measured peaks, and single-threaded steps request one core.
- `STARSOLO_ALIGN` runs on 16 threads above 1.5 GiB of FASTQ.
- `STARSOLO_INDEX` builds one index per run and sets `--genomeChrBinNbits` from the genome, which saves memory on fragmented assemblies.
- Faster: `PARSEBIO_CUSTOM_DEMUX` (about 18×), `MERGE_FASTQS` (concatenates gzipped lanes) and `SATURATION_TABLE` (uses all its cores).
- `MAPPING_STATS` and `MULTIQC` stage their inputs instead of scanning `outdir`.
- STARsolo's own `filtered/` matrices are kept under every `cellfilter_method`, for comparison.
- Submodules are fetched by Nextflow (`recurseSubmodules`), and PAVIAN uses a regular conda environment.
- Clearer dashboard labels and metric names. "Genes detected" counts a gene from one UMI.

### Bug fixes

- Doublet detection: a sample keeps one caller's calls when the other fails, instead of losing all. It is no longer skipped when `perform_cellsweep = false`, and the summary names the consensus method.
- `protocol = "seqspec"` failed with STARsolo (duplicate `--soloType`) and with alevin-fry.
- STAR could not read a gzipped `bc_whitelist`.
- Cell Ranger: `CR_PIPELINE` couldn't find its FASTQs and got an extra `versions.yml` argument. The reference is published as `genome/cellranger_ref/`.
- split-pipe: the reference is no longer copied into every sample's output.
- Parse: a well range crossing rows (e.g. `A1-B6`) selected the wrong wells. Malformed wells now stop the run.
- rRNA metrics: whole contigs named in `ref_gtf_addfeature` counted as rRNA instead of the annotated loci. `ref_fasta_addfeature` records that copy `ref_fasta` are dropped.
- alevin-fry: gene IDs containing `-` are handled, and `secondderiv_alevin.py` accepts every `alevin_usa_counts` value.
- The published `.h5ad` objects now carry gene symbols.
- The default Kraken database was never downloaded, and could not be found.
- `STARSOLO_INDEX` failed on Nextflow 26 and timed out on large genomes. `SATURATION_TABLE` timed out on large libraries.
- Thread counts now follow each task's CPU request.
- `MAPPING_STATS`: missing alevin-fry and Cell Ranger rows, and failures with a relative `outdir`.
- Dashboard: `PREPARE_DASHBOARD_INPUTS` failed for alevin-fry runs.
- Containers: some modules could not find their submodule or relied on tools from the host, `MTX_TO_H5AD` failed on numba's cache, and Apptainer pulls failed on Nextflow 25.10.0.
- `PERCELL_METRICS` used a threshold of 0 where `awk` is not GNU awk.

### Maintenance

- nf-schema is updated from 2.3.0 to 2.7.3, so `--help` and the `validation` settings work on every Nextflow version.
- The nf-core lint and pre-commit pass again, and every parameter in `nextflow_schema.json` has a description.
- `nextflow.config` no longer downloads nf-core/demo's config.
- GitHub Actions run on Node.js 24.
- Channel factories are called as `channel.` (no behaviour change).
- New checks: `alevin_usa`, `mt_rrna_metrics`, `addfeature_prepare`, `metacell_filtering`, `parse_collapse`, `parse_demux`, `well_ranges`, `seqspec_layout`, `barcode_check`, `sample_groups`, `containers`, `resource_efficiency` and `dynamic_resources`.
- `tests/bench/star_thread_sweep.sh` benchmarks `STARSOLO_ALIGN` thread counts on SLURM.

## [0.2.0] - 2026-02-24

- Added pavianCore as submodule and increased resources by @bonitavw in #41
- Set limitGenomeGenerateRAM and limitBAMsortRAM as custom config variables by @bonitavw in #42
- Converted to Seqera containers & Fixed permissions by @bonitavw in #43
- MAJOR refractoring of the workflows and subworkflows & added PavianCore by @bonitavw in #44
- Added possibility of specifying kit for split-pipe by @bonitavw in #45
- Remove strict conda dependency by @bonitavw in #46
- Increase resources withLabel process_high_memory by @bonitavw in #47
- Configuration updates and improvements by @bonitavw in #48
- Major refactoring of input and emit channels and fixed env GeneExt by @bonitavw in #49
- Added html dashboard generation and reporting workflow by @bonitavw in #50
- Modify emit channels preprocs subworkflows by @bonitavw in #51

## [0.1.2] - 2025-12-18

- Fix boolean expression to disable generation of BAM during STARsolo by @bonitavw in #36
- Fix STARsolo silent error because of missing pigz by @bonitavw in #37
- Fix empty MAPREADS variable by writing to file + fix seqera container by @bonitavw in #38
- Swithed to cellbender container exclusively by @bonitavw in #39
- Removed conda dependencies from geneext module by @bonitavw in #40

## [0.1.1] - 2025-12-16

- docs: Add troubleshooting guide by @bonitavw in #28
- Allow disabling BAM generation via star_generateBAM config by @bonitavw in #29
- Expanded well-range support (e.g. WT kit) for parse demultiplexing by @bonitavw in #30
- Update README.md & bump manifest version by @bonitavw in #31
- Added Seqera containers to custom modules by @bonitavw in #32
- Split 10x_saturate saturation_table process into two modules by @bonitavw in #33
- Added 10x v4 to seqtech parameters by @bonitavw in #34
- Missing out specification in subworkflow by @bonitavw in #35

## [0.1.0] - 2025-12-11

- Initial release.

[unreleased]: https://github.com/biodiversitycellatlas/bca_preprocessing/compare/v0.2.4...HEAD
[0.2.4]: https://github.com/biodiversitycellatlas/bca_preprocessing/compare/v0.2.3...v0.2.4
[0.2.0]: https://github.com/biodiversitycellatlas/bca_preprocessing/compare/v0.1.2...v0.2.0
[0.1.2]: https://github.com/biodiversitycellatlas/bca_preprocessing/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/biodiversitycellatlas/bca_preprocessing/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/biodiversitycellatlas/bca_preprocessing/releases/tag/v0.1.0
