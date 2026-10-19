# biodiversitycellatlas/bca_preprocessing: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.4] - 2026-10-19
### Breaking changes
- `ref_gtf_alt` and `perform_featurecounts` are deprecated. Setting either logs a warning and has no effect.
- `geneext_downstream` is new and defaults to `"geneext_only"`, so only GeneExt-annotated runs continue downstream. Use `"both"` for the old behaviour.
- `alevin_resolution` is new and defaults to `"cr-like"` (was `cr-like-em`).
- `validate_params = true` is now the default.
- Relative FASTQ paths in the samplesheet are read from the samplesheet's folder.
- A missing FASTQ, or a `star_index` that is set but missing, now stops the run.

### New features
- **Metacell filtering** (`perform_metacells = true`): Metacell2 metacells, an interactive `filtering_report.html`, and `APPLY_METACELL_FILTER` to rerun with a saved selection.
- **Mapping-rate check:** a sample below `min_mapping_rate` is skipped or aborts the pipeline.
- **`oak_v1` protocol** with i7/i5 demultiplexing, replacing `oak_seq`.
- **RNA velocity** (`perform_velocity`) for STARsolo and alevin-fry, with velocity-ready `.h5ad` files.
- **`test` profile** for tiny test data.
- **Combined `versions.yml`** in the output directory.
- **Nextflow 26.04 support** (also 25.04 and 25.10).
- New parameters: `h5ad_compression`, `geneext_bam_only` and `geneext_subsample_nreads`.
- alevin-fry can now be re-mapped against the GeneExt annotation.

### Improvements
- Read metrics: `CALC_READ_METRICS` replaces `CALC_MT_RRNA` and `CALC_ANTISENSE`. It makes a single BAM pass, adds called-cell and strand-split rows, and `PERCELL_METRICS` no longer reads the BAM.
- Antisense is now taken from STARsolo's `CellReads.stats`.
- alevin-fry uses the same reference as STARsolo, with the same counting rules.
- Input paths are resolved and validated the same way everywhere. Input files are staged into tasks, so relative and cloud paths work.
- GeneExt: Each sample is subsampled before pooling by default.
- Cell calling is one workflow for both mappers.
- Filtering: `MTX_TO_H5AD` now runs last, and the published `.h5ad` carries the doublet and CellSweep annotations.
- Memory: lower peak memory when reading matrices, and `dynamic_memory` now applies.
- `STARSOLO_INDEX` builds one index per run.
- `MAPPING_STATS` and `MULTIQC` stage their inputs instead of scanning `outdir`.
- Submodules are fetched by Nextflow (`recurseSubmodules`), and PAVIAN uses a regular conda environment.
- Clearer dashboard labels and metric names.

### Bug fixes
- alevin-fry: gene IDs containing `-` are handled.
- alevin-fry: `secondderiv_alevin.py` accepts every `alevin_usa_counts` value.
- The published `.h5ad` objects now carry gene symbols.
- Doublet detection is no longer skipped when `ambient_rna_remover = "none"`.
- The doublet filtering summary now names the consensus method.
- `CR_PIPELINE` couldn't find its FASTQs.
- The default Kraken database was never downloaded.
- `STARSOLO_INDEX` failed on Nextflow 26.
- `MAPPING_STATS`: missing alevin-fry and Cell Ranger rows, and failures with a relative `outdir`.
- Dashboard: `PREPARE_DASHBOARD_INPUTS` failed for alevin-fry runs.

### Maintenance
- Channel factories are called as `channel.` (no behaviour change).
- New checks: `alevin_usa`, `mt_rrna_metrics`, `metacell_filtering` and `resource_efficiency: dynamic_memory`.


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

[unreleased]: https://github.com/biodiversitycellatlas/bca_preprocessing/v0.2.0...HEAD
[0.2.0]: https://github.com/biodiversitycellatlas/bca_preprocessing/compare/v0.1.2...v0.2.0
[0.1.2]: https://github.com/biodiversitycellatlas/bca_preprocessing/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/biodiversitycellatlas/bca_preprocessing/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/biodiversitycellatlas/bca_preprocessing/releases/tag/v0.1.0
