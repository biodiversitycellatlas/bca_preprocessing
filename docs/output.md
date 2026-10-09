# Output

This page describes every directory and file the pipeline can produce, and how to interpret
them. Which directories appear depends on the options you enabled — see
[Configuration parameters](CONFIGURATION_PARAMETERS.md).

## Table of Contents

1. [Analytical runs and directory naming](#analytical-runs-and-directory-naming)
2. [Directory structure](#directory-structure)
3. [Run information](#run-information) — `pipeline_info/`
4. [Reference files](#reference-files) — `genome/`
5. [Read preparation](#read-preparation) — `demultiplex/`, `marsseq_reformat/`, `fastp/`
6. [Quality control](#quality-control) — `fastqc/`
7. [Mapping: STARsolo](#mapping-starsolo) — `mapping_STARsolo/`
8. [Mapping: alevin-fry](#mapping-alevin-fry) — `mapping_alevin/`
9. [Cell calling](#cell-calling) — `filtered_secondderiv/`
10. [Intronic counts and RNA velocity](#intronic-counts-and-rna-velocity) — `Velocyto/`, `anndata/<id>/velocity/`
11. [AnnData conversion](#anndata-conversion) — `anndata/`
12. [Ambient RNA removal](#ambient-rna-removal) — `cellsweep/`
13. [Doublet detection](#doublet-detection) — `doublet_filtering/`
14. [Metacells and filtering](#metacells-and-filtering) — `metacells/`, `filtering_report.html`, `metacell_filtering/`
15. [Optional analyses](#optional-analyses) — `saturation/`, `rRNA_mtDNA/`, `gene_ext/`, `kraken/`
16. [Summary reports](#summary-reports) — `summary_results/`, `dashboard.html`
17. [External pipelines](#external-pipelines)
18. [Which matrix should I use?](#which-matrix-should-i-use)

---

## Analytical runs and directory naming

Before reading the tree below, one concept explains most of it.

A single sample can be quantified more than once in one pipeline run — with STARsolo, with
alevin-fry, and again against a GeneExt-extended annotation. Each of these is an **analytical
run** with its own identifier, formed by appending a suffix to the sample name:

| Suffix | Produced when | Example directory name |
| ------ | ------------- | ---------------------- |
| `_starsolo` | STARsolo mapping (`mapping_software = "starsolo"` or `"both"`) | `sampleA_starsolo` |
| `_alevinfry` | alevin-fry mapping (`mapping_software = "alevin"` or `"both"`) | `sampleA_alevinfry` |
| `_geneext_starsolo` | STARsolo against the extended annotation (`perform_geneext = true`) | `sampleA_geneext_starsolo` |
| `_geneext_alevinfry` | alevin-fry against the extended annotation (`perform_geneext = true`, with a `mapping_software` that runs alevin-fry) | `sampleA_geneext_alevinfry` |
| `_subsampled_starsolo` | STARsolo on subsampled reads (`mapping_software = "alevin_subsampled_starsolo"`) | `sampleA_subsampled_starsolo` |
| `_geneext_subsampled_starsolo` | STARsolo on the same subsampled reads, against the extended annotation (`alevin_subsampled_starsolo` with `perform_geneext = true`) | `sampleA_geneext_subsampled_starsolo` |

Per-sample subdirectories throughout the output are named by the **analytical run**, not by the
sample. So a run of `sampleA` with `mapping_software = "both"`, `perform_geneext = true` and
`geneext_downstream = "both"` produces `sampleA_starsolo`, `sampleA_geneext_starsolo`,
`sampleA_alevinfry` and `sampleA_geneext_alevinfry` side by side. Nothing overwrites anything,
and the dashboard shows all of them together. With the default `geneext_downstream =
"geneext_only"` only the two `_geneext_` runs are made; `mapping_STARsolo/sampleA_starsolo/`
then holds just the alignment GeneExt read (BAM and logs), and appears nowhere else. See
[Which runs GeneExt adds](CONFIGURATION_PARAMETERS.md#which-runs-geneext-adds).

Both mappers are built on the same reference: `ref_fasta` and `ref_gtf` merged with
`ref_fasta_addfeature` and `ref_gtf_addfeature` when those are set, and the GeneExt annotation for
the `_geneext_` runs. 

## Directory structure

```
output_directory/
├── dashboard.html                  # Portable HTML dashboard — start here
├── filtering_report.html           # (Optional) Interactive metacell filtering report
├── versions.yml                    # Software versions of every process that ran
│
├── pipeline_info/                  # Run configuration, samplesheet, execution traces
├── genome/                         # Generated reference indices
│
├── demultiplex/                    # (Protocol-dependent) demultiplexed FASTQ files
├── marsseq_reformat/               # (MARS-seq only) reconstructed reads
├── fastp/                          # (Protocol-dependent) trimmed FASTQ files
├── fastqc/                         # Quality control reports for raw FASTQ files
│
├── mapping_STARsolo/               # STARsolo alignments, matrices and statistics
├── mapping_alevin/                 # (Optional) alevin-fry matrices and AlevinQC
├── anndata/                        # Count matrices converted to .h5ad
│   └── <analytical_run>/velocity/  # (Optional) Spliced/unspliced layers for RNA velocity
│
├── cellsweep/                      # (Optional) CellSweep ambient RNA removal
├── doublet_filtering/              # (Optional) Scrublet, scDblFinder and consensus calls
├── metacells/                      # (Optional) Metacell2 metacells and their QC summaries
├── metacell_filtering/             # (Optional) Filtered UMI matrices from the report's selection
│
├── saturation/                     # (Optional) Sequencing saturation analysis
├── rRNA_mtDNA/                     # rRNA, mtDNA and antisense reads (runs with a BAM)
├── gene_ext/                       # (Optional) GeneExt extended annotation, report and log
├── kraken/                         # (Optional) Taxonomic classification of unmapped reads
│
├── summary_results/                # MultiQC report, mapping_stats.tsv, per-cell metrics
├── containers/                     # (Optional) Downloaded Demuxafy image
│
├── CellRanger_pipeline/            # (Optional) External Cell Ranger outputs
├── ParseBio_pipeline/              # (Optional) External split-pipe outputs
└── BDrhapsody_pipeline/            # (Optional) External BD Rhapsody outputs
```

---

## Run information

<details markdown="1">
<summary><code>pipeline_info/</code></summary>

| File | Description |
| ---- | ----------- |
| `run_config_<timestamp>.txt` | The full resolved configuration of the run — every parameter as it was actually used. |
| `samplesheet_<timestamp>.csv` | A copy of the samplesheet used. |
| `execution_report_<timestamp>.html` | Nextflow's execution report: resource usage per process. |
| `execution_timeline_<timestamp>.html` | Nextflow's timeline of task execution. |
| `execution_trace_<timestamp>.txt` | Per-task trace: CPU, memory, I/O, exit status, retries. |
| `pipeline_dag_<timestamp>.html` | The workflow DAG. |
| `mapping_rate_check.tsv` | One row per sample with the verdict of the [mapping-rate check](CONFIGURATION_PARAMETERS.md#mapping-rate-check): mapper, stage (`mid-run` or `final`), reads processed, percentage mapped and threshold. Only written when the check ran, and overwritten by each run. |

Files are timestamp-suffixed rather than overwritten, so repeated runs into the same output
directory accumulate rather than replace.

> [!IMPORTANT]
> **Keep the `execution_trace_*.txt` files.** They are the input to
> `bin/resource_efficiency.py`, which reports what your cluster actually used per process and
> recommends adjusted resource requests. See [Resource tuning](RESOURCE_TUNING.md).

The run configuration and samplesheet are also embedded in `dashboard.html`, so a dashboard
found on its own remains self-describing.

The software versions are in `versions.yml` at the top of the output directory: one entry
per process that ran, listing the tools it used, plus a `Workflow` entry with the pipeline and
Nextflow versions. It is rewritten at the end of every run, so it describes the most recent one.

</details>

## Reference files

<details markdown="1">
<summary><code>genome/</code></summary>

Generated reference files, reusable across runs on the same species.

| Path | Produced by | Description |
| ---- | ----------- | ----------- |
| `star_index_<gtf_name>/` | `STARSOLO_INDEX` | STAR genome index, built once per run (`star_index_<gtf_name>_geneext/` for the GeneExt remap). |
| `splici_index/` | `SALMON_SPLICI` | Spliced + intronic (splici) reference (`*_splici*.fa`) and its three-column transcript-to-gene map (`*_t2g_3col.tsv`), used by alevin-fry (`splici_index_geneext/` for the GeneExt remap). |
| `salmon_index/` | `SALMON_INDEX` | Salmon index built on the splici reference (`salmon_index_geneext/` for the GeneExt remap). |
| `ref.fasta`, `ref.gtf` | `MERGE_REF_FASTA`, `MERGE_REF_GTF` | Reference merged with `ref_fasta_addfeature` / `ref_gtf_addfeature`, when those are set. |
| `ref_geneext.gtf` | `MERGE_REF_GTF_GENEEXT` | The GeneExt annotation merged with `ref_gtf_addfeature`, used by the GeneExt remap (`perform_geneext = true`). |
| `cellranger_ref/` | `CELLRANGER_MKREF` | Cell Ranger reference, only when `perform_cellranger = true`. |

> [!TIP]
> Building a STAR index is expensive. Point `star_index` at an existing index to skip
> `STARSOLO_INDEX` on subsequent runs.

</details>

## Read preparation

Which of these appear depends entirely on the protocol. See
[Protocol-specific steps & whitelists](PROTOCOLS_AND_WHITELISTS.md) for what each protocol
triggers.

<details markdown="1">
<summary><code>demultiplex/</code>, <code>marsseq_reformat/</code>, <code>fastp/</code></summary>

**`demultiplex/demux_custom/<sample>/`** — Parse Biosciences. Reads assigned to samples by
matching the first barcode against the split-well definition within a Hamming distance of 2.

**`demultiplex/demux_spipe/`** — Parse Biosciences, when demultiplexing with split-pipe instead
of the built-in demultiplexer.

**`demultiplex/<sample>/`** — OAK, written by `OAK_DEMUX`: the reads of one aliquot, selected
on its i7 (optional) and i5 index, as `<sample>_S1_L001_{R1,R2,I1,I2}_001.fastq.gz`, plus
`<sample>_oak_demux.log` with the read counts of the input and after each index step.

**`demux_reads/`** — sci-RNA-seq3, written by `SCIROCKET_DEMUX`:

| File | Description |
| ---- | ----------- |
| `<sample>_R1.fastq.gz` | The synthetic 48 nt barcode read (p7, p5, ligation, RT barcodes + UMI). |
| `<sample>_R2.fastq.gz` | cDNA. |
| `<sample>_R1_discarded*`, `<sample>_R2_discarded*` | Reads that could not be assigned. |
| `<sample>_whitelist_p5*`, `_p7*`, `_ligation*`, `_rt*` | Per-run barcode whitelists, generated from the observed data and passed straight to the mapper. |

**`marsseq_reformat/<sample>/`** — MARS-seq. Both reads rebuilt so that the batch barcode is
prepended to the cell barcode. The resulting barcode read is 14 nt for v1 (CB 1–10, UMI 11–14)
and 19 nt for v2 (CB 1–11, UMI 12–19).

**`fastp/`** — Trimmed FASTQ files (`trimmed_<original_name>`), plus fastp's own HTML and JSON
reports.

</details>

## Quality control

<details markdown="1">
<summary><code>fastqc/</code></summary>

Standard FastQC output for the raw FASTQ files: one `*_fastqc.html` and one `*_fastqc.zip` per
input file. The zips are also aggregated into the MultiQC report in `summary_results/`.

These are run on the **raw** input files, before any protocol-specific read manipulation.

</details>

## Mapping: STARsolo

<details markdown="1">
<summary><code>mapping_STARsolo/&lt;analytical_run&gt;/</code></summary>

| File | Description |
| ---- | ----------- |
| `<id>_Solo.out/` | STARsolo output tree — see below. |
| `<id>_Log.final.out` | Alignment summary: input reads, uniquely mapped %, multimapping %, unmapped % by reason. The source of most of the dashboard's Mapping tab. |
| `<id>_Log.out` | Verbose run log, including the resolved parameters STAR ran with. |
| `<id>_Log.progress.out` | Progress log. |
| `<id>_Aligned.sortedByCoord.out.bam` | Coordinate-sorted alignments, present when `star_generateBAM = true`. Includes unmapped reads (`--outSAMunmapped Within`), which is what makes the Kraken step possible. Under `run_method = "geneext_only"` with the default `geneext_bam_only` it is written for GeneExt alone: no `CR`/`UR`/`CB`/`UB` tags and no unmapped reads. |
| `<id>_Aligned.sortedByCoord.out.bam.bai` | BAM index. |
| `<id>_SJ.out.tab` | Splice junctions. |

Inside `<id>_Solo.out/` there is one directory per requested feature. The first two are produced
by default (`--soloFeatures "Gene GeneFull_Ex50pAS"`):

| Feature | Counts |
| ------- | ------ |
| `Gene/` | Exonic reads only. |
| `GeneFull_Ex50pAS/` | Exonic **and** intronic reads. This is the feature the pipeline uses downstream, and the one comparable to alevin-fry's collapsed USA matrix. |
| `Velocyto/` | Spliced, unspliced and ambiguous counts — the intronic matrix, for RNA velocity. **Only present when `perform_velocity = true`.** See [Intronic counts and RNA velocity](#intronic-counts-and-rna-velocity). |

`Gene/` and `GeneFull_Ex50pAS/` each contain:

| Path | Description |
| ---- | ----------- |
| `raw/` | Unfiltered matrix over all barcodes: `matrix.mtx`, `features.tsv`, `barcodes.tsv`. |
| `filtered/` | STARsolo's own cell-called matrix. **Only present when `cellfilter_method = "star_solocellfilter"`** — under the other methods it is removed after the summary is read, to avoid two competing "filtered" matrices. |
| `Summary.csv` | Per-sample summary: reads, saturation, cells, median UMIs and genes per cell. Always written from STARsolo's own filter, even when a different cell-calling method is used. |
| `CellReads.stats` | Per-barcode read statistics. |
| `UMIperCellSorted.txt` | UMI count per barcode, sorted descending. This is the barcode-rank curve, and the input to second-derivative cell calling. |
| `Features.stats` | Feature-level assignment statistics. |

Only `GeneFull_Ex50pAS/` gets a `filtered_secondderiv/` directory — it is the feature the
pipeline calls cells on. See [Cell calling](#cell-calling).

</details>

## Mapping: alevin-fry

<details markdown="1">
<summary><code>mapping_alevin/&lt;analytical_run&gt;/</code></summary>

| Path | Description |
| ---- | ----------- |
| `<id>_run/` | `salmon alevin` mapping output. |
| `<id>_run/aux_info/meta_info.json` | Mapping-rate statistics. |
| `<id>_out_permit/` | alevin-fry permit list: every whitelisted barcode with a read (`--unfiltered-pl`), or alevin-fry's knee where that does not apply. See [alevin-fry settings](CONFIGURATION_PARAMETERS.md#alevin-fry-settings). |
| `<id>_counts/` | alevin-fry quantification output. |
| `<id>_counts/alevin/` | The count matrix in USA mode: `quants_mat.mtx` (cells × columns), `quants_mat_rows.txt` (barcodes) and `quants_mat_cols.txt` (the gene IDs, then the same IDs with `-U`, then with `-A`). |
| `<id>_counts/quant.json` | Quantification summary. |
| `<id>_counts/cell_meta.tsv` | Per-cell summary statistics. |
| `<id>_counts/alevin/filtered_secondderiv/` | The pipeline's cell-called matrix. See [Cell calling](#cell-calling). |
| `<id>_alevinFry_QC.html` | AlevinQC report. |
| `gene_level_matrix/<datatype>/gene_level/` | The USA matrix collapsed to one count per gene — see below. |

> [!IMPORTANT]
> alevin-fry runs in **USA mode**, so `quants_mat.mtx` carries *three* columns per gene:
> spliced (S), unspliced (U) and ambiguous (A). `COLLAPSE_ALEVIN_USA` sums the blocks selected
> by `alevin_usa_counts` into a single count per gene and writes the result to
> `gene_level_matrix/`. **Use the collapsed matrix**, not the raw USA matrix, unless you
> specifically want the splicing breakdown.
>
> | `alevin_usa_counts` | Sums | STARsolo equivalent |
> | ------------------- | ---- | ------------------- |
> | `SUA` (default) | S + U + A | `GeneFull_Ex50pAS` |
> | `SA` | S + A | `Gene` |
> | `S` | S only | *(none)* |
> | `UA` | U + A | *(none)* |
> | `U` | U only | `Velocyto/` unspliced |
>
> `SUA` is the only setting under which the two mappers' main matrices are directly comparable.
> You do not need to set `alevin_usa_counts = "U"` to get an intronic matrix: with
> `perform_velocity = true` the unspliced block is extracted into its own
> `gene_level_matrix/<datatype>_unspliced/` directory alongside the main matrix. See
> [Intronic counts and RNA velocity](#intronic-counts-and-rna-velocity).

</details>

## Cell calling

<details markdown="1">
<summary><code>filtered_secondderiv/</code> (inside each mapping directory)</summary>

Present when `cellfilter_method = "second_derivative"` (the default) or `"manual_cutoff"`.
The location differs per mapper:

- STARsolo: `mapping_STARsolo/<id>/<id>_Solo.out/GeneFull_Ex50pAS/filtered_secondderiv/`
- alevin-fry: `mapping_alevin/<id>/<id>_counts/alevin/filtered_secondderiv/`

| File | Description |
| ---- | ----------- |
| `filtered/` | The cell-called count matrix: `matrix.mtx`, `features.tsv`, `barcodes.tsv`. **This is the matrix most analyses should start from.** |
| `<id>_cutoff.txt` | The UMI threshold that was applied. |
| `<id>_knee_data.json` | The barcode-rank curve and its second derivative, as plotted in the dashboard's Cell Calling tab. |
| `<id>_secondderiv_statistics.json` | Statistics recomputed against this cutoff: number of cells, UMI cutoff, mean reads per cell, median UMIs and genes per cell, total genes detected, noise (% of UMIs in non-cell barcodes), sequencing saturation, and reads required to reach `saturation_target`. |

The cutoff is placed at the most negative second derivative of the log-log cumulative UMI
curve, searched between 0.1× and 10× the sample's `expected_cells`. The retained cell count is
**not** forced towards `expected_cells`. For the full rationale, see
[Cell calling](CONFIGURATION_PARAMETERS.md#cell-calling).

> [!NOTE]
> `Summary.csv` in the STARsolo output is always written from STARsolo's own filter, even when
> second-derivative calling is used. When the two disagree, `<id>_secondderiv_statistics.json`
> is the one that describes the matrix in `filtered/`.

</details>

## Intronic counts and RNA velocity

<details markdown="1">
<summary><code>Velocyto/</code>, <code>gene_level_matrix/*_unspliced/</code>, <code>anndata/&lt;analytical_run&gt;/velocity/</code></summary>

Produced only when `perform_velocity = true`. Both mappers can split their counts by splicing
status, which is what RNA velocity needs — the unspliced matrix is the intronic count matrix.

**STARsolo** — `mapping_STARsolo/<id>/<id>_Solo.out/Velocyto/`

| Path | Description |
| ---- | ----------- |
| `raw/spliced.mtx` | Counts from UMIs compatible only with the mature mRNA. |
| `raw/unspliced.mtx` | Counts from UMIs overlapping an intron or an exon–intron boundary — **the intronic matrix**. |
| `raw/ambiguous.mtx` | Counts from UMIs that do not discriminate between the two. |
| `raw/barcodes.tsv`, `raw/features.tsv` | The single barcode and feature axis all three matrices share. |
| `filtered_secondderiv/filtered/` | The same three matrices subset to the called cells. Written under **every** `cellfilter_method`, not only the second-derivative ones. |

The three matrices are mutually exclusive and additive, so no count is double-reported.

There is no `Velocyto/filtered/` directory. STARsolo does write one, but it stops after the
first count column, so it holds `spliced.mtx` alone with no unspliced or ambiguous
counterpart — it is removed, and `filtered_secondderiv/filtered/` is derived from `raw/`
instead, whichever cell-calling method is in force.

**alevin-fry** — `mapping_alevin/<id>/gene_level_matrix/<datatype>_unspliced/`

alevin-fry already quantifies against a splici reference in USA mode, so its unspliced block
is extracted with the same collapse step that produces the main matrix, run with `--counts U`.
This is the counterpart of STARsolo's `unspliced.mtx`.

**Both mappers** — `anndata/<analytical_run>/velocity/<id>_<datatype>_velocity.h5ad`

A velocity-ready AnnData object of the called cells (`filtered`; alevin-fry's `full` under
`cellfilter_method = "star_solocellfilter"`, where that is its cell call), which is what
scVelo, velocyto and CellRank actually expect. Each mapper's layers follow its own
convention, recorded in `adata.uns["velocity_layers"]`:

| Field | STARsolo | alevin-fry |
| ----- | -------- | ---------- |
| `adata.layers["spliced"]` | Velocyto spliced | S + A, as alevin-fry's velocity tutorial and pyroe define it |
| `adata.layers["unspliced"]` | Velocyto unspliced (the intronic matrix) | U |
| `adata.layers["ambiguous"]` | Velocyto ambiguous | A, for reference: already inside `spliced` |
| `adata.X` | The sum of the three layers | S + U + A |

Compare the two mappers' `spliced` and `unspliced` layers; their ambiguous classes are
defined differently.

> [!IMPORTANT]
> The velocity matrices take their cell set from the `GeneFull_Ex50pAS` cell call rather than a
> cutoff of their own, so `Velocyto/filtered_secondderiv/filtered/barcodes.tsv` is identical to
> `GeneFull_Ex50pAS/filtered_secondderiv/filtered/barcodes.tsv`. All published matrices for a
> sample therefore describe the same cells.

> [!NOTE]
> `adata.X` is **not** expected to equal the `GeneFull_Ex50pAS` matrix. The two features assign
> reads to genes by different rules, so the totals are close but not identical. For the same
> reason, subtracting `Gene` from `GeneFull_Ex50pAS` is *not* a substitute for this output: the
> two are separate counting passes, not a superset and a subset, so the difference is negative
> for some gene/cell pairs and is not an intron-only count.

> [!WARNING]
> These counts cannot be recovered from an existing run — they require a fresh STARsolo pass.
> `run_method = "post_mapping"` will re-stage them from a previous run that set
> `perform_velocity`, but cannot add them to one that did not.

</details>

## AnnData conversion

<details markdown="1">
<summary><code>anndata/&lt;analytical_run&gt;/&lt;datatype&gt;/</code></summary>

Count matrices converted to `.h5ad` for downstream analysis in scanpy. `<datatype>` is one of:

| `datatype` | Matrix |
| ---------- | ------ |
| `raw` | STARsolo unfiltered matrix (all barcodes). |
| `filtered` | Cell-called matrix. |
| `full` | alevin-fry full matrix (the alevin-fry counterpart of `raw`): every barcode on alevin-fry's permit list. |

These are the **last** step of the filtering workflow, not the first: ambient RNA removal
(which uses `raw`/`full`, since it models the empty droplets that cell calling removes) and
doublet detection (which uses the cell-called matrix) both work on the `mtx`/`barcodes`/`features`
triplet straight out of the mapper, and their results are merged in here. So each object
carries not only the counts but everything that was called on them:

| Field | Written when |
|-------|--------------|
| `X` | always — the counts, with the consensus doublets already removed under `perform_doublet_filtering` |
| `obs["sample_id"]` | always |
| `obs["doublet_status"]` | `perform_doublet_detection` — `singlet` / `doublet`. The empty droplets of a `raw`/`full` matrix were never evaluated and read `singlet` |
| `layers["cellsweep"]` | `ambient_rna_remover = "cellsweep"`, `raw`/`full` only — the denoised counts, aligned to `X` |
| `obs["is_empty"]`, `obs["celltype"]`, `obs["alpha_hat"]` | as above — empty-droplet call, Leiden grouping, per-cell ambient fraction |
| `var["ambient_hat"]` | as above — per-gene ambient profile |

The velocity objects in `anndata/<analytical_run>/velocity/` are a separate, terminal output —
they are deliberately not fed to ambient RNA removal or doublet detection, since neither is
meaningful on unspliced counts. See
[Intronic counts and RNA velocity](#intronic-counts-and-rna-velocity).

</details>

## Ambient RNA removal

Controlled by `ambient_rna_remover` (`"cellsweep"` by default, or `"none"`).
It operates on the **unfiltered** matrix.

<details markdown="1">
<summary><code>cellsweep/&lt;analytical_run&gt;/</code></summary>

| File | Description |
| ---- | ----------- |
| `<id>_cs_filtered.h5ad` | Decontaminated matrix, cells only. |
| `<id>_cs_full.h5ad` | Decontaminated matrix, all barcodes. |
| `<id>_ambient_hat_histogram.png` | Distribution of estimated ambient contribution per cell. |
| `<id>_top_ambient_genes.csv` | Genes contributing most to the ambient profile — worth inspecting, since a dominant ambient gene often points at a specific problem. |
| `<id>_umap_comparison.png` | UMAP before vs. after decontamination. |

All four visual outputs are reproduced in the dashboard's Cell Filtering tab.

</details>

## Doublet detection

<details markdown="1">
<summary><code>doublet_filtering/&lt;analytical_run&gt;/</code></summary>

Present when `perform_doublet_detection = true` (the default). Runs on the **cell-called**
matrix.

| Path | Description |
| ---- | ----------- |
| `scrublet/<id>_scrublet_results.tsv` | Per-barcode Scrublet doublet scores and calls. |
| `scrublet/<id>_doublet_score_histogram.png` | Doublet score distribution. A clearly bimodal distribution means the automatic threshold is trustworthy. |
| `scdblfinder/<id>_scDblFinder_results.tsv` | Per-barcode scDblFinder scores and calls. |
| `scdblfinder/<id>_scDblFinder_summary.tsv` | Summary counts. |
| `combined/<id>_combined_doublets_w_combined_assignments.tsv` | The consensus call per barcode, combining both methods under `doublet_consensus_method`. |
| `combined/<id>_combined_doublets_summary.tsv` | Consensus summary counts. |
| `<datatype>/<id>_<datatype>_doublet_filtered/` | Count matrix with doublets **removed**, as an `mtx`/`barcodes`/`features` triplet — only when `perform_doublet_filtering = true`. This is what ambient RNA removal then denoises and what `anndata/` publishes. |
| `<datatype>/<id>_<datatype>_doublet_filter_summary.txt` | What was removed. |
| `<datatype>/<id>_<datatype>_doublet_filtering_summary.png` | Visual summary of the filtering. |
| `10x_export/<id>_10x/` | Cell-called matrix exported in 10x format, used as input to the doublet callers. |

> [!NOTE]
> By default doublets are **annotated, not removed** — set `perform_doublet_filtering = true`
> to remove them. Without it, no `*_doublet_filtered/` matrix is written, `obs["doublet_status"]`
> in `anndata/` flags them, and the calls in `combined/` are yours to apply.
>
> Detection is independent of `ambient_rna_remover`: the calls reach `anndata/` whether or not
> CellSweep runs.
>
> Doublet detection can legitimately fail on sparse libraries: scDblFinder's k-nearest-neighbour
> step needs more cells than a poor library retains, and Scrublet needs a bimodal score
> distribution. Both steps are configured to skip rather than abort, so a missing subdirectory
> here means the sample was left unannotated — check the run log for
> `[WARNING] No doublet calls for ...`. See
> [When doublet detection cannot run](CONFIGURATION_PARAMETERS.md#when-doublet-detection-cannot-run).

</details>

## Metacells and filtering

Present when `perform_metacells = true`. Each analytical run's cell-called AnnData object is
grouped into Metacell2 metacells, and `filtering_report.html` turns your choices into a
selection. See
[Metacells and the filtering report](CONFIGURATION_PARAMETERS.md#metacells-and-the-filtering-report).

<details markdown="1">
<summary><code>metacells/</code></summary>

| Path | Description |
| ---- | ----------- |
| `gene_table/gene_table.tsv` | One row per GTF gene: `gene_id`, `gene_name`, `chrom`, `biotype`, `is_mito` (on `mt_contig`), `is_rrna` (biotype matches `grep_rrna`) and `pfam` (domains from `gene_annotation`). |
| `gene_table/gene_table_stats.json` | How many genes are mitochondrial or rRNA, and how many annotation rows resolved to a gene. |
| `<id>/<id>_mc2_cells.h5ad` | Every called cell. `X` holds the raw counts, `obs["metacell"]` its group, plus per-cell `total_umis`, `pct_mito`, `pct_rrna`, the doublet call, CellSweep's `alpha_hat` and `layers["cellsweep"]` when available. `uns["mc_fingerprint"]` identifies the assignment. |
| `<id>/<id>_mc2_metacells.h5ad` | One row per group, with summed raw UMIs over every gene and the group's QC in `obs`: cells, UMIs, mito %, doublet %, mean `alpha_hat`, UMAP coordinates and top markers. |
| `<id>/<id>_mc_summary.json` | The report's input for this sample. For a sample with too few cells it holds `status: skipped` and the reason, and the two h5ad files are absent. |

Besides the metacells (named as Metacell2 names them, e.g. `M12.34`), `obs["metacell"]` uses two
pseudo-groups: `__outliers__` for cells Metacell2 left out of every metacell, and `__excluded__`
for cells it never grouped (outside `mc2_min_cell_umis`/`mc2_max_cell_umis`).

</details>

<details markdown="1">
<summary><code>filtering_report.html</code></summary>

A single HTML file covering every analytical run, selected with a dropdown. Plotly loads from
its CDN, so the first view needs an internet connection.

| Tab | Contents |
| --- | -------- |
| Cell-level filtering | Histograms of mitochondrial %, doublet % and mean ambient fraction per metacell, with a draggable threshold line, Metacell2's UMAP (click or lasso to blacklist) beside the current blacklist, and a sortable metacell table. |
| Gene-level filtering | Ranked total UMIs per gene and their histogram, with a draggable minimum on both, PFAM domains by tick-box, preset or regular expression, mitochondrial genes, a gene blacklist beside the genes those settings exclude, and a sortable table of all genes with their status. |
| Export selection | *Apply settings to all samples* and *Reset all samples to defaults*, then three steps: a table of every sample with the status of its two levels (*Applied* or *Pending*) and what it keeps, the download of `metacell_selection.json`, and the rerun that applies it. Below them, *Load selection* continues from an earlier export. |

Each sample has its own cell- and gene-level settings. A level is *pending* until you open its
tab with the sample selected, which applies it, or copy another sample's settings with *Apply
settings to all samples* (which leaves each sample's blacklists alone). Only samples with both
levels applied go into the selection.

The report remembers your choices in the browser between visits. *Reset all samples to
defaults* clears them.

</details>

<details markdown="1">
<summary><code>metacell_filtering/&lt;analytical_run&gt;/&lt;id&gt;_final/</code></summary>

Written by the `-resume --metacell_selection metacell_selection.json` rerun.

| Path | Description |
| ---- | ----------- |
| `matrix.mtx`, `barcodes.tsv`, `features.tsv` | The kept cells × kept genes, raw UMIs, in the 10x orientation. |
| `<id>_final.h5ad` | The same cells, with `obs["metacell"]`, the per-cell QC and annotations, and `layers["cellsweep"]` when present. Under `metacell_gene_mode = "flag"` every gene is kept, with `var["pass_filter"]` and `var["exclusion_reasons"]`. |
| `<id>_final_metacells.h5ad` | The kept metacells × genes, summed raw UMIs. |
| `<id>_filter_summary.json` | What was kept and from what (cells, metacells, genes, UMIs), with the rules the selection recorded. |

</details>

## Optional analyses

<details markdown="1">
<summary><code>saturation/&lt;analytical_run&gt;/</code></summary>

Produced when saturation analysis is enabled. Reads are progressively downsampled and the
library re-quantified at each depth.

| File | Description |
| ---- | ----------- |
| `saturation_output.tsv` | The saturation table: reads sampled against genes and UMIs detected. |
| `<id>_saturation.png` | Saturation curve, with `saturation_target` marked. |
| `<id>_saturation_residuals.png` | Residuals of the curve fit. |
| `<id>_saturation.log` | Fit log, including the estimated reads needed to reach the target. |

Both images appear in the dashboard's Saturation tab.

</details>

<details markdown="1">
<summary><code>rRNA_mtDNA/</code></summary>

Produced for every STARsolo run when `star_generateBAM = true`, by one pass over its BAM
(`CALC_READ_METRICS`).

| File | Description |
| ---- | ----------- |
| `<id>_mt_rrna_metrics.txt` | mtDNA and rRNA read counts and percentages, for the whole library and for the called cells, plus the split of mtDNA reads by strand. A CSV file whose metric names are quoted, because they contain commas. |
| `<id>_antisense_metrics.txt` | Uniquely mapped reads sense and antisense to genes, from STARsolo's `CellReads.stats`, and the antisense percentage of the two together, also split into exonic and intronic. For the whole library and for the called cells. Not produced when `star_soloStrand = "Unstranded"`, or without `CellReads.stats`. |
| `<id>_barcode_reads.tsv.gz` | Mapped, mtDNA and rRNA reads (primary alignment) of every barcode, from the same pass. A gzipped TSV with the columns `CB`, `MappedReads`, `MTReads` and `rRNAReads`. The per-cell metrics are computed from it. |

**Reads and alignments.** STAR writes one BAM record per alignment, so a read mapped to six
loci has six records. Every metric named "reads" counts each read once, at its primary
alignment. Only the rows named "alignments (all alignments)" count every record. Those rows are
diagnostics: they show where multimapper alignments land, such as rRNA repeats and nuclear
copies of mtDNA (NUMTs).

**mtDNA reads** are reads whose primary alignment is on a contig listed in `mt_contig`, i.e.
on the mitochondrial genome. The name describes where the read maps, not the molecule it came
from. In RNA-seq nearly all of these reads come from mitochondrially encoded transcripts. This
is the read-level counterpart of `percent.mt` / `pct_counts_mt`.

**rRNA reads** are reads overlapping any feature whose biotype contains "rRNA" (rRNA, Mt_rRNA,
rRNA_pseudogene, …), plus every read on a contig of `ref_gtf_addfeature`. Mt_rRNA (MT-RNR1/2)
is polyadenylated and is often most of the rRNA signal in poly(A)-captured libraries. A read on
it is both an mtDNA and an rRNA read, so the two percentages overlap and must not be summed.

The metrics, each as a count and a percentage:

| Denominator | Meaning |
| ----------- | ------- |
| `of mapped reads, primary alignment` | The headline value. Nuclear rRNA comes from multi-copy arrays and multimaps, so this is the fraction of reads that are rRNA. |
| `of uniquely mapped reads` | Uniquely mapped reads only (NH = 1). A lower bound for rRNA. For mtDNA, the gap to the mapped-reads value shows how much of the signal is ambiguous with NUMTs. |
| `of multimapped reads, primary alignment` | Multimappers only, each at its primary alignment. |
| `of multimapped read alignments, all alignments` | Every alignment of every multimapper (diagnostic). |
| `…, called cells` | The mapped and uniquely mapped values again, restricted to the barcodes of the filtered matrix. This is the cell set from cell calling, before doublet removal and CellSweep. `N/A` when no filtered matrix exists. |

Rows without a scope cover the whole BAM: empty droplets, ambient RNA and barcodes outside the
whitelist included. The two scopes can differ a lot. In single-nucleus data mitochondrial
transcripts concentrate in the ambient fraction, so a library value well above the called-cell
value points to ambient contamination.

**mtDNA reads by origin.** The mtDNA reads are split into those sense to a mitochondrial gene,
antisense to one, and outside every mitochondrial gene, as percentages of the mtDNA reads. The
genes are the exon rows of the main annotation on `mt_contig`. Transcripts land sense to the
genes. Fragments of mtDNA itself land on both strands about equally and also cover
untranscribed stretches. The mitochondrial genome is transcribed from both strands almost end
to end (L-strand transcripts, 7S RNA from the control region), so this is not a measurement of
DNA contamination but an upper bound on it. The useful comparison is the antisense share here
against the genome-wide antisense share in `<id>_antisense_metrics.txt` (the exonic one is the
closer match, since the mitochondrial split counts exons only): a markedly higher
value on the mitochondrial genome suggests DNA-derived reads. In poly(dT)-primed protocols DNA
is captured only through internal priming, so values near the genome-wide background are
expected. When `star_soloStrand = "Unstranded"`, only the share outside the genes is reported.

**Reads vs UMIs.** These are read percentages. They include PCR duplicates, so they differ from
UMI-based fractions such as the metacell report's mitochondrial UMI percentage, or a
`percent.mt` threshold set downstream.

High rRNA is usually a library preparation issue. High mtDNA in called cells usually indicates
stressed or damaged cells. In single-nucleus data it indicates cytoplasmic or ambient
contamination instead, since nuclei hold no mitochondria. Both are surfaced in the dashboard's
Mapping tab and in the per-cell metrics.

**Antisense reads.** STARsolo counts the strand set by `star_soloStrand`, and flags every
uniquely mapped read as exonic, intronic, exonicAS or intronicAS against the `GeneFull_Ex50pAS`
gene model (`CellReads.stats`, written with `soloCellReadStats Standard`, with or without
`perform_velocity`). The antisense percentage is `(exonicAS + intronicAS)` over all four, so it
uses the same alignments, gene model and strand as the count matrix and needs no pass over the
BAM. The library value sums every barcode, `CBnotInPasslist` included; the called-cell value
sums the filtered matrix's barcodes. `GeneFull_Ex50pAS` leaves fully exonic antisense reads out
of the matrix; intronic antisense reads can still be counted in it. A high antisense share
points to a wrong `star_soloStrand`, to internal priming on the opposite strand, or to
DNA-derived reads. The library value appears in the Quality Metrics card of the dashboard's
Mapping tab, and the metacell report shows the called-cell value.

</details>

<details markdown="1">
<summary><code>gene_ext/</code></summary>

Produced when `perform_geneext = true` or `run_method = "geneext_only"`.

| File | Description |
| ---- | ----------- |
| `geneext.gtf` | The extended gene annotation. |
| `geneext.gtf.Report.html` | GeneExt's own standalone report, including the per-gene extension table. |
| `geneext.gtf.GeneExt.log` | GeneExt's run log. |

GeneExt is run **once** on the merged alignments of all samples in the run, not per sample —
pooling gives enough 3′ coverage to support extensions that a single sample would not.

Because they are pooled, every sample is first capped at `geneext_subsample_nreads` reads
(50M by default), so a deeply sequenced sample cannot dominate the MACS2 peaks that decide
each gene's extension. Samples already under the cap are merged untouched, and a run with a
single sample is not subsampled at all. This is not the
same as `geneext_subsamplebam`, which is GeneExt's own cap on the *merged* BAM and downsamples
every sample by the same fraction, leaving the imbalance in place.

The headline statistics from `geneext.gtf.Report.html` — how many genes were extended and by
how much, how the MACS2 peaks were filtered, and which `--maxdist` was used — are summarised
in the dashboard's **Gene Extension** tab. The per-gene extension table is not duplicated
there; open `geneext.gtf.Report.html` for it.

GeneExt always reads STARsolo alignments of the full data. With `mapping_software = "alevin"`
or `"alevin_subsampled_starsolo"`, or with `geneext_downstream = "geneext_only"`, that is a
STARsolo pass mapped for GeneExt alone, stripped to the BAM.

With `perform_geneext = true`, the pipeline then rebuilds the index and re-maps against
`geneext.gtf` with every mapper of `mapping_software`: `_geneext_starsolo` (or
`_geneext_subsampled_starsolo`, on the subsample) and `_geneext_alevinfry` analytical runs.
By default (`geneext_downstream = "geneext_only"`) these are the only runs reported and
analysed downstream. With `geneext_downstream = "both"` the standard runs are kept as well,
and comparing `<sample>_starsolo` against `<sample>_geneext_starsolo` in the dashboard shows
exactly what the extension changed.

With `run_method = "geneext_only"` only that STARsolo pass is mapped and the pipeline stops here, which is the mode to use when you
want to generate a reference improvement once and reuse it across many runs — pass it back as
`ref_gtf` on subsequent runs. In that mode `mapping_STARsolo/` carries the logs and the BAM
but no count matrices, no cell calls and no saturation curves: the steps that produce them
feed only the reporting that this mode never reaches. Set `geneext_bam_only = false` to map
as a `standard` run would and still stop after GeneExt — see
[The `geneext_only` run method](CONFIGURATION_PARAMETERS.md#the-geneext_only-run-method).

</details>

<details markdown="1">
<summary><code>kraken/</code></summary>

Produced when `perform_kraken = true`. Requires `star_generateBAM = true`, since the unmapped
reads are extracted from the BAM.

| File | Description |
| ---- | ----------- |
| `<analytical_run>.k2report` | Kraken 2 report: reads assigned per taxon, with the standard rank hierarchy. |
| `<analytical_run>*.sankey.html` | Pavian Sankey diagram of the classification. |
| `kraken_db/kraken_db_path.txt` | The resolved path of the database used. |

The Sankey data are embedded into the dashboard's Taxonomic Classification tab, so you do not
need to open these files separately.

</details>

<details markdown="1">
<summary><code>containers/</code></summary>

`containers/demuxafy/Demuxafy.sif` — the Demuxafy singularity image (~7.5 GB), downloaded
automatically and MD5-verified when `demuxafy_sif` is not set. It is published with
`overwrite: false`, so it is downloaded once and reused. Set `demuxafy_sif` to this path on
subsequent runs to skip the download entirely.

</details>

## Summary reports

<details markdown="1">
<summary><code>dashboard.html</code></summary>

**The place to start.** A single self-contained HTML file — no server, no installation, no
network access needed beyond loading the plotting library. Open it in any browser, or send it
to a collaborator.

Seven tabs, plus one that appears only when the run produced an extended annotation:

| Tab | Content |
| --- | ------- |
| **Mapping** | Read counts, Q30 rates, uniquely/multi/unmapped percentages, intronic read percentage, rRNA and mtDNA read percentages (library and called cells, and mtDNA reads by strand; see `rRNA_mtDNA/`), cells, saturation, median UMIs and genes. Includes a cross-sample overview table. |
| **Cell Calling** | Barcode-rank knee curves and second-derivative curves per sample, with the chosen cutoff marked, alongside `expected_cells` and the cutoffs each method would have chosen. |
| **Saturation** | Saturation and residual curves. |
| **Per-Cell Visualizations** | Per-cell scatter over intronic %, mtDNA %, rRNA %, mapped reads and cell/non-cell status. |
| **Taxonomic Classification** | The Kraken 2 Sankey diagram. |
| **Cell Filtering** | CellSweep ambient-contribution histogram, top ambient genes and UMAP comparison. |
| **Gene Extension** | Only with `perform_geneext = true`. How many genes GeneExt extended and by how much, the 3′-extension length distribution, the MACS2 peak-filtering flow and coverage threshold, the `--maxdist` that was applied and any annotation fixes GeneExt made. |
| **Configurations** | The run configuration merged with each sample's samplesheet row. |

All analytical runs of a sample appear together, so STARsolo, alevin-fry and GeneExt results
are compared side by side. The Gene Extension tab is the exception: GeneExt runs once per
pipeline run on the merged alignments of every sample, so that tab describes the extended
annotation as a whole and does not follow the sample selector.

> [!NOTE]
> The dashboard embeds curves as data rather than images, so its size scales with barcode
> count. Curves are therefore thinned under fixed budgets (3,000 points for knee and
> second-derivative curves, 50,000 non-cell barcodes), uniformly in log10(rank) so the
> rendering is visually indistinguishable from the complete curve on log axes. **Every reported
> number — cell counts, medians, thresholds, saturation — is computed upstream on complete data
> and is unaffected.**

The dashboard can also be regenerated from a finished results directory without re-running the
pipeline:

```bash
bin/generate_dashboard.py --result-dir /path/to/output_directory --output dashboard.html
```

</details>

<details markdown="1">
<summary><code>summary_results/</code></summary>

| Path | Description |
| ---- | ----------- |
| `mapping_stats.tsv` | One row per analytical run with the headline mapping and quantification statistics. The read-region columns (exonic, intronic, intergenic, antisense, mtDNA) are percentages of uniquely mapped reads. For STARsolo they are summed over the called cells from `CellReads.stats`; for sci-rocket they come from its per-sample totals. |
| `multiqc_report.html` | MultiQC report aggregating FastQC, STAR (`Log.final.out`), Kraken2 and Salmon (alevin-fry mapping) outputs. |
| `multiqc_data/` | The underlying MultiQC data tables. |
| `R_images/UMI_dist_*.png` | UMI distribution plots (runs that mapped with STARsolo). |
| `R_images/cells_genes_*.png` | Cell and gene count plots (runs that mapped with STARsolo). |
| `per-cell_metrics/<id>_metrics.csv` | Per-cell metrics table, one row per barcode. The columns are described below. |
| `per-cell_metrics/<id>_metrics.json` | The same metrics as embedded in the dashboard. |
| `per-cell_metrics/*.png` | Per-cell metric plots. |

The per-cell read columns come from `rRNA_mtDNA/<id>_barcode_reads.tsv.gz`, so they use the same
definitions as `rRNA_mtDNA/`, and the called cells add up to its called-cell rows. Percentages
run from 0 to 100.

| Column | Meaning |
| ------ | ------- |
| `MappedReads` | Reads with this barcode, each counted once at its primary alignment. |
| `MTPercent` | mtDNA reads, as a % of `MappedReads`. |
| `rRNAPercent` | rRNA reads, as a % of `MappedReads`. |
| `IntronicPercent` | Reads STARsolo flags intronic, as a % of the barcode's uniquely mapped reads (`CellReads.stats`, `intronic / genomeU`). This is the same classification as the dashboard's library-level intronic percentage. |
| `UnsplicedPercent` | Only with `perform_velocity`. Unspliced UMIs as a % of all Velocyto UMIs (spliced + unspliced + ambiguous). Velocyto classifies each UMI once against the transcript models of the `Gene` feature, so this is a UMI-level cross-check of `IntronicPercent`, not expected to equal it. |
| `IsCell` | 1 for a barcode of the filtered matrix. |

</details>

## External pipelines

Produced only when the corresponding vendor pipeline is enabled. These run alongside the BCA
pipeline on the same input, as a validation reference — their outputs do **not** feed into any
BCA pipeline channel.

| Directory | Pipeline | Notes |
| --------- | -------- | ----- |
| `CellRanger_pipeline/` | 10x Genomics Cell Ranger | Containerised; no manual installation needed. Standard `outs/` tree per sample. |
| `ParseBio_pipeline/<sample>/` | Parse Biosciences split-pipe | Requires a local installation. |
| `BDrhapsody_pipeline/` | BD Rhapsody Sequence Analysis Pipeline | Requires a local installation. |

See [Installing external pipelines](INSTALLATION_EXTERNAL_PIPELINES.md).

## Which matrix should I use?

The pipeline deliberately produces several matrices. For a standard analysis:

1. **Start from** `mapping_STARsolo/<sample>_starsolo/<sample>_starsolo_Solo.out/GeneFull_Ex50pAS/filtered_secondderiv/filtered/`
   — the cell-called, exonic+intronic matrix.
   For alevin-fry, the equivalent is
   `mapping_alevin/<sample>_alevinfry/<sample>_alevinfry_counts/alevin/filtered_secondderiv/filtered/`.
2. **If you enabled ambient RNA removal**, prefer `cellsweep/<sample>_starsolo/<sample>_starsolo_cs_filtered.h5ad`,
   which is decontaminated and cell-called. To keep the raw counts alongside the denoised ones,
   take `anndata/<sample>_starsolo/raw/` instead: it carries both, plus every annotation
   (see [AnnData conversion](#anndata-conversion)).
3. **Doublet calls are already in** `anndata/.../obs["doublet_status"]`. The per-barcode tables
   behind them are in `doublet_filtering/<sample>_starsolo/combined/`, and with
   `perform_doublet_filtering = true` the doublets are gone from every published matrix.
4. **If you ran GeneExt**, take the `_geneext_starsolo` analytical run (or `_geneext_alevinfry`
   for alevin-fry) instead. By default these are the only runs; with
   `geneext_downstream = "both"` the dashboard's Mapping tab shows them next to the standard
   runs, so you can check that the extension improved gene detection.
5. **For RNA velocity**, start from `anndata/<sample>_starsolo/velocity/` instead — it carries
   the same cells with the splicing breakdown as layers. Requires `perform_velocity = true`.
6. **If you filtered metacells** with `filtering_report.html`, use
   `metacell_filtering/<sample>_starsolo/<sample>_starsolo_final/`. It holds the cells and genes
   the selection kept, with each cell's metacell attached.

Use the exonic-only `Gene/` matrix rather than `GeneFull_Ex50pAS/` only if you specifically want
to exclude intronic reads. For single-nucleus data, `GeneFull_Ex50pAS` is usually the
right choice.

The three matrices in short:

| Want | Use |
| ---- | --- |
| Exonic counts only | `Gene/` |
| Exonic + intronic counts (the default for most analyses) | `GeneFull_Ex50pAS/` |
| Intronic counts on their own, or spliced/unspliced for velocity | `Velocyto/`, or the velocity `.h5ad` |

Do **not** try to obtain intronic counts by subtracting `Gene/` from `GeneFull_Ex50pAS/`. They
are two independent counting passes, not a superset and a subset, so the difference goes
negative for some gene/cell pairs and mixes true intronic reads with reads that were merely
unassignable under exon-only rules. Set `perform_velocity = true` instead.

> [!TIP]
> Every matrix a sample produces describes the same set of cells, because cell calling is done
> once on `GeneFull_Ex50pAS` and the other matrices are subset to that result. You can move
> between them without re-intersecting barcodes.
