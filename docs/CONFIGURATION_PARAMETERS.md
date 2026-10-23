# Custom Configuration Variables


Within each custom configuration file the following variables can be defined:

## Table of Contents

1. [Base Variables](#base-variables) (and [How paths are resolved](#how-paths-are-resolved))
2. [STAR Variables](#star-variables) (and [STARsolo multi-gene reads](#starsolo-multi-gene-reads))
3. [Alevin-fry Variables](#alevin-fry-variables)
4. [Mapping-rate check](#mapping-rate-check)
5. [FeatureCounts Variables](#featurecounts-variables)
6. [Cell calling](#cell-calling)
7. [Gene Extension Variables](#gene-extension-variables)
8. [Taxonomic-classification Variables](#taxonomic-classification-variables)
9. [Filtering Variables](#filtering-variables)
10. [Metacells and the filtering report](#metacells-and-the-filtering-report)
11. [External Pipeline Variables](#external-pipeline-variables)
12. [Execution Reports](#execution-reports)

---


## Base Variables

| Variable               | Required/Optional | Description |
|------------------------|-------------------|-------------|
| `input`                | __Required__      | Path to the samplesheet. |
| `outdir`               | __Required__      | Path to the results/output directory; created if missing. |
| `bc_whitelist`         | __Required__      | Path or link to the barcode whitelist file(s), multiple ones separated by whitespace. If links are given, they are automatically downloaded (and unzipped if applicable) for any protocol. Not used by `"marsseq_v1"`/`"marsseq_v2"`, which run without a whitelist.|
| `protocol`              | __Required__     | Specifies the sequencing technology used (must be one of the following: `"oak_v1"`, `"10xv1"`, `"10xv2"`, `"10xv3"`, `"10xv4"`, `"parse_biosciences_WT_mini"` or `"parse_biosciences_WT"`,     `"bd_rhapsody_v1"`, `"bd_rhapsody_enhancedbeads"`, `"sciRNAseq3"` , `"ultima_genomics"`, `"marsseq_v1"`, `"marsseq_v2"` or `"seqspec"`). |
| `ref_fasta`            | __Required__      | Path to the genome FASTA file used for mapping reads. |
| `ref_gtf`              | __Required__      | Path to the GTF/GFF file formatted for STARsolo. |
| `ref_gtf_alt`          | Optional          | **Deprecated, no effect.** The Cell Ranger and split-pipe references are built from `ref_gtf`, like the STARsolo one. Setting it logs a warning. |
| `run_method`             | Optional          | Method of running the pre-processing pipeline, demonstrated in the [pipeline diagram](img/Preprocs_Pipeline.png), one of `"standard"`, `"geneext_only"`, `"external_pipeline_only"` or `"post_mapping"`. Default is set to `"standard"`. `"post_mapping"` picks a finished run up at the point mapping ended and redoes everything after it. |
| `previous_outdir`      | Optional          | Only used with `run_method = "post_mapping"`: the results directory the published mapping results are read back from. Defaults to `outdir`, i.e. re-running in place. |
| `mapping_software`     | Optional          | Software used to map reads (must be one of the following: `"starsolo"`, `"alevin"`, `"both"/"alevin_starsolo"` or `"alevin_subsampled_starsolo"`). Default set to `"starsolo"`. |
| `cellfilter_method`    | Optional          | How cells are called: `"star_solocellfilter"` (each mapper's own cell filter: STARsolo's `--soloCellFilter`, alevin-fry's knee), `"second_derivative"` or `"manual_cutoff"`. Default set to `"second_derivative"`. The latter two apply to both mappers: they filter the count matrix on a UMI cutoff and recompute every cell-dependent statistic in the report on the retained cells, differing only in where that cutoff comes from. See [Cell calling](#cell-calling). |
| `perform_demultiplexing` | Optional        | Boolean flag to enable or disable demultiplexing of the FASTQ files, where applicable. Default is `true`. |
| `perform_trimming`     | Optional          | Boolean flag to enable or disable trimming of adapters and low-quality bases with fastp, where applicable (sci-RNA-seq3 and seqspec). Default is `true`. |
| `oak_demux_max_mismatches` | Optional      | (OAK only) Mismatches allowed between an index read and the aliquot's i5/i7 barcode. Default is `1` (set in `conf/seqtech_parameters.config`). |
| `seqspec_file`         | Optional          | Path to the seqspec file. |
| `subsample_nreads`     | Optional          | The size (number of reads) of the subset used to map to STARsolo, in case the parameter `mapping_software = alevin_subsampled_starsolo`. Default set to `100000000` reads. |
| `alevin_usa_counts`    | Optional          | Which of alevin-fry's USA blocks are summed into each gene's count before ambient-RNA removal and doublet detection: `"SUA"` (spliced + unspliced + ambiguous), `"SA"`, `"S"`, `"UA"` or `"U"`. Default set to `"SUA"`. See [USA counts](#usa-counts). |
| `alevin_resolution`    | Optional          | alevin-fry's UMI resolution: `"cr-like"` discards UMIs that stay ambiguous between genes, as Cell Ranger and STARsolo with `star_soloMultiMappers = "Unique"` do; `"cr-like-em"` splits them between the genes by expectation-maximisation, giving fractional counts, the counterpart of STARsolo's default `"EM"`. Default set to `"cr-like"`. See [alevin-fry settings](#alevin-fry-settings). |
| `perform_velocity`     | Optional          | Boolean flag to additionally produce intronic counts for RNA velocity: STARsolo's `Velocyto` feature (spliced/unspliced/ambiguous) and alevin-fry's unspliced block, plus a velocity-ready `.h5ad` per mapper. Default is `false`. See [RNA velocity](#rna-velocity). |
| `validate_params`      | Optional          | Check every parameter and the samplesheet against `nextflow_schema.json` and `assets/schema_input.json` before anything runs: types, allowed values, and that given input files and directories exist. Default is `true`; set `false` to skip it. |


### How paths are resolved

- **Parameters** (`input`, `outdir`, `ref_fasta`, `ref_gtf`, `bc_whitelist`, `star_index`, ...): a relative path is read from the directory `nextflow run` is launched from. Use absolute paths, or build them in the config from `${projectDir}` (the pipeline's folder) or `${launchDir}`, so a config works from any launch directory.
- **FASTQs in the samplesheet**: a relative path is read from the samplesheet's own folder, so a samplesheet can sit next to its data. Samplesheets written for older versions, with paths relative to the launch directory, still work: such a path is used when it is not found next to the samplesheet, with a warning. Globs (e.g. `name_I*` in `fastq_indices`) are expanded.
- Absolute paths and cloud URIs (`s3://`, `gs://`, `az://`) are used as they are; in the samplesheet, `~/` is expanded to the home directory. Every input file is staged into the tasks, so it does not have to be visible at the same path on the compute nodes. The exceptions are the external installations (`rhapsody_installation`, `splitpipe_demultiplex_script`, `splitpipe_conda_env`), the Kraken database and `demuxafy_sif`, which are used where they are.
- A given path that does not exist stops the run before anything is scheduled. A FASTQ named in the samplesheet but missing does too, except with `run_method = "post_mapping"`, which reads only the sample names.



## STAR Variables

| Variable               | Required/Optional | Description |
|------------------------|-------------------|-------------|
| `star_index`           | Optional          | Path to the pre-generated STAR index. By default the STAR index is created within the pipeline.|
| `star_genomeSAindexNbases` | Optional         | Lenght of the SA pre-indexing string in STAR. See [protocol-specific defaults](../conf/seqtech_parameters.config) set in the seqtech_paramaters.config file. |
| `star_genomeSAsparseD`    | Optional       | Suffix array sparsity in STAR.  See [protocol-specific defaults](../conf/seqtech_parameters.config) set in the seqtech_paramaters.config file. |
| `star_solocellfilter`  | Optional          | Cell filtering type and parameters used by STAR. See [protocol-specific defaults](../conf/seqtech_parameters.config) set in the seqtech_paramaters.config file.|
| `star_soloFeatures`    | Optional         | Genomic features for which the UMI counts per Cell Barcode are collected. See [protocol-specific defaults](../conf/seqtech_parameters.config) set in the seqtech_paramaters.config file. `Velocyto` is appended automatically when `perform_velocity = true`, so it does not need to be listed here; see [RNA velocity](#rna-velocity). |
| `star_outFilterScoreMin` | Optional       | Alignment will be output only if its score is higher than or equal to this value normalized by read length. See [protocol-specific defaults](../conf/seqtech_parameters.config) set in the seqtech_paramaters.config file. |
| `star_outSAMunmapped`     | Optional      | Output of unmapped reads in the SAM format. See [protocol-specific defaults](../conf/seqtech_parameters.config) set in the seqtech_paramaters.config file. |
| `star_outSAMattributes` | Optional        | String of desired SAM attributes. See [protocol-specific defaults](../conf/seqtech_parameters.config) set in the seqtech_paramaters.config file. |
| `star_generateBAM`      | Optional        | Boolean to determine if a BAM file should be generated. Will automatically adjust the star_outSAMattributes. Default is set to `true`. |
| `star_soloTypestring` | Optional          | String of defining soloType and barcode structure. See [protocol-specific defaults](../conf/seqtech_parameters.config) set in the seqtech_paramaters.config file. |
| `star_soloCBmatchWLtype` | Optional       | Type of matching the cell barcodes to the barcode whiteList. See [protocol-specific defaults](../conf/seqtech_parameters.config) set in the seqtech_paramaters.config file. |
| `star_soloUMIfiltering` | Optional        | Type of UMI filtering. See [protocol-specific defaults](../conf/seqtech_parameters.config) set in the seqtech_paramaters.config file. |
| `star_soloMultiMappers` | Optional        | Counting method for reads mapping to multiple genes, and so which STARsolo matrix the pipeline analyses: `"Unique"` drops them (`raw/matrix.mtx`); `"EM"`, `"Uniform"`, `"PropUnique"` or `"Rescue"` adds them to the unique counts (`raw/UniqueAndMult-<type>.mtx`, fractional counts). With several types listed, the first one other than `Unique` is analysed. Every protocol defaults to `"EM"` (see [protocol-specific defaults](../conf/seqtech_parameters.config)). See [STARsolo multi-gene reads](#starsolo-multi-gene-reads). |
| `star_soloUMIdedup` | Optional         | Type of UMI deduplication algorithm. See [protocol-specific defaults](../conf/seqtech_parameters.config) set in the seqtech_paramaters.config file. |
| `star_clipAdapterType` | Optional         | Type of adapter clipping. See [protocol-specific defaults](../conf/seqtech_parameters.config) set in the seqtech_paramaters.config file. |
| `star_limitBAMsortRAM` | Optional         | STAR's own ceiling on the BAM sort buffer, in bytes. Default `0`, which makes the pipeline derive it from `task.memory` minus the resident genome index. Set a number only to pin it, e.g. to the figure a STAR error asks for. Note this is *not* the scheduler allocation — see [Resource Tuning](RESOURCE_TUNING.md#tool-internal-memory-ceilings). |
| `star_limitGenomeGenerateRAM` | Optional  | STAR's own ceiling for index construction, in bytes. Default `null`, which derives it from `task.memory`. A pinned value neither grows with a larger resource tier nor shrinks with a smaller one. |


### STARsolo multi-gene reads

A multi-gene read aligns equally well to two or more genes. It either maps uniquely to a
spot where genes overlap, or maps to several loci that belong to different genes
(paralogues, duplicated genes, multi-copy arrays). Cell Ranger, and STARsolo with
`"Unique"`, discard them. With a multi-gene type, STARsolo adds them to the unique counts:
`"EM"` distributes each such UMI among its genes by maximum likelihood over the other UMIs
of the same cell (Kaminow et al. 2021). STAR always writes the unique counts to
`raw/matrix.mtx`; `raw/UniqueAndMult-<type>.mtx` holds the unique plus the multi-gene
counts, and is the matrix the pipeline analyses whenever `star_soloMultiMappers` lists such
a type:

- **Cell calling.** With `second_derivative` or `manual_cutoff`, the cutoff is read off a
  curve derived from that matrix, `filtered_secondderiv/<id>_EM_UMIperCellSorted.txt`
  (STAR's own `UMIperCellSorted.txt` counts unique-gene UMIs only), and filters it; the
  filtered matrix is written as `filtered_secondderiv/filtered/matrix.mtx`. With
  `star_solocellfilter`, STAR's own cells are kept (STAR calls them on unique counts), and
  the multi-gene matrix is subset to them in `GeneFull_Ex50pAS/filtered_UniqueAndMult-EM/`;
  STAR's own `filtered/` keeps the unique matrix. The dashboard's median UMIs, median genes
  and genes detected are those of the subset matrix
  (`GeneFull_Ex50pAS/<id>_UniqueAndMult-EM_statistics.json`); the cells, the cutoff and the
  read-level statistics stay STAR's.
- **Downstream.** CellSweep, doublet detection, the `.h5ad` files and Metacell2 receive the
  multi-gene counts. They are fractional: `X` is a float matrix, and Metacell2 rounds them
  for its grouping only. A gene counts as *detected* in a cell from one UMI on, so the
  small fractions EM spreads over a gene's paralogues are not counted as genes.
- **Unchanged.** `Summary.csv` (STAR's own statistics), the read-level statistics from
  `CellReads.stats` (sequencing saturation, fraction of reads in cells) and the Velocyto
  matrices count unique-gene reads only.


## Alevin-fry Variables

alevin-fry is set up to count the same molecules as the default STARsolo run, following
the alevin-fry documentation where it gives a recommendation. Settings marked *fixed* are
not configurable. Geometries use salmon's custom syntax, `"<read>[<ranges>]"`: 1-based
inclusive ranges, concatenated in the order listed. Read `1` is the CB/UMI FASTQ and read
`2` the cDNA FASTQ, whichever of R1/R2 carries the cDNA.

| Variable | Default | Description |
|----------|---------|-------------|
| `alevin_bc_geometry` | [protocol-specific](../conf/seqtech_parameters.config) | Position(s) of the cell barcode, e.g. `"1[1-16]"` or `"1[1-9,14-22,27-35]"`. |
| `alevin_umi_geometry` | [protocol-specific](../conf/seqtech_parameters.config) | Position of the UMI, e.g. `"1[17-28]"`. |
| `alevin_read_geometry` | `"2[1-end]"` (every protocol) | Position of the biological sequence. |
| `alevin_resolution` | `"cr-like-em"` | `"cr-like"` is Cell Ranger's rule (STARsolo's `"Unique"`): a UMI ambiguous between genes is discarded. `"cr-like-em"` splits it by EM instead, the counterpart of STARsolo's default `"EM"`. |
| `alevin_usa_counts` | `"SUA"` | Which spliced/unspliced/ambiguous blocks are summed into each gene's count. See [USA counts](#usa-counts). |
| Orientation | *fixed*, follows `star_soloStrand` | `generate-permit-list -d`: `Forward` → `fw`, `Reverse` → `rc`, `Unstranded` → `both`. Both mappers count the same strand; `GeneFull_Ex50pAS` leaves antisense reads out too. |
| Permit list | *fixed*, follows `cellfilter_method` | `--unfiltered-pl`, `--min-reads 1` | Every whitelisted barcode with a read is kept, a raw matrix like STARsolo's `raw`, for the pipeline's own cell call, so the ambient-RNA step sees the same empty droplets (alevin-fry's default is 10 reads). Falls back to alevin-fry's knee (`-k`) with no whitelist (MARS-seq), with `cellfilter_method = "star_solocellfilter"`, or when a multi-part whitelist expands to more than 10 million barcodes. The task log says which permit list was used and why. |
| Multi-segment barcodes (Parse, BD Rhapsody, sci-RNA-seq3) | *fixed*, from `bc_whitelist` | Every combination of the per-segment whitelist files, joined in segment order. alevin-fry corrects the joined barcode as a whole: one edit over the whole barcode, where STARsolo allows one mismatch per segment. |


### USA counts

alevin-fry runs in USA mode: `quants_mat.mtx` holds each gene three times, in three equal
blocks — spliced (S, bare gene ID), unspliced (`<gene>-U`) and ambiguous (`<gene>-A`). Left
as is, every gene would reach ambient-RNA removal and doublet detection as three correlated
features, distorting their highly-variable-gene selection and PCA. The pipeline therefore
splits the blocks by position (as pyroe and simpleaf do, so gene IDs ending in `-A` or `-U`
are safe; any other layout stops the run) and sums the ones chosen by `alevin_usa_counts`
into one column per gene. The result is published under
`mapping_alevin/<sample>/gene_level_matrix/`.

| Value   | Counts                          | STARsolo counterpart  | When to use |
|---------|---------------------------------|-----------------------|-------------|
| `"SUA"` (default) | spliced + unspliced + ambiguous | `GeneFull_Ex50pAS` | **The only setting comparable with STARsolo**, whose analysed matrix is fixed to `GeneFull_Ex50pAS` (exons *and* introns). Introns are counted for single-cell and single-nucleus data alike, as Cell Ranger now does by default |
| `"SA"`  | spliced + ambiguous             | `Gene`                | The conventional exonic count (as in `fishpond::loadFry`), when wanted for its own sake; it no longer lines up with STARsolo. It, not `"S"`, is the counterpart of `Gene`, because A holds reads that fit an exon and its neighbouring intron's flank equally well, mostly exonic |
| `"S"`   | spliced only                    | --                    | Strictly spliced counts, with the same caveat as `"SA"` |
| `"UA"`  | unspliced + ambiguous           | --                    | Rarely: it makes unspliced counts the matrix every downstream step and statistic runs on. For velocity, set `perform_velocity = true` instead, which publishes the unspliced block *alongside* the main matrix; see [RNA velocity](#rna-velocity) |
| `"U"`   | unspliced only                  | `Velocyto/` unspliced | As `"UA"` |

**Consistency with cell calling.** With `"second_derivative"`, the per-cell UMI totals the
cutoff is derived from and the alevin statistics the report shows against it (cells, median
UMIs and genes per cell, total genes detected) count the same blocks, so a cutoff of *n* UMIs
means *n* UMIs of the analysed kind. With `"star_solocellfilter"` the report shows
alevin-fry's own `quant.json` and `cell_meta.tsv` numbers instead, which will not match the
analysed matrix under a non-default `alevin_usa_counts`.



## Mapping-rate check

A sample whose reads mostly fail to map (contamination, wrong species or reference, or a wrong protocol/barcode geometry) would otherwise spend hours in the mapper and then run every downstream step on almost nothing. The pipeline therefore cancels any sample that maps below `min_mapping_rate`.

| Variable                   | Required/Optional | Description |
|----------------------------|-------------------|-------------|
| `skip_mapping_rate_check`  | Optional          | Boolean flag to turn the check off, so that no sample is ever cancelled for a low mapping rate. Default is `false`. |
| `min_mapping_rate`         | Optional          | Percentage below which a sample is cancelled: *uniquely* mapped reads for STARsolo, mapped reads (`percent_mapped`) for salmon/alevin-fry. `0` also disables the check. Default set to `10`. |
| `mapping_rate_action`      | Optional          | What a failed check does: `"abort"` stops the whole pipeline, `"skip_sample"` drops only that sample (for every mapper) and lets the other samples finish. Default set to `"abort"`. |
| `mapping_rate_check_reads` | Optional          | STARsolo only: check the rate from `Log.progress.out` while STAR is still mapping, once this many reads are processed. `0` checks only when STAR finishes. Default set to `20000000`. |
| `mapping_rate_poll_secs`   | Optional          | How often a running mapper is checked, in seconds. STAR refreshes `Log.progress.out` once a minute. Default set to `60`. |

**When the rate is checked**

- **STARsolo** writes its running statistics to `Log.progress.out` once a minute. After `mapping_rate_check_reads` reads have been processed, the "Mapped unique" percentage is checked and STAR is stopped straight away if it is below the threshold. Samples with fewer reads than that are checked from `Log.final.out` when STAR finishes.
- **Alevin-fry** has no progress file that can be parsed. Its rate is read from salmon's `aux_info/meta_info.json` as soon as `salmon alevin --justAlign` finishes, *before* the alevin-fry permit-list, collate and quant steps.

**When both mappers run** (`both`, `alevin_starsolo` or `alevin_subsampled_starsolo`), the first mapper to reach a verdict decides for the sample. If it fails, the other mapper is stopped at its next poll (or never starts mapping). If it passes, the other mapper skips its own check. With `alevin_subsampled_starsolo`, this usually means the subsampled STARsolo run decides for the full-depth alevin-fry run. The GeneExt re-mapping (`_geneext_starsolo` and `_geneext_alevinfry`) reuses the verdict of the first pass.

**What you see.** The cancelled task exits with code `42` and prints a `MAPPING RATE CHECK FAILED` block in its log. At the end of the run, the pipeline prints a summary of every cancelled sample and writes all verdicts to `pipeline_info/mapping_rate_check.tsv`. Both messages recommend finding out what the unmapped reads are with Kraken2. The pipeline's own Kraken step reads unmapped reads from the STARsolo BAM, which a cancelled sample never produces, so rerun the sample with the check disabled:

```bash
--skip_mapping_rate_check true --perform_kraken true --star_generateBAM true
```

or run `kraken2` directly on a subsample of the cDNA FASTQ.

The mappers share their verdicts through files under `<workDir>/.mapping_rate_check/`, so the check needs a work directory that all tasks can see, as is the case on a shared-filesystem cluster. `run_method = "post_mapping"` does not map and is not checked.


## FeatureCounts Variables

| Variable               | Required/Optional | Description |
|------------------------|-------------------|-------------|
| `perform_featurecounts`  | Optional        | **Deprecated, no effect.** The mtDNA, rRNA and antisense read percentages are now computed for every STARsolo run with a BAM (`star_generateBAM = true`), in the same pass that gives the per-cell metrics their read counts. Setting it logs a warning. |
| `mt_contig`            | Optional          | Name(s) of the mitochondrial contig(s) in the reference annotation, used to calculate mtDNA content. Multiple contigs can be given separated by whitespace, in which case reads on any of them count as mitochondrial. The genes on these contigs are also the mitochondrial genes of the metacell report and of CellSweep's MT filter before clustering. Default set to `"chrM M MT"`. |
| `grep_rrna`            | Optional          | String used to grep ribosomal RNA (rRNA) reads from annotations. Default set to `"rRNA"`|



## Cell calling

| Variable               | Required/Optional | Description |
|------------------------|-------------------|-------------|
| `cellfilter_method`    | Optional          | How cells are separated from empty droplets: `"star_solocellfilter"`, `"second_derivative"` or `"manual_cutoff"` (see the methods below). Default set to `"second_derivative"`. |
| `star_solocellfilter`  | Optional          | STARsolo's `--soloCellFilter`; only decides the cells with `cellfilter_method = "star_solocellfilter"`. See [STAR Variables](#star-variables). |
| `expected_cells` (samplesheet column) | Optional | `"second_derivative"` only: the knee is searched for between 0.1x and 10x this value, which keeps the search out of a deep ambient tail. It bounds the search but does not pull the cell count towards it (300 real cells with `expected_cells = 3000` still call ~300). It is capped to the number of barcodes above 100 UMIs, and the report flags the sample when the cap applies. Left blank, the search runs unbounded and the report warns. Ignored under `"manual_cutoff"`. |
| `manual_cutoff` (samplesheet column) | __Required__ with `"manual_cutoff"` | The UMI threshold, one value per sample. Every sample needs one under `cellfilter_method = "manual_cutoff"`, or the run is rejected at start-up; under any other method the column is ignored with a warning. |

**Cell calling methods**

| `cellfilter_method`  | Where the cells come from |
|----------------------|---------------------------|
| `"star_solocellfilter"` | Each mapper's own cell call: STARsolo's `--soloCellFilter` and alevin-fry's knee-based permit list. |
| `"second_derivative"` (default) | The most negative second derivative of the log-log cumulative UMI curve, within the `expected_cells` window. The curve is first resampled onto a grid uniform in log10(rank), so the integer quantisation of the tail does not bury the knee. Adapts to knees too shallow for the mappers' built-in filters. |
| `"manual_cutoff"`    | The sample's `manual_cutoff` column; the search is skipped. Use it when a curve has no knee the search finds, or when the cutoff was chosen by eye from a previous run's report. Everything else is as with `"second_derivative"`. |


## Gene Extension Variables

| Variable               | Required/Optional | Description |
|------------------------|-------------------|-------------|
| `perform_geneext`      | Optional          | Boolean flag to enable or disable the gene extension step after mapping. GeneExt is built from full-data STARsolo alignments, then every mapper of `mapping_software` is re-mapped against the extended annotation (`_geneext_starsolo` or `_geneext_subsampled_starsolo`, and `_geneext_alevinfry`). Works with every `mapping_software`; see [Which runs GeneExt adds](#which-runs-geneext-adds). Default is `false`. |
| `geneext_downstream`   | Optional          | Only consulted with `perform_geneext`. `"geneext_only"` maps, reports and analyses downstream (doublets, CellSweep, AnnData, metacells) only the runs against the extended annotation; `"both"` also the standard-annotation runs, for a side-by-side comparison at roughly twice the cost. Default is `"geneext_only"`. |
| `geneext_gtf`          | Optional          | Only used with `run_method = "post_mapping"`: the extended annotation to pick up instead of running GeneExt again. Defaults to `<previous_outdir>/gene_ext/geneext.gtf`. |
| `geneext_subsample_nreads` | Optional      | Per-sample cap on the reads reaching GeneExt, applied to each sample's BAM *before* they are merged. Default is `50000000`; `0` disables it. See [Balancing the samples](#balancing-the-samples-that-reach-geneext). |
| `geneext_subsamplebam` | Optional          | GeneExt's own `--subsamplebam`. Applied to the **merged** BAM, so it caps GeneExt's total input but downsamples every sample by the same fraction and therefore does **not** balance them against each other. Unset by default. |
| `geneext_bam_only`     | Optional          | Strips a STARsolo pass that is mapped only for GeneExt down to the BAM GeneExt reads: with `run_method = "geneext_only"`, and with `perform_geneext` wherever the standard `_starsolo` run is not mapped (see below). Default is `true`. See [The `geneext_only` run method](#the-geneext_only-run-method). |

### Which runs GeneExt adds

GeneExt always reads STARsolo alignments of the **full** data, never of a subsample. Where
the standard `<sample>_starsolo` run is mapped anyway (`mapping_software = "starsolo"` or
`"both"`, with `geneext_downstream = "both"`) its BAMs are used. Otherwise a STARsolo pass is
mapped for GeneExt alone, also named `<sample>_starsolo`. It is stripped to the BAM (see
`geneext_bam_only`), and it is neither shown in the dashboard nor carried downstream.

Every mapper of `mapping_software` is then re-mapped against the extended annotation, and
`geneext_downstream` decides whether the standard-annotation runs are mapped as well:

| `mapping_software` | `geneext_downstream = "geneext_only"` (default) | Added with `"both"` |
|--------------------|-------------------------------------------------|---------------------|
| `starsolo` | `_geneext_starsolo` | `_starsolo` |
| `alevin` | `_geneext_alevinfry` | `_alevinfry` |
| `both` / `alevin_starsolo` | `_geneext_starsolo`, `_geneext_alevinfry` | `_starsolo`, `_alevinfry` |
| `alevin_subsampled_starsolo` | `_geneext_subsampled_starsolo`, `_geneext_alevinfry` | `_subsampled_starsolo`, `_alevinfry` |

Every run in the table is reported in the dashboard and MultiQC and goes through doublet
detection, ambient-RNA removal, the AnnData export and, with `perform_metacells`, the
metacells. With the default, the dashboard therefore has no standard-annotation run to
compare each `_geneext_` run against; GeneExt's own statistics are still in its
**Gene Extension** tab. Set `geneext_downstream = "both"` when that comparison is wanted.

`_geneext_subsampled_starsolo` maps the same subsample of `subsample_nreads` reads as
`_subsampled_starsolo`, so `mc2_skip_subsampled` skips its metacells too.

### Balancing the samples that reach GeneExt

GeneExt runs **once** on the merged alignments of every sample in the run, because pooling
gives enough 3′ coverage to support extensions a single sample would not. The side effect is
that a deeply sequenced sample dominates the MACS2 peaks that decide each gene's extension,
so one tissue can end up setting the annotation used for all of them.

`geneext_subsample_nreads` caps every sample at the same read count before the merge, which
is what actually equalises their contribution. Samples already under the cap are passed
through untouched, and the seed is fixed, so a re-run subsamples identically. A run with a
single sample is never subsampled, since there is nothing to balance it against and the cap
would only throw away reads.

Note that this applies to **both** `perform_geneext = true` and
`run_method = "geneext_only"`, since the pooling is the same in either case. Set
`geneext_subsample_nreads = 0` to merge the samples exactly as sequenced.

### The `geneext_only` run method

`run_method = "geneext_only"` exists to produce one artefact — `gene_ext/geneext.gtf` — to be
reused as `ref_gtf` on later runs. Whatever `mapping_software` says, only the full-data
STARsolo pass GeneExt reads (`<sample>_starsolo`) is mapped. Nothing after mapping runs in this mode, so with
`geneext_bam_only` left at its default the pipeline stops paying for the outputs those steps
would have consumed:

- STARsolo counts only the cheap `Gene` feature (STAR has no `--soloFeatures None`), and its
  raw matrix is discarded. `Summary.csv` and `Barcodes.stats` are kept as QC on the alignment.
- `--soloCellFilter None` and `--soloCellReadStats None`, so no cell calling — neither
  STARsolo's own nor the pipeline's second-derivative call.
- The BAM is written without the `CR`/`UR`/`CB`/`UB` (and `GX`/`GN`/`sM`) tags and without
  unmapped reads, since GeneExt reads 3′ coverage rather than tags. It is correspondingly
  smaller and faster to sort, merge and read.
- The BAM inspection stack — saturation curve, read metrics, Kraken — is skipped.

Set `geneext_bam_only = false` to map exactly as a `standard` run would and still stop after
GeneExt, which is the mode to use when the count matrices and QC of the pre-extension mapping
are wanted alongside the annotation. The BAM is generated either way: `geneext_only` with
`geneext_bam_only = false` **and** `star_generateBAM = false` is rejected up front, since it
would leave GeneExt with nothing to read.


## 10x_saturate Variables

| Variable               | Required/Optional | Description |
|------------------------|-------------------|-------------|
| `perform_10x_saturate` | Optional          | Boolean flag to enable or disable the 10x_saturate step after mapping. Default is `true`. |
| `saturation_target`    | Optional          | The saturation target fraction used to predict the input reads needed. Default set to `0.7`. |



## Taxonomic-classification Variables

| Variable               | Required/Optional | Description |
|------------------------|-------------------|-------------|
| `perform_kraken`       | Optional          | Boolean flag to enable or disable Kraken2 classification of unmapped reads. Default is `false`. |
| `kraken_db_path`       | Optional          | Path to the Kraken2 database used for taxonomic classification of unmapped reads, if empty, a default database will be installed. |



## Filtering Variables

| Variable               | Required/Optional | Description |
|------------------------|-------------------|-------------|
| `perform_doublet_detection` | Optional          | Boolean flag to enable or disable consensus doublet detection (Scrublet + scDblFinder, via [Demuxafy](https://demultiplexing-doublet-detecting-docs.readthedocs.io/)). Detection runs on the cell-called matrix, and the calls reach `.obs["doublet_status"]` of every published AnnData object — independently of `ambient_rna_remover`, though CellSweep additionally projects them onto its UMAP comparison plot when it runs. Default is `true`. |
| `perform_doublet_filtering` | Optional          | Boolean flag to filter the detected (consensus) doublets out of the count matrices instead of only annotating them. The filtered matrix replaces the original for everything after it, so CellSweep denoises and `<outdir>/anndata/` publishes a matrix the doublets are gone from. Requires `perform_doublet_detection`; setting it without aborts the run at launch. Default is `false`. |
| `doublet_consensus_method`  | Optional          | Method passed to Demuxafy's `Combine_Results.R --method` to combine the Scrublet and scDblFinder calls. Default is `"AnySinglet"` (a cell is only called a doublet if *all* tools agree -- i.e. an intersection of doublet calls). See the [Demuxafy combining-results docs](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/CombineResults.html) for other options (`MajoritySinglet`, `AtLeastHalfSinglet`, `AnyDoublet`). |
| `demuxafy_sif`              | Optional          | Path to a pre-installed Demuxafy singularity image (`Demuxafy.sif`). Demuxafy is only distributed as a prebuilt `.sif` (~7.5GB, see the [Demuxafy installation docs](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/Installation.html)), not as a pullable Docker image. If unset, it is downloaded and md5-checked automatically the first time it's needed. Provide a shared/pre-installed path to skip re-downloading it on every fresh run. |
| `ambient_rna_remover`  | Optional          | Specifies the ambient RNA filtering software, must be one of the following: `"cellsweep"` or `"none"`. Default is set to `"cellsweep"`. |
| `h5ad_compression`     | Optional          | HDF5 compression for the objects in `<outdir>/anndata/`: `"gzip"` (default), `"lzf"` or `"none"`. Compression dominates `MTX_TO_H5AD`'s runtime on a raw matrix; `"lzf"` writes several times faster for roughly 15% more disk, but is only readable through h5py — fine for scanpy, not for R-side readers. |
| `gpu_cluster_options`  | Optional          | Scheduler flags used to request a GPU, applied only when running with `-profile gpu`. Defaults to the SLURM syntax `--partition=gpu --gres=gpu:1g.10gb:1`; set to `null` on schedulers that take no such flags, or override with your site's equivalent. |


## Metacells and the filtering report

With `perform_metacells = true`, each sample's cell-called AnnData object is grouped into
[Metacell2](https://github.com/tanaylab/metacells) metacells, and `<outdir>/filtering_report.html`
lets you decide, metacell by metacell and gene by gene, what goes into downstream analysis. The
report exports that decision as a `metacell_selection.json`, which a second run applies.

| Variable | Required/Optional | Description |
|----------|-------------------|-------------|
| `perform_metacells` | Optional | Group the called cells into metacells and build `filtering_report.html`. Needs `run_method` `"standard"` or `"post_mapping"`. Default is `false`. |
| `mc2_target_metacell_size` | Optional | Metacell2's target number of cells per metacell. Default is `96`. |
| `mc2_random_seed` | Optional | Seed for every random step of Metacell2. It is fixed so that a rerun reproduces the assignment a selection was made on. Default is `42`. |
| `mc2_min_cell_umis` / `mc2_max_cell_umis` | Optional | Cells outside this UMI range are not grouped. They are reported as the `__excluded__` group. Defaults are `200` and `null` (no upper bound). |
| `mc2_max_excluded_genes_fraction` | Optional | Metacell2's cap on the fraction of a cell's UMIs in excluded genes. The mitochondrial genes are excluded genes, so this cap would drop the most mito-rich cells *before* grouping and hide them from the report. Default is `null` (off). |
| `mc2_excluded_gene_patterns` | Optional | Extra genes kept out of the grouping (space-separated regular expressions). The mitochondrial genes are always excluded. Default is `""`. |
| `mc2_lateral_gene_patterns` | Optional | Genes marked *lateral* (e.g. cell-cycle genes). They do not drive the grouping but are still used to detect outliers. Default is `""`. |
| `mc2_min_cells` | Optional | Fewest groupable cells a sample needs. A sample with fewer is shown in the report as skipped, with the reason. Default is `300`. |
| `mc2_skip_subsampled` | Optional | Skip the `_subsampled_starsolo` and `_geneext_subsampled_starsolo` runs. Default is `true`. |
| `gene_annotation` | Optional | TSV with a header line, giving protein domains (PFAM) per gene, transcript or protein id. eggNOG-mapper's `.emapper.annotations` works as it is (see below). Without it, the report's PFAM panel is hidden. |
| `gene_annotation_id_col` | Optional | Identifier column of `gene_annotation`. A leading `#` in the header is ignored. Default is `"gene_id"`. |
| `gene_annotation_pfam_col` | Optional | Domain column of `gene_annotation`. Domains are separated by `,`, `;`, `\|` or spaces, and `-` means none. Default is `"pfam"`. |
| `mcfilter_max_mito_pct` | Optional | The report's starting threshold on mitochondrial UMIs per metacell (%). `null` starts it switched off. Default is `20`. |
| `mcfilter_max_doublet_pct` | Optional | Starting threshold on consensus doublets per metacell (%). Default is `50`. |
| `mcfilter_max_alpha_hat` | Optional | Starting threshold on the mean CellSweep ambient fraction per metacell. Default is `null`. |
| `mcfilter_min_gene_umis` | Optional | Starting minimum of total UMIs per gene. Default is `10`. |
| `mcfilter_pfam_presets` | Optional | Domain regular expressions offered as one-click presets in the report, separated by `;`. Default is `"^Ribosomal;^Histone"`. |
| `metacell_selection` | Optional | `metacell_selection.json` exported from the report. Rerun with `-resume` to write the filtered matrices. Needs `perform_metacells`. |
| `metacell_gene_mode` | Optional | `"drop"` removes the excluded genes from the final matrices. `"flag"` keeps every gene and records the decision in `var["pass_filter"]` and `var["exclusion_reasons"]`, which is convenient when samples are integrated later. Default is `"drop"`. |


## External Pipeline Variables

| Variable               | Required/Optional | Description |
|------------------------|-------------------|-------------|
| `perform_cellranger`   | Optional          | Boolean flag to enable or disable the CellRanger pipeline. Default is `false`. |
| `rhapsody_installation` | Optional          | Path to the BD-Rhapsody pipeline installation folder. |
| `splitpipe_installation`| Optional        | Path to the split-pipe installation folder. |
| `splitpipe_conda_env`   | Required (if `splitpipe_installation` is provided)         | Path to the split-pipe conda environment created by following these [instructions](docs/INSTALLATION_EXTERNAL_PIPELINES.md), required if running the split-pipe pipeline. |


## Execution Reports

Enabled for every run in [`nextflow.config`](../nextflow.config); there is nothing to
configure. Each run writes the following into `<outdir>/pipeline_info/`, suffixed with
the run's timestamp so repeated runs into the same directory accumulate rather than
overwrite each other:

| File | Contents |
|------|----------|
| `execution_trace_<timestamp>.txt`    | One row per task: requested and used cpus, memory and walltime, exit status, attempt, and bytes read |
| `execution_report_<timestamp>.html`  | Nextflow's own single-run resource report |
| `execution_timeline_<timestamp>.html`| Nextflow's own single-run timeline |
| `run_config_<timestamp>.txt`         | Every parameter the run used, written by `SAVE_RUN_CONFIG` |
| `samplesheet_<timestamp>.csv`        | Copy of the input samplesheet |

The trace's `fields` list is written out explicitly rather than left to Nextflow's
default, because the default records only the resources a task *used* — efficiency is
used ÷ requested, so the requested columns must be present too.

**Keep these files.** They are the input to `bin/resource_efficiency.py`, which turns
them into per-process efficiency reports and retuned resource requests. The more runs
accumulate, the better it can predict resource use on larger or smaller datasets. See
[Resource Tuning](RESOURCE_TUNING.md).

To disable reporting for a particular run, add to your custom configuration file:

```groovy
trace.enabled    = false
report.enabled   = false
timeline.enabled = false
```
