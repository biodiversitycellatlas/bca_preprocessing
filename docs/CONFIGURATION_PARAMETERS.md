# Custom Configuration Variables


Within each custom configuration file the following variables can be defined:

## Table of Contents

1. [Base Variables](#base-variables)
2. [STAR Variables](#star-variables)
3. [FeatureCounts Variables](#featurecounts-variables)
4. [Gene Extension Variables](#gene-extension-variables)
5. [Taxonomic-classification Variables](#taxonomic-classification-variables)
6. [Filtering Variables](#filtering-variables)
7. [External Pipeline Variables](#external-pipeline-variables)

---


## Base Variables

| Variable               | Required/Optional | Description |
|------------------------|-------------------|-------------|
| `input`                | __Required__      | Path to the samplesheet. |
| `outdir`               | __Required__      | Path to the results/output directory; must exist before running. |
| `bc_whitelist`         | __Required__      | Path or link to the barcode whitelist file(s), multiple ones separated by whitespace. If links are given, they are automatically downloaded (and unzipped if applicable) for any protocol. Not used by `"marsseq_v1"`/`"marsseq_v2"`, which run without a whitelist.|
| `protocol`              | __Required__     | Specifies the sequencing technology used (must be one of the following: `"oak_seq"`, `"10xv1"`, `"10xv2"`, `"10xv3"`, `"10xv4"`, `"parse_biosciences_WT_mini"` or `"parse_biosciences_WT"`,     `"bd_rhapsody_v1"`, `"bd_rhapsody_enhancedbeads"`, `"sciRNAseq3"` , `"ultima_genomics"`, `"marsseq_v1"`, `"marsseq_v2"` or `"seqspec"`). |
| `ref_fasta`            | __Required__      | Path to the genome FASTA file used for mapping reads. |
| `ref_gtf`              | __Required__      | Path to the GTF/GFF file formatted for STARsolo. |
| `ref_gtf_alt`          | Optional          | Path to the GTF/GFF file formatted specifically for analysis with Parse Biosciences / CellRanger pipeline. Defaults to the same path as `ref_gtf`. |
| `run_method`             | Optional          | Method of running the pre-processing pipeline, demonstrated in the [pipeline diagram](img/Preprocs_Pipeline.png), one of `"standard"`, `"geneext_only"`, `"external_pipeline_only"` or `"post_mapping"`. Default is set to `"standard"`. `"post_mapping"` picks a finished run up at the point mapping ended and redoes everything after it. |
| `previous_outdir`      | Optional          | Only used with `run_method = "post_mapping"`: the results directory the published mapping results are read back from. Defaults to `outdir`, i.e. re-running in place. |
| `mapping_software`     | Optional          | Software used to map reads (must be one of the following: `"starsolo"`, `"alevin"`, `"both"/"alevin_starsolo"` or `"alevin_subsampled_starsolo"`). Default set to `"starsolo"`. |
| `cellfilter_method`    | Optional          | How cells are called: `"star_solocellfilter"` (each mapper's own cell filter: STARsolo's `--soloCellFilter`, alevin-fry's knee), `"second_derivative"` or `"manual_cutoff"`. Default set to `"second_derivative"`. The latter two apply to both mappers: they filter the count matrix on a UMI cutoff and recompute every cell-dependent statistic in the report on the retained cells, differing only in where that cutoff comes from. See [Cell calling](#cell-calling). |
| `perform_demultiplexing` | Optional        | Boolean flag to enable or disable demultiplexing of the FASTQ files, where applicable. Default is `true`. |
| `seqspec_file`         | Optional          | Path to the seqspec file. |
| `subsample_nreads`     | Optional          | The size (number of reads) of the subset used to map to STARsolo, in case the parameter `mapping_software = alevin_subsampled_starsolo`. Default set to `100000000` reads. |
| `alevin_usa_counts`    | Optional          | Which of alevin-fry's USA blocks are summed into each gene's count before ambient-RNA removal and doublet detection: `"SUA"` (spliced + unspliced + ambiguous), `"SA"` or `"S"`. Default set to `"SUA"`. See [USA counts](#usa-counts). |



## STAR Variables

| Variable               | Required/Optional | Description |
|------------------------|-------------------|-------------|
| `star_index`           | Optional          | Path to the pre-generated STAR index. By default the STAR index is created within the pipeline.|
| `star_genomeSAindexNbases` | Optional         | Lenght of the SA pre-indexing string in STAR. See [protocol-specific defaults](../conf/seqtech_parameters.config) set in the seqtech_paramaters.config file. |
| `star_genomeSAsparseD`    | Optional       | Suffix array sparsity in STAR.  See [protocol-specific defaults](../conf/seqtech_parameters.config) set in the seqtech_paramaters.config file. |
| `star_solocellfilter`  | Optional          | Cell filtering type and parameters used by STAR. See [protocol-specific defaults](../conf/seqtech_parameters.config) set in the seqtech_paramaters.config file.|
| `star_soloFeatures`    | Optional         | Genomic features for which the UMI counts per Cell Barcode are collected. See [protocol-specific defaults](../conf/seqtech_parameters.config) set in the seqtech_paramaters.config file. |
| `star_outFilterScoreMin` | Optional       | Alignment will be output only if its score is higher than or equal to this value normalized by read length. See [protocol-specific defaults](../conf/seqtech_parameters.config) set in the seqtech_paramaters.config file. |
| `star_outSAMunmapped`     | Optional      | Output of unmapped reads in the SAM format. See [protocol-specific defaults](../conf/seqtech_parameters.config) set in the seqtech_paramaters.config file. |
| `star_outSAMattributes` | Optional        | String of desired SAM attributes. See [protocol-specific defaults](../conf/seqtech_parameters.config) set in the seqtech_paramaters.config file. |
| `star_generateBAM`      | Optional        | Boolean to determine if a BAM file should be generated. Will automatically adjust the star_outSAMattributes. Default is set to `true`. |
| `star_soloTypestring` | Optional          | String of defining soloType and barcode structure. See [protocol-specific defaults](../conf/seqtech_parameters.config) set in the seqtech_paramaters.config file. |
| `star_soloCBmatchWLtype` | Optional       | Type of matching the cell barcodes to the barcode whiteList. See [protocol-specific defaults](../conf/seqtech_parameters.config) set in the seqtech_paramaters.config file. |
| `star_soloUMIfiltering` | Optional        | Type of UMI filtering. See [protocol-specific defaults](../conf/seqtech_parameters.config) set in the seqtech_paramaters.config file. |
| `star_soloMultiMappers` | Optional        | Counting method for reads mapping to multiple genes. See [protocol-specific defaults](../conf/seqtech_parameters.config) set in the seqtech_paramaters.config file. |
| `star_soloUMIdedup` | Optional         | Type of UMI deduplication algorithm. See [protocol-specific defaults](../conf/seqtech_parameters.config) set in the seqtech_paramaters.config file. |
| `star_clipAdapterType` | Optional         | Type of adapter clipping. See [protocol-specific defaults](../conf/seqtech_parameters.config) set in the seqtech_paramaters.config file. |
| `star_limitBAMsortRAM` | Optional         | STAR's own ceiling on the BAM sort buffer, in bytes. Default `0`, which makes the pipeline derive it from `task.memory` minus the resident genome index. Set a number only to pin it, e.g. to the figure a STAR error asks for. Note this is *not* the scheduler allocation — see [Resource Tuning](RESOURCE_TUNING.md#tool-internal-memory-ceilings). |
| `star_limitGenomeGenerateRAM` | Optional  | STAR's own ceiling for index construction, in bytes. Default `null`, which derives it from `task.memory`. A pinned value neither grows with a larger resource tier nor shrinks with a smaller one. |


## Alevin-fry Variables

The geometries use salmon's custom geometry syntax, `"<read>[<ranges>]"`, with 1-based inclusive ranges that are concatenated in the order they are listed. Read `1` is the CB/UMI FASTQ and read `2` the cDNA FASTQ, regardless of which of R1/R2 carries the cDNA for the protocol.

| Variable               | Required/Optional | Description |
|------------------------|-------------------|-------------|
| `alevin_bc_geometry`   | Optional          | Position(s) of the cell barcode, e.g. `"1[1-16]"` or `"1[1-9,14-22,27-35]"`. See [protocol-specific defaults](../conf/seqtech_parameters.config) set in the seqtech_paramaters.config file. |
| `alevin_umi_geometry`  | Optional          | Position of the UMI, e.g. `"1[17-28]"`. See [protocol-specific defaults](../conf/seqtech_parameters.config) set in the seqtech_paramaters.config file. |
| `alevin_read_geometry` | Optional          | Position of the biological sequence, `"2[1-end]"` for every protocol. See [protocol-specific defaults](../conf/seqtech_parameters.config) set in the seqtech_paramaters.config file. |


## FeatureCounts Variables

| Variable               | Required/Optional | Description |
|------------------------|-------------------|-------------|
| `perform_featurecounts`  | Optional        | Boolean flag to enable or disable calculation of mtDNA & rRNA percentages. Default is `false`. |
| `mt_contig`            | Optional          | Name(s) of the mitochondrial contig(s) in the reference annotation, used to calculate mtDNA content. Multiple contigs can be given separated by whitespace, in which case reads on any of them count as mitochondrial. Default set to `"chrM M MT"`. |
| `grep_rrna`            | Optional          | String used to grep ribosomal RNA (rRNA) reads from annotations. Default set to `"rRNA"`|



## Gene Extension Variables

| Variable               | Required/Optional | Description |
|------------------------|-------------------|-------------|
| `perform_geneext`      | Optional          | Boolean flag to enable or disable the gene extension step after mapping. Default is `false`. |


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
| `perform_doublet_detection` | Optional          | Boolean flag to enable or disable consensus doublet detection (Scrublet + scDblFinder, via [Demuxafy](https://demultiplexing-doublet-detecting-docs.readthedocs.io/)). Detection runs on the cell-called matrix, and the calls are handed to CellSweep, which annotates them in `.obs["doublet_status"]` and projects them onto the UMAP comparison plot when cellsweep is enabled. Default is `true`. |
| `perform_doublet_filtering` | Optional          | Boolean flag to filter the detected (consensus) doublets out of the count matrices instead of only annotating them. Requires `perform_doublet_detection`; setting it without aborts the run at launch. Default is `false`. |
| `doublet_consensus_method`  | Optional          | Method passed to Demuxafy's `Combine_Results.R --method` to combine the Scrublet and scDblFinder calls. Default is `"AnySinglet"` (a cell is only called a doublet if *all* tools agree -- i.e. an intersection of doublet calls). See the [Demuxafy combining-results docs](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/CombineResults.html) for other options (`MajoritySinglet`, `AtLeastHalfSinglet`, `AnyDoublet`). |
| `demuxafy_sif`              | Optional          | Path to a pre-installed Demuxafy singularity image (`Demuxafy.sif`). Demuxafy is only distributed as a prebuilt `.sif` (~7.5GB, see the [Demuxafy installation docs](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/Installation.html)), not as a pullable Docker image. If unset, it is downloaded and md5-checked automatically the first time it's needed. Provide a shared/pre-installed path to skip re-downloading it on every fresh run. |
| `perform_cellbender`   | Optional          | Specifies the ambient RNA filtering software, must be one of the following: `"cellsweep"`, `"cellbender"` or `"none"`. Default is set to `"cellsweep"`. |
| `cellbender_extraargs` | Optional          | Provide extra arguments to the CellBender function as a string. Refer to the [CellBender manual](https://cellbender.readthedocs.io/en/latest/reference/index.html) for options. |
| `gpu_cluster_options`  | Optional          | Scheduler flags used to request a GPU, applied only when running with `-profile gpu`. Defaults to the SLURM syntax `--partition=gpu --gres=gpu:1g.10gb:1`; set to `null` on schedulers that take no such flags, or override with your site's equivalent. |


## External Pipeline Variables

| Variable               | Required/Optional | Description |
|------------------------|-------------------|-------------|
| `perform_cellranger`   | Optional          | Boolean flag to enable or disable the CellRanger pipeline. Default is `false`. |
| `rhapsody_installation` | Optional          | Path to the BD-Rhapsody pipeline installation folder. |
| `splitpipe_installation`| Optional        | Path to the split-pipe installation folder. |
| `splitpipe_conda_env`   | Required (if `splitpipe_installation` is provided)         | Path to the split-pipe conda environment created by following these [instructions](docs/INSTALLATION_EXTERNAL_PIPELINES.md), required if running the split-pipe pipeline. |
