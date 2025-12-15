# Custom Configuration Variables


Within each custom configuration file the following variables can be defined:

## Table of Contents

1. [Base Variables](#base-variables)
2. [Mapping Variables](#mapping-variables)
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
| `bc_whitelist`         | __Required__      | Path or link to the barcode whitelist file(s). If it's a link, it will be automatically downloaded and unzipped if applicable.|
| `protocol`              | __Required__     | Specifies the sequencing technology used (must be one of the following: `"oak_seq"`, `"10xv3"`, `"parse_biosciences_WT_mini"` or `"parse_biosciences_WT"`,     `"bd_rhapsody"`, `"sciRNAseq3"` , `"ultima_genomics"` or `"seqspec"`). |
| `ref_fasta`            | __Required__      | Path to the genome FASTA file used for mapping reads. |
| `ref_gtf`              | __Required__      | Path to the GTF/GFF file formatted for STARsolo. |
| `ref_gtf_alt`          | Optional          | Path to the GTF/GFF file formatted specifically for analysis with Parse Biosciences / CellRanger pipeline. Defaults to the same path as `ref_gtf`. |
| `run_method`             | Optional          | Method of running the pre-processing pipeline, demonstrated in the [pipeline diagram](img/Preprocs_Pipeline.png), currently either `"standard"` or `"geneext_only"`. Default is set to `"standard"`. |
| `perform_demultiplexing` | Optional        | Boolean flag to enable or disable demultiplexing of the FASTQ files, where applicable. Default is `true`. |
| `seqspec_file`         | Optional          | Path to the seqspec file. |



## Mapping Variables

| Variable               | Required/Optional | Description |
|------------------------|-------------------|-------------|
| `mapping_software`     | Optional          | Software used to map reads (must be one of the following: `"starsolo"`, `"alevin"` or `"both"`). Default set to `"starsolo"`. |
| `mt_contig`            | Optional          | Name of the mitochondrial contig in the reference annotation, used to calculate mtDNA content. Default set to `"chrM M MT"`. |
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


## FeatureCounts Variables

| Variable               | Required/Optional | Description |
|------------------------|-------------------|-------------|
| `perform_featurecounts`  | Optional        | Boolean flag to enable or disable calculation of mtDNA & rRNA percentages. Default is `false`. |
| `grep_rrna`            | Optional          | String used to grep ribosomal RNA (rRNA) reads from annotations. Default set to `"rRNA"`|



## Gene Extension Variables

| Variable               | Required/Optional | Description |
|------------------------|-------------------|-------------|
| `perform_geneext`      | Optional          | Boolean flag to enable or disable the gene extension step after mapping. Default is `false`. |


## 10x_saturate Variables

| Variable               | Required/Optional | Description |
|------------------------|-------------------|-------------|
| `perform_10x_saturate`      | Optional          | Boolean flag to enable or disable the 10x_saturate step after mapping. Default is `true`. |



## Taxonomic-classification Variables

| Variable               | Required/Optional | Description |
|------------------------|-------------------|-------------|
| `perform_kraken`       | Optional          | Boolean flag to enable or disable Kraken2 classification of unmapped reads. Default is `false`. |
| `kraken_db_path`       | Optional          | Path to the Kraken2 database used for taxonomic classification of unmapped reads, if empty, a default database will be installed. |



## Filtering Variables

| Variable               | Required/Optional | Description |
|------------------------|-------------------|-------------|
| `perform_cellbender`   | Optional          | Boolean flag to enable or disable removal of ambient RNA using CellBender. Default is `false`. |



## External Pipeline Variables

| Variable               | Required/Optional | Description |
|------------------------|-------------------|-------------|
| `perform_cellranger`   | Optional          | Boolean flag to enable or disable the CellRanger pipeline. Default is `false`. |
| `rhapsody_installation` | Optional          | Path to the BD-Rhapsody pipeline installation folder. |
| `splitpipe_installation`| Optional        | Path to the split-pipe installation folder. |
| `splitpipe_conda_env`   | Required (if `splitpipe_installation` is provided)         | Path to the split-pipe conda environment created by following these [instructions](docs/INSTALLATION_EXTERNAL_PIPELINES.md), required if running the split-pipe pipeline. |
| `scalerna_installation`  | Optional          | Path to the base of the ScaleRna pipeline installation folder. |
| `scalerna_scalePlex`     | Required (if `scalerna_installation` is provided)  | Boolean flag if the data is a ScalePlex dataset. |
| `scalerna_libStructure`  | Required (if `scalerna_installation` is provided) | Scale scRNA assay version, should be one of the following: "libV1.json", "libV1.1.json", or "libQuantumV1.0.json" |
