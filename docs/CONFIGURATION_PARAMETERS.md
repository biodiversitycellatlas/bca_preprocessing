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
| `protocol`              | __Required__     | Specifies the sequencing technology used (must be one of the following: `"oak_seq"`, `"10xv3"`, `"parse_biosciences_WT_mini"` or `"parse_biosciences_WT"`,     `"bd_rhapsody"`, `"sciRNAseq3"` , `"ultima_genomics"` or `"seqspec"`). |
| `ref_fasta`            | __Required__      | Path to the genome FASTA file used for mapping reads. |
| `ref_gtf`              | __Required__      | Path to the GTF/GFF file formatted for STARsolo. |
| `ref_gtf_alt`          | Optional          | Path to the GTF/GFF file formatted specifically for analysis with Parse Biosciences / CellRanger pipeline. Defaults to the same path as `ref_gtf`. |
| `seqspec_file`         | Optional          | Path to the seqspec file. |



## Mapping Variables

| Variable               | Required/Optional | Description |
|------------------------|-------------------|-------------|
| `mapping_software`     | Optional          | Software used to map reads (must be one of the following: `"starsolo"`, `"alevin"` or `"both"`). Default set to `"starsolo"`. | 
| `mt_contig`            | Optional          | Name of the mitochondrial contig in the reference annotation, used to calculate mtDNA content. Default set to `"chrM M MT"`. |
| `star_index`           | Optional          | |
| `star_solocellfilter`  | Optional          | |



## FeatureCounts Variables

| Variable               | Required/Optional | Description |
|------------------------|-------------------|-------------|
| `perform_featurecounts`  | Optional        | Boolean flag to enable or disable calculation of mtDNA & rRNA percentages. Default is `true`. |
| `grep_rrna`            | Optional          | String used to grep ribosomal RNA (rRNA) reads from annotations. Default set to `"rRNA"`| 



## Gene Extension Variables

| Variable               | Required/Optional | Description |
|------------------------|-------------------|-------------|
| `perform_geneext`      | Optional          | Boolean flag to enable or disable the gene extension step in preprocessing. Default is `false`. |



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
| `splitpipe_installation`| Optional         | Path to the split-pipe installation folder, that can be used as a control. |
| `splitpipe_conda_env`   | Optional         | Path to the split-pipe conda environment created by following these [instructions](docs/INSTALLATION_EXTERNAL_PIPELINES.md), required if running the split-pipe pipeline. |
| `rhapsody_installation`| Optional          | Path to the BD-Rhapsody pipeline installation folder, that can be used as a control. |
