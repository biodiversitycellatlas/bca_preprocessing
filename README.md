# BCA Pre-Processing Pipeline

---

## Table of Contents

1. [Overview](#overview)
2. [Installation](#installation)
3. [Setup](#setup)
3. [Usage](#usage)
4. [Output & Logs](#output--logs)

---

## Overview
This nextflow pipeline is designed to pre-process single-cell and single-nucleus RNA-seq data. It accepts FASTQ files from multiple sequencing platforms, at the moment being: 
- Parse Biosciences
- BD Rhapsody
- OAK-seq 
- 10x Genomics
- sci-RNA-seq3
- and others, when providing a [seqspec](https://github.com/pachterlab/seqspec) file

Depending on the chosen sequencing technique, it handles the processing of the FASTQ files accordingly. For an in-depth explenation of the different steps for each of the sequencing methods, you can read [this page](subworkflows/README.md). Whenever possible, we compared our results to a commercial pre-processing pipeline for that sequencing technique. For example, comparing our Parse Biosciences results to the official split-pipe pipeline from Parse Biosciences. While we cannot provide this commercial software directly, you can install it yourself (e.g. by following [these instructions](assets/README.md)), and provide a path where the installation is located. This way, it will be executed alongside of the BCA pre-processing pipeline.

The pipeline will produce the following output files:
- Raw & Filtered count matrices (intronic, exonic & full gene) from Mapping step​
- Filtered count matrices (h5 files) ​from Filtering step​
- Summary report (MultiQC) & Mapping statistics table

![pipeline](/img/Preprocs_Pipeline.png)

---

## Installation

1. **Clone the bca_preprocessing repository:**
```
git clone https://github.com/biodiversitycellatlas/bca_preprocessing.git
```

2. **Conda & Nextflow**

In order to run the pipeline, you must have [Conda](https://anaconda.org/) and [Nextflow](https://www.nextflow.io/docs/latest/index.html) installed. When working on a HPC, there might be a module system available to use instead, simplifying the use of different software. 
To see which modules are available and how to load them: 
```
# List currently loaded modules
module list

# Check for available Nextflow modules
module avail Nextflow

# Load Nextflow module (make sure to change version accordingly!)
module load Nextflow/24.04.3
```


To test the installations, please try:
```
# Verify conda installation
conda -h

# Verify nextflow installation
nextflow -h
``` 

---

## Setup

### 1. Create a samplesheet

| sample                 |  fastq_cDNA | fastq_BC_UMI | expected_cells | p5                 |  p7 | rt |
|------------------------|-------------|--------------|-------|-----------------|-------------|--------------|


### 2. Edit (or create) the Nextflow configuration file

To set the custom parameters for each run, a nextflow configuration file is created that extends the general `nextflow.config` file in the base of this repository. The custom configuration files can be stored within the `/conf/` directory, and should be based upon the example config file in [`conf/example.config`](https://github.com/biodiversitycellatlas/bca_preprocessing/blob/main/conf/example.config). If you have a multiple-species experiment, one configuration file per species must be created in order to analyze the data with the corresponding genome annotation files. 

Within each custom configuration file the following variables can be defined: 


| Variable               | Required/Optional | Description |
|------------------------|-------------------|-------------|
| `input`                | __Required__         | Path to the samplesheet. |
| `output_dir`           | __Required__          | Path to the results/output directory; must exist before running. |
| `protocol`              | __Required__          | Specifies the sequencing technology used (must be one of the following: `"oak_seq"`, `"10xv3"`, `"parse_biosciences"`,     `"bd_rhapsody"`, `"sciRNAseq3"` or `"seqspec"`). |
| `parsebio_groups`       | Optional          | Required for Parse Biosciences data, to split the FASTQ files by well. Should be a nested list, with their desired name and range of wells (example: `[['group1', 'A1-A3'], ['group2', 'A4-A7'], ...]`)  |
| `ref_fasta`            | __Required__          | Path to the genome FASTA file used for mapping reads. |
| `ref_gtf`              | __Required__          | Path to the GTF/GFF file formatted for STARsolo. |
| `ref_parse_gtf`        | Optional              | Path to the GTF/GFF file formatted specifically for analysis with Parse Biosciences pipeline. Defaults to the same path as `ref_gtf`. |
| `seqspec_file`         | Optional              | Path to the seqspec file. |
| `mt_contig`            | Optional          | Name of the mitochondrial contig in the reference annotation, used to calculate mtDNA content. Default set to `"^MT"` |
| `grep_rrna`            | Optional          | String used to grep ribosomal RNA (rRNA) reads from annotations. Default set to `"rRNA"`|
| `mapping_software`     | Optional          | Software used to map reads (must be one of the following: `"starsolo"`, `"alevin"` or `"both"`). Default set to `"starsolo"`. | 
| `perform_geneext`      | Optional          | Boolean flag to enable or disable the gene extension step in preprocessing. Default is `true`. |
| `perform_featurecounts`  | Optional          | Boolean flag to enable or disable calculation of mtDNA & rRNA percentages. Default is `true`. |
| `perform_kraken`       | Optional          | Boolean flag to enable or disable Kraken2 classification of unmapped reads. Default is `false`. | 
| `perform_cellbender`   | Optional          | Boolean flag to enable or disable removal of ambient RNA using CellBender. Default is `false`. | 
| `kraken_db_path`       | Optional          | Path to the Kraken2 database used for taxonomic classification of unmapped reads, if empty, a default database will be installed. |
| `cellranger_dir`       | Optional          | Path to the Cell Ranger software directory (used for 10x Genomics & OAK-seq data). |
| `bdrhap_pipeline_dir`  | Optional          | Path to the BD Rhapsody pipeline directory. |
| `parsebio_pipeline_dir`| Optional          | Path to the Parse Biosciences pipeline directory. |



### 3. Optional: Add custom configuration file as a new profile

After creating a new configuration file, it can be added as a profile in the [`nextflow.config`](https://github.com/biodiversitycellatlas/bca_preprocessing/blob/main/nextflow.config) file. Define an unique name and set the path, for example within the /conf/ directory. 

```
## file: nextflow.config

profiles {
    slurm {
        process {
            queue = 'genoa64'
            executor = "slurm"
            clusterOptions = '--qos=short'
        }
    }
    /* ======= !!! EDIT BELOW TO INCLUDE YOUR CUSTOM CONFIG FILES !!! ======= */
    'custom_parameters' {
        includeConfig 'conf/custom_parameters.config'
    }
    /* ======= !!! EDIT ABOVE TO INCLUDE YOUR CUSTOM CONFIG FILES !!! ======= */
}
```


---

## Usage

### Pre-requisites:
- [ ] Configured the custom config file (config/custom_parameters.config)
- [ ] Added custom config as profile in the main config file (nextflow.config)

### Running the Pipeline

> [!WARNING]
> The pipeline must be run in the conda base environment, it cannot activate the different environments properly with any prior environment activation. It should have access to run `nextflow` and `conda` in the commandline.

```
# ((optional: load Nextflow module OR have local Nextflow installation)) 
module load Nextflow/24.04.3

# Submit pipeline to SLURM queue
sbatch submit_nextflow.sh main.nf
```
