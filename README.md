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

Depending on the chosen sequencing technique, it handles the FASTQ files accordingly. 
Parse Biosciences data will be demultiplexed depending on the groups parameter, which seperates possible different techniques or samples within the same plate based on the wells. After demultiplexing, it is mapped using the official split-pipe pipeline from Parse Biosciences, to offer a comparation between their data processing platform and the results of our pipeline. The same comparisons are provided for BD Rhapsody (using BD-Rhapsody pipeline), OAKseq and 10x data (using Cell Ranger).   

![pipeline](/img/Preprocs_Pipeline.png)

---

## Installation

### Basic Installation

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


3. **Sequencing technique-specific requirements**


When working with OAK-seq or 10x Genomics data, the CellRanger whitelist must be decompressed:
```
gunzip /seq_techniques/oak_seq/barcodes_R2.txt.gz
```


### Optional: Installing commercial pipelines

#### Cell Ranger
Followed the installation guide on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-in), and downloaded Cell Ranger version 9.0.1.
```
# Downloading Cell Ranger using wget
wget -O cellranger-9.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.gz?Expires=1742511989&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=cr3Sw2Q~GjIzwmorEmEGgU7eKfFqTfz14Gd5Dt329DuPX549cfXVEsLJo7oq3xcijzJpNIFbbxDL7JyP0-2LA4GyQhyyKvEjoBuHNjDkdB8qo8lQ4yJ57oThwz8kTPvc3NVBy8jHQYfA8ywcz4dWrrt0--K5bnp4OMEi6A0QTFaUUfINjPjaFqUoOrMuYtJdOSBcj11h9xI~eVbF4d~-bB4zKiqwjnKlW0aJ2K6jWt7Ho22V8YXSN3o70hJOPnklf2hVpfvKmKqfGFC1IeHGEOvoEcBAPUI9qnmoIKAn3FYvIhnqtkOrAH42naArGmFfiSGl8iaOgBOBYCIQUohMmA__"

# Unpack
tar -xzvf cellranger-9.0.1.tar.gz
```

#### BD Rhapsody

First, download the ['cwl' folder](https://bitbucket.org/CRSwDev/cwl/src/master/) of the BD Single-Cell Multiomics Software, hosted on Bitbucket. This repository contains the CWL and YML files required to run the BD Rhapsody pipeline locally. Secondly, pull the docker image using Docker or Apptainer (see code for Apptainer example). 

```
# Pull docker image using Apptainer
apptainer pull docker://bdgenomics/rhapsody
```

#### Parse Biosciences
In order to incorporate the commercial Parse Biosciences pipeline (also called split-pipe), we first have to install the code from the personal account page on the Parse Biosciences website, and then create a new conda environment. 

```
# Move to new split-pipe installation
cd /installation_pipeline/

# Create a new conda environment
conda create –name spipe 
bash ./install_dependencies_conda.sh -i -y 
pip install . --no-cache-dir 

# Activate environment
conda activate spipe
```

We encountered a few errors durring the creation of the spipe environment. Below, there are a few fixes to common set up problems of spipe: 
```
# Pip not installed
conda install anaconda::pip=23.3.1 

# Failed building wheel for louvain
conda install -c conda-forge python-igraph 
pip install cmake 

# AttributeError: module 'numpy' has no attribute 'NAN'. Did you mean: 'nan'? 
# In the file 'utils.py', in line 'def report_percent_str(num, den=1, round_to=2, zero=np.NAN, perchar=True): ' 
# replace NAN by nan. 
# After saving, rerun:
pip install . --no-cache-dir 
```

To test if the installation of split-pipe was successful:
```
split-pipe -h 
```

---

## Setup

### 1. Create a new Nextflow configuration file

To set the custom parameters for each run, a nextflow configuration file is created that extends the general `nextflow.config` file in the base of this repository. The custom configuration files can be stored within the `/config/` directory, and should be based upon the example config file in [`/examples/config/`](https://github.com/biodiversitycellatlas/bca_preprocessing/blob/main/examples/config/example.config). If you have a multiple-species experiment, one configuration file per species must be created in order to analyze the data with the corresponding genome annotation files. 

Within each custom configuration file the following variables will be defined: 


| Variable                | Description |
|------------------------|-------------|
| `baseDir`              | Path to the BCA-preprocessing codebase; used to access sequencing-specific scripts and pipeline components. |
| `resDir`               | Path to the results/output directory containing FASTQ files and reference genome files; must exist before running. |
| `seqTech`              | Specifies the sequencing technology used (e.g., `"oak_seq"`, `"parse_biosciences"`, `"bd_rhapsody"`). |
| `annot_type`           | Annotation file type used for reference (either `"GTF"` or `"GFF"`). |
| `ref_star_gtf`         | Path to the GTF/GFF file formatted for STARsolo. |
| `ref_parse_gtf`        | Path to the GTF/GFF file formatted specifically for Parse Biosciences requirements. |
| `ref_fasta`            | Path to the genome FASTA file used for mapping reads. |
| `perform_geneext`      | Boolean flag to enable or disable the gene extension step in preprocessing. |
| `mt_contig`            | Name of the mitochondrial contig in the reference annotation, used to calculate mtDNA content. |
| `grep_rrna`            | String used to grep ribosomal RNA (rRNA) reads from annotations. |
| `perform_kraken`       | Boolean flag to enable Kraken2 classification of unmapped reads. |
| `kraken_db_path`       | Path to the Kraken2 database used for taxonomic classification of unmapped reads, if empty, a new database will be installed. |
| `cellranger_dir`       | (optional) Path to the Cell Ranger software directory (used for 10x Genomics & OAK-seq data). |
| `bdrhap_pipeline_dir`  | (optional) Path to the BD Rhapsody pipeline directory. |
| `parsebio_pipeline_dir`| (optional) Path to the Parse Biosciences pipeline directory. |


### 2. Add custom configuration file as a new profile

After creating a new configuration file, it has to be added as a profile in the [`nextflow.config`](https://github.com/biodiversitycellatlas/bca_preprocessing/blob/main/nextflow.config) file. Define an unique name and set the path, for example within the /config/ directory. 
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
    'test_config' {
        includeConfig 'config/test.config'
    }
    ...
}
```

### 3. Edit submit_nextflow.sh

```
## file: submit_nextflow.sh

# Include both slurm & custom config profile
nextflow run -profile slurm,test_config -ansi-log false "$@" & pid=$! 
```

---

## Usage

### Required files & parameters
- [ ] Nextflow configuration file
- [ ] /fastq/ folder with raw FASTQ files.
- [ ] /genome/ folder with genome annotation files (FASTA and GFF/GTF files)

### Optional files & parameters
- [ ] spec.yaml - created using seqspec

### Running the Pipeline

> [!WARNING]
> The pipeline must be run in the conda base environment, it cannot activate the different environments properly with any prior environment activation. It should have access to run `nextflow` and `conda` in the commandline.

```
# ((optional: load Nextflow module OR have local Nextflow installation)) 
module load Nextflow/24.04.3

# Submit pipeline to SLURM queue
sbatch submit_nextflow.sh main.nf
```
