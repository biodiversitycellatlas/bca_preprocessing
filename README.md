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

Clone the bca_preprocessing repository:
```
git clone https://github.com/biodiversitycellatlas/bca_preprocessing.git
```
In order to run the pipeline, you must have [Conda](https://anaconda.org/) installed.

When working with OAK-seq or 10x Genomics data, the CellRanger whitelist must be uncompressed:
```
gunzip /seq_techniques/oak_seq/barcodes_R2.txt.gz
```

### Conda environments
The main environment of the pipeline is 'bca_env', which must be created and activated before running the pre-processing pipeline. This is the only conda environment that needs to be created prior to running the pipeline, as the other environments within the /env/ directory are automatically created and activated by Nextflow. 

To create & activate the 'bca_env' conda environment:
```
# Create environment
conda env create -n bca_env -f ./env/bca_env.yml

# Activate environment
conda activate bca_env
```

### Installing commercial pipelines

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

1. Create a new Nextflow configuration file

To set the custom parameters for each run, a nextflow configuration file is created that extends the general nextflow.config file in the base of this repository. The custom configuration files can be stored within the /config/ directory, and will be based upon the example config file in'./examples/config/'. If you have a multiple-species experiment, one configuration file per species must be created in order to analyze the data with the corresponding genome annotation files. 

2. Add custom configuration file as a new profile

After creating a new configuration file, it has to be added as a profile in the 'nextflow.config' file. Define an unique name and set the path, for example within the /config/ directory. 
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

3. Edit submit_nextflow.sh

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
```
# Activate environment
conda activate bca_env

# Submit pipeline to SLURM queue
sbatch submit_nextflow.sh main.nf
```
