# BCA Pre-Processing Pipeline

---

## Table of Contents

1. [Overview](#overview)
2. [Installation & Setup](#installation--setup)
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
Parse Biosciences data will be demultiplexed depending on the groups parameter, which seperates possible different techniques or samples within the same plate. After demultiplexing, it is mapped using the official split-pipe code from Parse Biosciences, to offer a comparation between their data processing platform and the results of our pipeline. BD Rhapsody does not require demultiplexing, and is therefore sent straight to the mapping using STARsolo. 

![pipeline](/img/Preprocs_Pipeline.png)

---

## Installation & Setup

In order to run the pipeline, you must have [Nextflow](https://www.nextflow.io/) and [Conda](https://anaconda.org/) installed.

### Conda environments
The default conda environment is called 'bca_env', which is automatically created and activated using the .yml file upon running the pipeline. For a few subprocesses, different conda environments were created as they were conflicting with certain versions of packages. All of them, except for spipe (see below) are installed and activated automatically upon execution using the same approach. 

To manually activate any of the conda environments:
```
# Create environment
conda env create -n bca_env -f bca_env.yaml

# Activate environment
conda activate bca_env
```


### Installing external/commercial packages

### Cell Ranger
Followed the installation guide on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-in), and downloaded Cell Ranger version 9.0.1.
```
# Downloading Cell Ranger using wget
wget -O cellranger-9.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.gz?Expires=1742511989&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=cr3Sw2Q~GjIzwmorEmEGgU7eKfFqTfz14Gd5Dt329DuPX549cfXVEsLJo7oq3xcijzJpNIFbbxDL7JyP0-2LA4GyQhyyKvEjoBuHNjDkdB8qo8lQ4yJ57oThwz8kTPvc3NVBy8jHQYfA8ywcz4dWrrt0--K5bnp4OMEi6A0QTFaUUfINjPjaFqUoOrMuYtJdOSBcj11h9xI~eVbF4d~-bB4zKiqwjnKlW0aJ2K6jWt7Ho22V8YXSN3o70hJOPnklf2hVpfvKmKqfGFC1IeHGEOvoEcBAPUI9qnmoIKAn3FYvIhnqtkOrAH42naArGmFfiSGl8iaOgBOBYCIQUohMmA__"

# Unpack
tar -xzvf cellranger-9.0.1.tar.gz

# Add Cell Ranger to PATH
export PATH=/path/to/cellranger-9.0.1:$PATH
```

### BD Rhapsody

First, download the ['cwl' folder](https://bitbucket.org/CRSwDev/cwl/src/master/) of the BD Single-Cell Multiomics Software, hosted on Bitbucket. This repository contains the CWL and YML files required to run the BD Rhapsody pipeline locally. Then, create and activate the conda environment 'bd_pipe', add the cwlref-runner to the PATH variable and finally pull the docker image. In our case, we pulled the container using Apptainer. 

```
# Add cwlref-runner to PATH
export PATH=$PATH:/path/to/bd_pipe/env/lib/python3.13/site-packages

# Pull docker image using Apptainer
apptainer pull docker://bdgenomics/rhapsody
```

#### Parse Biosciences
```
conda create –name spipe 
bash ./install_dependencies_conda.sh -i -y 
pip install . --no-cache-dir 
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

## Usage

### Required files & parameters
- [ ] /fastq/ folder with raw FASTQ files.
- [ ] Genome annotation files (FASTA and GFF/GTF files)
- [ ] Configuration file

### Optional files & parameters
- [ ] spec.yaml - created using seqspec

### Running the Pipeline
```
sbatch submit_nextflow.sh main.nf
```
