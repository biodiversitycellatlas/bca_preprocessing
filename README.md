# Nextflow pipeline for pre-processing single-cell and single-nucleus data

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

Depending on the chosen sequencing technique, it handles the FASTQ files accordingly. 
Parse Biosciences data will be demultiplexed depending on the groups parameter, which seperates possible different techniques or samples within the same plate. After demultiplexing, it is mapped using the official split-pipe code from Parse Biosciences, to offer a comparation between their data processing platform and the results of our pipeline. BD Rhapsody does not require demultiplexing, and is therefore sent straight to the mapping using STARsolo. 

![pipeline](/img/Preprocs_Pipeline.png)

---

## Installation & Setup

In order to run the pipeline, you must have Nextflow and Conda installed. Here is the official Nextflow GitHub page and a link to their documentation. The conda installation can be either Anaconda or Miniconda. 

### Conda environments
The default conda environment is called 'bca_int', which is automatically created and activated using the .yml file upon running the pipeline. For a few subprocesses, different conda environments were created as they were conflicting with certain versions of packages. All of them, except for spipe (see below) are installed and activated automatically upon execution using the same approach. 

To manually activate the conda environment:
```
# Create environment
conda env create -n bca_int -f bca_int_env.yaml

# Activate environment
conda activate bca_int
```


### Installing external packages

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
- [ ] Accession list (to download directly from SRA) or a /data/fastq folder with raw FASTQ files.
- [ ] Genome annotation files (FASTA and GFF/GTF files)
- [ ] split-pipe (spipe) environment

### Optional files & parameters
- [ ] spec.yaml - created using seqspec
- [ ] Configuration file

### Running the Pipeline
```
sbatch submit_nextflow.sh main.nf
```
