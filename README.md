# RNA-seq Pipeline for preprocessing single-cell and single-nucleus data

---

## Table of Contents

1. [Overview](#overview)
2. [Tools](#tool-overview)
3. [Requirements](#requirements)
4. [Installation & Setup](#installation--setup)
5. [Usage](#usage)
6. [Example Execution](#example-execution)
7. [Output & Logs](#output--logs)

---

## Overview

### Data Types
- **Single-cell RNA-seq** (scRNA-seq)
- **Single-nucleus RNA-seq** (snRNA-seq)

### Sequencing Technologies
- **Parse Biosciences**
- **BD Rhapsody**

This pipeline reads input samples from an accession list and automates:
1. Preprocessing
2. Quality control (FastQC, MultiQC)
3. Alignment
4. Optional gene extension and reindexing
5. Mapping statistics & saturation metrics
6. Ambient RNA removal

---

## Installation & Setup

### Conda environment
The packages for this project were installed using conda, in the environment called 'bca_int'. To run the Nextflow pipeline, you can replicate our conda environment by running the following line of code:
```
conda create --name bca_int --file bca_int_environment.txt
```

### Installing external packages


#### Parse Biosciences
```
conda create –name spipe 
bash ./install_dependencies_conda.sh -i -y 
pip install . --no-cache-dir 
```

You might get the following error during execution of the previous steps, here's how we handled the following errors:
```
# Pip not installed
conda install anaconda::pip=23.3.1 

# Failed building wheel for louvain
conda install -c conda-forge python-igraph 
pip install cmake 

# AttributeError: module 'numpy' has no attribute 'NAN'. Did you mean: 'nan'? 
# In the file 'utils.py', in line 'def report_percent_str(num, den=1, round_to=2, zero=np.NAN, perchar=True): ' 
# replace NAN by nan. 
After saving, rerun:
pip install . --no-cache-dir 
```

To test if the installation of split-pipe was successful:
```
split-pipe -h 
```


#### GeneExt
```
# clone repository
git clone https://github.com/sebepedroslab/GeneExt.git

# create environment
conda env create -n geneext -f environment.yaml

# activate environment
conda activate geneext
```

### Library structure
```
 ┌─ code/ ────┌─ integrated_pipe/
 │            └─ seperate_sbatchs/
 │
 │            ┌─ accession_lists/
 ├─ data/ ────┼─ experiment/
 │            └─ ..
 │
 ├─ seq_techniques/ ────┌─ bd_rhapsody/
 │                      └─ parse_biosciences/
 │    
 └─ logs/
```


## Usage

### Prerequisites
- [ ] Accession list (example) when wanting to download data from the SRA or a /data/fastq folder with raw data.
- [ ] spec.yaml - created using seqspec
- [ ] Genome annotation files
- [ ] Downloaded required external packages using this README

### Running the Pipeline
```
sbatch submit_nextflow.sh main.nf
```
