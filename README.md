# BCA_preprocessing

## Installation

### Conda environment
The packages for this project were installed using conda, and installed in the environment called 'bca_int'. To run the Nextflow pipeline, you can replicate our conda environment by running the following line of code:
```
conda create --name bca_int --file bca_int_environment.txt
```
This will create a new conda environment using the mentioned txt file, which contains an (explicit) list of packages, and their download/installation links. 

This environment is automatically activated when running the pipeline, from within the 'submit_nextflow.sh' script.

### Installing external packages

#### Parse Biosciences
```
conda create –name spipe 
bash ./install_dependencies_conda.sh -i -y 
pip install . --no-cache-dir 
```

You might get the following error during execution of the previous steps, here's how we handled the following errors:

- Pip not installed:
```
conda install anaconda::pip=23.3.1 
```

-     Failed building wheel for louvain
```
conda install -c conda-forge python-igraph 
pip install cmake 
```

- AttributeError: module 'numpy' has no attribute 'NAN'. Did you mean: 'nan'? 
In the file 'utils.py', in line 'def report_percent_str(num, den=1, round_to=2, zero=np.NAN, perchar=True): ' replace NAN by nan. After saving, rerun:
```
pip install . --no-cache-dir 
```

To test if the installation of split-pipe was successful:
```
split-pipe -h 
```

#### 10x_saturation
```
```

#### GeneExt
```
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
sbatch submit_nextflow.sh preprocs_pipeline.nf
```
