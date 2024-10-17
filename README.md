# BCA_preprocessing


## TODO

- [ ] Retrieve sequencing technique directly from seqspec file 
- [ ] Loop over given species 
- [ ] Automatically create an accession_file if fastq folder is found 
- [ ] Integrate seqspec into mapping (https://github.com/pachterlab/seqspec/blob/main/docs/UNIFORM.md) 

***

## Installation

### Conda environment
The packages for this project were installed using conda, and installed in the environment called 'bca'. To run the Nextflow pipeline, you can replicate our conda environment by running the following line of code:
```
conda create --name bca --file bca_environment.txt
```
This will create a new conda environment using the mentioned txt file, which contains an (explicit) list of packages, and their download/installation links. 

This environment is automatically activated when running the pipeline, from within the 'submit_nextflow.sh' script.


### Library structure
```
 ┌─ code/ ────┌─
 │            ├─
 │            └─
 │            ┌─ accession_lists/
 ├─ data/ ────┼─ {abbreviation species, e.g.: spol}
 │            └─ ...      
 └─ logs/
```


## Usage

### Prerequisites
- [ ] Accession list (example) when wanting to download data from the SRA or a /data/fastq folder with raw data.
- [ ] spec.yaml - created using seqspec
- [ ] Genome annotation files
- [ ] Cell Barcode file(s)

### Running the Pipeline
```
sbatch submit_nextflow.sh preprocs_pipeline.nf
```
