# Installing Commercial Software (optional)

---

## Table of Contents

1. [Cell Ranger (10x Genomics & OAK-seq)](#cell-ranger-10x-genomics--oak-seq)
2. [BD Single-Cell Multiomics Software (BD Rhapsody)](#bd-single-cell-multiomics-software-bd-rhapsody)
3. [split-pipe (Parse Biosciences)](#split-pipe-parse-biosciences)
4. [sci-rocket (sci-RNA-seq3)](#sci-rocket-sci-rna-seq3)

---

### Cell Ranger (10x Genomics & OAK-seq)
Follow the installation guide on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-in), and downloaded Cell Ranger using wget for the latest version on the [downloads](https://www.10xgenomics.com/support/software/cell-ranger/downloads) page.

---

### BD Single-Cell Multiomics Software (BD Rhapsody)

First, download the ['cwl' folder](https://bitbucket.org/CRSwDev/cwl/src/master/) of the BD Single-Cell Multiomics Software, hosted on Bitbucket. This repository contains the CWL and YML files required to run the BD Rhapsody pipeline locally. Secondly, pull the docker image using Docker or Apptainer (see code for Apptainer example). 

```
# Pull docker image using Apptainer
apptainer pull docker://bdgenomics/rhapsody
```

---

### split-pipe (Parse Biosciences)
In order to incorporate the commercial Parse Biosciences pipeline (also called split-pipe), first create an account and install the code from the personal account page on the Parse Biosciences website. For this project, we tested version 1.3.1, 

Instructions on creating the split-pipe conda environment:

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
conda install anaconda::pip 

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


### sci-rocket (sci-RNA-seq3)

Github page: https://github.com/odomlab2/sci-rocket/tree/main
Documentation: https://odomlab2.github.io/sci-rocket/

Version(s) tested: 
```
# Installation
git clone https://github.com/odomlab2/sci-rocket

## created conda environment with snakemake
## created config & samplesheet

# Run sci-rocket
cd sci-rocket/workflow/
snakemake --use-conda --configfile <config file> --cores 1
```

Using the software out-of-the-box, we encountered two issues:
- Error: clock skew problem
- Error: creating STAR index



#### 1. Error: clock skew problem

Shortened error message: 
```
Output read_1.fastq.gz has older modification time than input read_2.fastq.gz.
This could indicate a clock skew problem in your network and would trigger a rerun of this job in the next execution and should therefore be fixed on system level.
```

Which could be fixed by modifying the file sci-rocket/workflow/rules/step1_bcl2fastq.smk:
```
rule merge_sequencing_runs:
    shell:
        """
-        # If only one sequencing run, then just hardlink it (and remove the original).
-        if [ {params.total_sequencing_runs} -eq 1 ]; then
-            ln {input.R1} {output.R1}
-            ln {input.R2} {output.R2}
-       else
            cat {input.R1} > {output.R1}
            cat {input.R2} > {output.R2}
-       fi
        """
```



#### 2. Error: creating STAR index

In the documentation it mentions that in the config file, if you set the variable star_index: False, automatically a STAR index will be generated. This seems to be the true, but does raise an error when passing this newly created STAR index onto the next process. So in general, it is better to create a STAR index yourself beforehand, and set the path manually.

Shortened error message: 
```
MissingOutputException in rule generate_index_STAR in file "sci-rocket/workflow/rules/step3_alignment.smk", line 42:
Job 86  completed successfully, but some output files are missing. Missing files after 5 seconds. This might be due to filesystem latency. If that is the case, consider to increase the wait time with --latency-wait
```


