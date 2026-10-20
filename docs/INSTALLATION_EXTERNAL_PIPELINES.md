# Installing external software (optional)

For some protocols, see table below, users can run external pipelines simultaneaously with the BCA pre-processing pipeline. Depending on the sequencing technique, you only have to provide the path to the installation as within the [`conf/custom_parameters.config`](conf/custom_parameters.config) file or a boolean flag to enable/disable the execution.

This option is limited to the following sequencing technologies:
| Sequencing technology | External pipeline | Requires manual installation |
|-----------------------|-------------------|-------------------|
| 10x Genomics | [Cell Ranger](https://github.com/10XGenomics/cellranger) | :x: |
| OAK | [Cell Ranger](https://github.com/10XGenomics/cellranger) | :x: |
| sci-RNA-seq3 | [sci-rocket](https://github.com/odomlab2/sci-rocket/tree/main) | :heavy_check_mark: |
| Parse Biosciences | [split-pipe](https://support.parsebiosciences.com/hc/en-us/articles/27066395947412-How-Do-I-Analyze-my-Parse-Biosciences-Data) | :heavy_check_mark: |
| BD-Rhapsody | [BD Rhapsody™ Sequence Analysis Pipeline](https://www.bdbiosciences.com/en-us/products/software/rhapsody-sequence-analysis-pipeline) | :heavy_check_mark: |

---


## sci-rocket (sci-RNA-seq3)

[Github page](https://github.com/odomlab2/sci-rocket/tree/main) & [Documentation](https://odomlab2.github.io/sci-rocket/)


```
# Installation
git clone https://github.com/odomlab2/sci-rocket
```

## BD Single-Cell Multiomics Software (BD Rhapsody)

First, download the ['cwl' folder](https://bitbucket.org/CRSwDev/cwl/src/master/) of the BD Single-Cell Multiomics Software, hosted on Bitbucket. This repository contains the CWL and YML files required to run the BD Rhapsody pipeline locally. 

Second, pull the docker image using Docker or Apptainer:

```
# Pull docker image using Apptainer
apptainer pull docker://bdgenomics/rhapsody
```


## split-pipe (Parse Biosciences)
In order to incorporate the commercial Parse Biosciences pipeline (also called split-pipe), first create an account and install the code from the personal account page on the official Parse Biosciences website. 

*Important to note, we only tested version 1.3.1.*

Instructions on creating the split-pipe conda environment:

```
# Create conda environment from file
conda env create -p <destination_path>/<env_name> -f assets/split-pipe_conda_env.yml
conda activate <destination_path>/<env_name>

# Move to your split-pipe installation directory
cd /<splitpipe_dir>/

# Install dependencies
pip install . --no-cache-dir

## declare this new environment in your custom_config file as 'splitpipe_conda_env = <destination_path>/<env_name>'
```

To test if the installation of split-pipe was successful:
```
split-pipe -h
```

