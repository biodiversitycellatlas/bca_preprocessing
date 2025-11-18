# BCA Pre-Processing Pipeline

---

## Table of Contents

1. [Overview](#overview)
2. [Installation](#installation)
3. [Setup](#setup)
4. [Usage](#usage)
5. [Output](#output)
6. [Citations](#citations)

---

## Overview

This nextflow pipeline is designed to pre-process single-cell and single-nucleus RNA-seq data. It accepts FASTQ files from multiple sequencing platforms, at the moment being:

- Parse Biosciences
- BD Rhapsody
- 10x Genomics
- OAK-seq
- Ultima Genomics
- sci-RNA-seq3
- and others, when providing a [seqspec](https://github.com/pachterlab/seqspec) file

Depending on the chosen sequencing technique, it handles the processing of the FASTQ files accordingly. Whenever possible, we compared our results to a commercial pre-processing pipeline for that sequencing technique. For example, comparing our Parse Biosciences results to the official split-pipe pipeline from Parse Biosciences. While we cannot provide this commercial software directly, you can install it yourself (e.g. by following [these instructions](docs/INSTALLATION_EXTERNAL_PIPELINES.md)), and provide a path where the installation is located in the configuration file. This way, it will be executed alongside of the BCA pre-processing pipeline.

The pipeline will produce the following output files:

- Raw & Filtered count matrices (exonic, exonic & intronic) from Mapping step​
- Filtered count matrices (h5 files) ​from (optional) Filtering step​
- MultiQC report

![pipeline](/img/Preprocs_Pipeline.png)

---

## Installation

1. **Clone the bca_preprocessing repository:**

```
git clone https://github.com/biodiversitycellatlas/bca_preprocessing.git
```

2. **Clone submodules**

The pipeline uses two git submodules, [10x_saturate](https://github.com/zolotarovgl/10x_saturate) to plot the saturation curve and [GeneExt](https://github.com/zolotarovgl/GeneExt) for an extended gene annotation file. These are not included automatically when performing git clone, so have to be updated explicitly by:

```
git submodule init
git submodule update
```

3. **Conda & Nextflow**

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

4. **(Optional) Installing external pipelines as validation**

After following these [installation instructions](docs/INSTALLATION_EXTERNAL_PIPELINES.md) for some of the sequencing technologies, see table below, users can run external pipelines simultaneaously with the BCA pre-processing pipeline. Depending on the sequencing technique, you only have to provide the path to the installation as within the [`conf/custom_parameters.config`](conf/custom_parameters.config) file or a boolean flag to enable/disable the execution, see Setup explenation below, and it will automatically start.

This option is limited to the following sequencing technologies:
| Sequencing technology | External pipeline | Requires manual installation |
|-----------------------|-------------------|-------------------|
| 10x Genomics | [10x Genomics Cell Ranger](https://github.com/10XGenomics/cellranger) | :x: |
| OAKseq | [10x Genomics Cell Ranger](https://github.com/10XGenomics/cellranger) | :x: |
| Ultima Genomics | [10x Genomics Cell Ranger](https://github.com/10XGenomics/cellranger) | :x: |
| ScaleBio | [ScaleRna](https://github.com/ScaleBio/ScaleRna) | :x: |
| Parse Biosciences | [split-pipe](https://support.parsebiosciences.com/hc/en-us/articles/27066395947412-How-Do-I-Analyze-my-Parse-Biosciences-Data) | :heavy_check_mark: |
| BD-Rhapsody | [BD Rhapsody™ Sequence Analysis Pipeline](https://www.bdbiosciences.com/en-us/products/software/rhapsody-sequence-analysis-pipeline) | :heavy_check_mark: |

---

## Setup

### 1. Create a samplesheet

The samplesheet should be a comma-seperated file, specifying the names and locations of the files and details necessary for pipeline execution. Depending on the chosen sequencing technique the order of the FASTQ files is altered, R1 might contain the cDNA while in other cases this might contain the Cell barcode & UMI's, check the available documentation or do a manual inspection.

When analysing sci-RNA-seq3 data, it is necessary to also provide the index to the p5, p7 and rt's for each sample to analyse this data successfully. These indexes are defined in the [sci-RNA-seq3 barcode whitelist file](assets/sciRNAseq3_bwl.txt), where you can match the sequences to the corresponding barcode, which is used in the samplesheet. It is also important to note that for sci-RNA-seq3 data, the p5 + p7 should be present in the headers of the FASTQ files, see [this bcl2fastq script](assets/bcl2fastq2.sh).

In the case of Parse Biosciences data, the rt column should be filled with the group-well definitions, where:

> Wells are specified in blocks, ranges, or individually like this:<br/> >&nbsp;&nbsp;&nbsp;&nbsp;'A1:C6' specifies a block as [top-left]:[bottom-right]; A1-A6, B1-B6, C1-C6.<br/> >&nbsp;&nbsp;&nbsp;&nbsp;'A1-B6' specifies a range as [start]-[end]; A1-A12, B1-6.<br/> >&nbsp;&nbsp;&nbsp;&nbsp;'C4' specifies a single well.<br/>
> Multiple selections are joined by commas (no space), e.g. 'A1-A6,B1:D3,C4'

In the table below, the available variables are summarized:
| Variable | Required/Optional | Description |
| -------------- | ----------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| sample | **Required** | Must be unique unless you want the FASTQ files to be merged after demultiplexing. |
| fastq_cDNA | **Required** | Path to the FASTQ file containing cDNA |
| fastq_CB_UMI | **Required** | Path to the FASTQ file containing the cell barcode & UMI |
| fastq_indices | Optional | Path to the FASTQ index file(s), to provide both I1 and I2, use an asterisk to the path like /path/name_I\* |
| expected_cells | **Required** | Number of expected cells |
| p5 | Optional | Only required for sci-RNA-seq3 (p5-associated barcode from [sci-RNA-seq3 barcode whitelist file](assets/sciRNAseq3_bwl.txt)|
| p7 | Optional | Only required for sci-RNA-seq3 (p7-associated barcode from [sci-RNA-seq3 barcode whitelist file](assets/sciRNAseq3_bwl.txt) |
| rt | Optional | Only required for sci-RNA-seq3 (rt-associated barcode from [sci-RNA-seq3 barcode whitelist file](assets/sciRNAseq3_bwl.txt) & Parse Biosciences (group-well definition) |

To illustrate how the samplesheet would be filled across the different sequencing techniques, the table below is given, as well as an example samplesheet published [here](conf/example_samplesheet.csv). The names within the parenthesis of the p5 and rt column indicate the official names, often referenced like this in the official documentation of this protocol.

| sample                    | fastq_cDNA | fastq_BC_UMI | fastq_indices | expected_cells | p5             | p7  | rt            |
| ------------------------- | ---------- | ------------ | ------------- | -------------- | -------------- | --- | ------------- |
| sciRNAseq3_example        | R2         | R1           |               | expected_cells | p5             | p7  | rt            |
| bd_rhapsody_example       | R2         | R1           |               | expected_cells |                |     |               |
| parse_biosciences_example | R1         | R2           |               | expected_cells |                |     | rt (wells)    |
| 10xv3_example             | R2         | R1           |               | expected_cells |                |     |               |
| oak_seq_example           | R2         | R1           |               | expected_cells |                |     |               |
| ultima_genomics_example   | R2         | R1           |               | expected_cells |                |     |               |
| scalebio_example          | R2         | R1           | I1 & I2       | expected_cells | p5 (libIndex2) |     | p7 (barcodes) |
| scalebio_quantum_example  | R1         | R2           | I1 & I2       | expected_cells | p5 (libIndex2) |     | p7 (barcodes) |

### 2. Edit (or create) a custom configuration file

To set the custom parameters for each run, the easiest solution is to fill in the fiels in [`conf/custom_parameters.config`](conf/custom_parameters.config). This custom configuration file extends the general [`nextflow.config`](nextflow.config) file in the base of this repository. [`conf/custom_parameters.config`](conf/custom_parameters.config) contains the minimum variables in order to run this pipeline, and are described in the table below in more detail.

To customize the run, you can add other (optional) variables to the [`conf/custom_parameters.config`](conf/custom_parameters.config) file. If you have a multi-species experiment, one configuration file per species must be created in order to analyze the data with the corresponding genome annotation files.

Within each custom configuration file the following variables can be defined:

| Variable                | Required/Optional | Description                                                                                                                                                                                                                       |
| ----------------------- | ----------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `input`                 | **Required**      | Path to the samplesheet.                                                                                                                                                                                                          |
| `outdir`                | **Required**      | Path to the results/output directory; must exist before running.                                                                                                                                                                  |
| `protocol`              | **Required**      | Specifies the sequencing technology used (must be one of the following: `"oak_seq"`, `"10xv3"`, `"parse_biosciences_WT_mini"` or `"parse_biosciences_WT"`, `"bd_rhapsody"`, `"sciRNAseq3"` , `"ultima_genomics"` or `"seqspec"`). |
| `bc_whitelist`          | **Required**      | Path or link to the barcode whitelist file(s). If it's a link, it will be automatically downloaded and unzipped if applicable.                                                                                                    |
| `ref_fasta`             | **Required**      | Path to the genome FASTA file used for mapping reads.                                                                                                                                                                             |
| `ref_gtf`               | **Required**      | Path to the GTF/GFF file formatted for STARsolo.                                                                                                                                                                                  |
| `run_method`            | Optional          | Method of running the pre-processing pipeline, demonstrated in the [pipeline diagram](img/Preprocs_Pipeline.png), currently either `"standard"`, `"geneext_only"` or `"exteral_pipeline_only"`. Default is set to `"standard"`.   |
| `seqspec_file`          | Optional          | Path to the seqspec file.                                                                                                                                                                                                         |
| `mapping_software`      | Optional          | Software used to map reads (must be one of the following: `"starsolo"`, `"alevin"` or `"both"`). Default set to `"starsolo"`.                                                                                                     |
| `perform_geneext`       | Optional          | Boolean flag to enable or disable the gene extension step in preprocessing. Default is `false`.                                                                                                                                   |
| `perform_featurecounts` | Optional          | Boolean flag to enable or disable calculation of mtDNA & rRNA percentages. Default is `false`.                                                                                                                                    |
| `perform_kraken`        | Optional          | Boolean flag to enable or disable Kraken2 classification of unmapped reads. Default is `false`.                                                                                                                                   |
| `perform_cellbender`    | Optional          | Boolean flag to enable or disable removal of ambient RNA using CellBender. Default is `false`.                                                                                                                                    |
| `kraken_db_path`        | Optional          | Path to the Kraken2 database used for taxonomic classification of unmapped reads, if empty, a default database will be installed.                                                                                                 |

To modify the behaviour of certain processes or enable external pipelines, additional variables can be added to the configuration file. An overview of the extended custom parameters is listed [here](docs/CONFIGURATION_PARAMETERS.md).

---

## Usage

### Pre-requisites:

- [ ] Created a samplesheet in CSV format (see [conf/example_samplesheet.csv](conf/example_samplesheet.csv))
- [ ] Edited ([`conf/custom_parameters.config`](conf/custom_parameters.config)) or created a new custom configuration file
- [ ] Conda & Nextflow available in base environment

### Running the Pipeline

> [!WARNING]
> The pipeline must be run in the conda base environment, it cannot activate the different environments properly with a prior environment activated. It should have access to run both `nextflow` and `conda` in the commandline.

```

nextflow run
    -profile <institution_config>,conda
    -c </path/to/custom_parameters.config>
    -w </path/to/workdir>

# OR - submitting pipeline through a bash script
sbatch submit_nextflow.sh main.nf
```

---

## Output

The pipeline will produce the following output files:

- Raw & Filtered count matrices (exonic, exonic & intronic) from Mapping step​
- Filtered count matrices (h5 files) ​from (optional) Filtering step​
- MultiQC report

---

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [CITATIONS.md](CITATIONS.md) file.
