# BCA Pre-Processing Pipeline

[![website][]][website-link]
[![nextflow][]][nextflow-link]
[![conda][]][conda-link]
[![docker][]][docker-link]
[![singularity][]][singularity-link]

[website]: https://img.shields.io/badge/website-biodiversitycellatlas.org-blue
[website-link]: https://biodiversitycellatlas.org
[nextflow]: https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io
[nextflow-link]: https://www.nextflow.io/
[conda]: http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda
[conda-link]: https://docs.conda.io/en/latest/
[docker]: https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker
[docker-link]: https://www.docker.com/
[singularity]: https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000
[singularity-link]: https://sylabs.io/docs/

## Table of Contents

1. [Overview](#overview)
2. [Installation](#installation)
3. [Setup](#setup)
4. [Usage](#usage)
5. [Output](#output)
6. [Citations](#citations)

## Overview

This nextflow pipeline is designed to pre-process single-cell and single-nucleus RNA-seq data. It accepts FASTQ files from multiple sequencing platforms, at the moment being:

- Parse Biosciences
- BD Rhapsody
- 10x Genomics
- OAK-seq
- Ultima Genomics
- sci-RNA-seq3
- MARS-seq (v1 and v2)
- and others, when providing a [seqspec](https://github.com/pachterlab/seqspec) file

Depending on the chosen sequencing technique, it handles the processing of the FASTQ files accordingly. Whenever possible, we compared our results to a commercial pre-processing pipeline for that sequencing technique. For example, comparing our Parse Biosciences results to the official split-pipe pipeline from Parse Biosciences. While we cannot provide this commercial software directly, you can install it yourself (e.g. by following [these instructions](docs/INSTALLATION_EXTERNAL_PIPELINES.md)), and provide a path where the installation is located in the configuration file. This way, it will be executed alongside of the BCA pre-processing pipeline.

The pipeline will produce the following output files:

- Quality control reports for raw FASTQ files
- Raw & Filtered count matrices (exonic, exonic & intronic) from STARsolo and/or Alevin-fry
- (Optional) mtDNA and rRNA statistics
- (Optional) Sequencing saturation analysis
- (Optional) Extended gene annotation file
- (Optional) Filtered count matrices (h5 files) ​from CellBender step
- (Optional) Kraken2 taxonomic classification
- MultiQC report
- Portable HTML dashboard displaying results

![pipeline](/img/Preprocs_Pipeline.png)

---

## Installation

1. **Download the latest release**

Go to the [Releases](https://github.com/biodiversitycellatlas/bca_preprocessing/releases) page and download the `.zip` or `.tar.gz` file for the version you want.

2. **Unpack the files**

```
unzip bca_preprocessing-<version>.zip
cd bca_preprocessing-<version>
```

or for tar.gz:

```
tar -xzf bca_preprocessing-<version>.tar.gz
cd bca_preprocessing-<version>
```

3. **Download submodules**

The pipeline uses two submodules, [10x_saturate](https://github.com/zolotarovgl/10x_saturate) to plot the saturation curve and [GeneExt](https://github.com/zolotarovgl/GeneExt) for an extended gene annotation file. These are not included automatically, so have to be installed explicitly by running:

```
bash fetch_submodules.sh
```

4. **Conda & Nextflow**

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

5. **(Optional) Installing external pipelines as validation**

After following these [installation instructions](docs/INSTALLATION_EXTERNAL_PIPELINES.md) for some of the sequencing technologies, see table below, users can run external pipelines simultaneaously with the BCA pre-processing pipeline. Depending on the sequencing technique, you only have to provide the path to the installation as within the [`conf/custom_parameters.config`](conf/custom_parameters.config) file or a boolean flag to enable/disable the execution, see Setup explenation below, and it will automatically start.

This option is limited to the following sequencing technologies:
| Sequencing technology | External pipeline | Requires manual installation |
|-----------------------|-------------------|-------------------|
| 10x Genomics | [10x Genomics Cell Ranger](https://github.com/10XGenomics/cellranger) | :x: |
| OAK-seq | [10x Genomics Cell Ranger](https://github.com/10XGenomics/cellranger) | :x: |
| Ultima Genomics | [10x Genomics Cell Ranger](https://github.com/10XGenomics/cellranger) | :x: |
| Parse Biosciences | [split-pipe](https://support.parsebiosciences.com/hc/en-us/articles/27066395947412-How-Do-I-Analyze-my-Parse-Biosciences-Data) | :heavy_check_mark: |
| BD-Rhapsody | [BD Rhapsody™ Sequence Analysis Pipeline](https://www.bdbiosciences.com/en-us/products/software/rhapsody-sequence-analysis-pipeline) | :heavy_check_mark: |

## Setup

### 1. Create a samplesheet

The samplesheet should be a comma-seperated file, specifying the names and locations of the files and details necessary for pipeline execution. Depending on the chosen sequencing technique the order of the FASTQ files is altered, R1 might contain the cDNA while in other cases this might contain the Cell barcode & UMI's, check the available documentation or do a manual inspection.

Some protocols need more than the FASTQ paths: sci-RNA-seq3 needs the `p5`, `p7` and `rt` columns filled in, Parse Biosciences needs the group-well definition in `rt`, and Parse Biosciences and MARS-seq put the cDNA in read 1 rather than read 2. These per-protocol requirements, together with the read layouts and the pre-processing steps they trigger, are described in [Protocol-specific steps & whitelists](docs/PROTOCOLS_AND_WHITELISTS.md).

In the table below, the available variables are summarized:
| Variable | Required/Optional | Description |
| -------------- | ----------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| sample | **Required** | Must be unique unless you want the FASTQ files to be merged after demultiplexing. |
| fastq_cDNA | **Required** | Path to the FASTQ file containing cDNA. Uncompressed (`.fastq`/`.fq`) and gzipped (`.fastq.gz`/`.fq.gz`) files are both accepted, and may be mixed within a sample. |
| fastq_CB_UMI | **Required** | Path to the FASTQ file containing the cell barcode & UMI. Uncompressed and gzipped files are both accepted. |
| fastq_indices | Optional | Path to the FASTQ index file(s), to provide both I1 and I2, use an asterisk to the path like /path/name_I\*. Uncompressed and gzipped files are both accepted. |
| expected_cells | **Required** | Number of expected cells |
| p5 | Optional | Only required for sci-RNA-seq3 |
| p7 | Optional | Only required for sci-RNA-seq3 |
| rt | Optional | Only required for sci-RNA-seq3 & Parse Biosciences (group-well definition) |

To illustrate how the samplesheet would be filled across the different sequencing techniques, the table below is given, as well as an example samplesheet published [here](conf/example_samplesheet.csv). The names within the parenthesis of the p5 and rt column indicate the official names, often referenced like this in the official documentation of this protocol.

| sample                    | fastq_cDNA | fastq_BC_UMI | fastq_indices | expected_cells | p5  | p7  | rt         |
| ------------------------- | ---------- | ------------ | ------------- | -------------- | --- | --- | ---------- |
| sciRNAseq3_example        | R2         | R1           |               | expected_cells | p5  | p7  | rt         |
| bd_rhapsody_example       | R2         | R1           |               | expected_cells |     |     |            |
| parse_biosciences_example | R1         | R2           |               | expected_cells |     |     | rt (wells) |
| 10xv3_example             | R2         | R1           |               | expected_cells |     |     |            |
| oak_seq_example           | R2         | R1           |               | expected_cells |     |     |            |
| ultima_genomics_example   | R2         | R1           |               | expected_cells |     |     |            |
| marsseq_example           | R1         | R2           |               | expected_cells |     |     |            |

### 2. Edit (or create) a custom configuration file

To set the custom parameters for each run, the easiest solution is to fill in the fiels in [`conf/custom_parameters.config`](conf/custom_parameters.config). This custom configuration file extends the general [`nextflow.config`](nextflow.config) file in the base of this repository. [`conf/custom_parameters.config`](conf/custom_parameters.config) contains the minimum variables in order to run this pipeline, and are described in the table below in more detail.

To customize the run, you can add other (optional) variables to the [`conf/custom_parameters.config`](conf/custom_parameters.config) file. If you have a multi-species experiment, one configuration file per species must be created in order to analyze the data with the corresponding genome annotation files.

Within each custom configuration file the following variables can be defined:

| Variable                 | Required/Optional | Description                                                                                                                                                                                                                                                        |
| ------------------------ | ----------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `input`                  | **Required**      | Path to the samplesheet.                                                                                                                                                                                                                                           |
| `outdir`                 | **Required**      | Path to the results/output directory; must exist before running.                                                                                                                                                                                                   |
| `protocol`               | **Required**      | Specifies the sequencing technology used (must be one of the following: `"oak_seq"`, `"10xv1"`, `"10xv2"`, `"10xv3"`, `"10xv4"`, `"parse_biosciences_WT_mini"` or `"parse_biosciences_WT"`, `"bd_rhapsody_v1"`, `"bd_rhapsody_enhancedbeads"`, `"sciRNAseq3"` , `"ultima_genomics"`, `"marsseq_v1"`, `"marsseq_v2"` or `"seqspec"`). |
| `bc_whitelist`           | **Required**      | Path or link to the barcode whitelist file(s), multiple ones separated by whitespace. If links are given, they are automatically downloaded (and unzipped if applicable) for any protocol. Not used by `"marsseq_v1"`/`"marsseq_v2"`, which run without a whitelist. |
| `ref_fasta`              | **Required**      | Path to the genome FASTA file used for mapping reads.                                                                                                                                                                                                              |
| `ref_gtf`                | **Required**      | Path to the GTF/GFF file formatted for STARsolo.                                                                                                                                                                                                                   |
| `run_method`             | Optional          | Method of running the pre-processing pipeline, demonstrated in the [pipeline diagram](img/Preprocs_Pipeline.png), one of `"standard"`, `"geneext_only"`, `"external_pipeline_only"` or `"post_mapping"`. Default is set to `"standard"`. `"post_mapping"` redoes everything after mapping on a finished run; re-calling cells with a different `expected_cells` or a manual cutoff, without mapping again; see [Re-running after mapping](docs/CONFIGURATION_PARAMETERS.md#re-running-after-mapping). |
| `previous_outdir`        | Optional          | Only with `run_method = "post_mapping"`: the finished run's results directory to read the mapping results back from. Defaults to `outdir`.                                                                                                                        |
| `perform_demultiplexing` | Optional          | Boolean flag to enable or disable demultiplexing of the FASTQ files, where applicable. Default is `true`.                                                                                                                                                          |
| `seqspec_file`           | Optional          | Path to the seqspec file.                                                                                                                                                                                                                                          |
| `mapping_software`       | Optional          | Software used to map reads (must be one of the following: `"starsolo"`, `"alevin"` or `"both"`). Default set to `"starsolo"`.                                                                                                                                      |
| `perform_geneext`        | Optional          | Boolean flag to enable or disable the gene extension step in preprocessing. Default is `false`.                                                                                                                                                                    |
| `perform_featurecounts`  | Optional          | Boolean flag to enable or disable calculation of mtDNA & rRNA percentages. Default is `false`.                                                                                                                                                                     |
| `perform_kraken`         | Optional          | Boolean flag to enable or disable Kraken2 classification of unmapped reads. Default is `false`.                                                                                                                                                                    |

To modify the behaviour of certain processes or enable external pipelines, additional variables can be added to the configuration file. An overview of the extended custom parameters is listed [here](docs/CONFIGURATION_PARAMETERS.md).

## Usage

### Pre-requisites:

- [ ] Created a samplesheet in CSV format (see [conf/example_samplesheet.csv](conf/example_samplesheet.csv))
- [ ] Edited ([`conf/custom_parameters.config`](conf/custom_parameters.config)) or created a new custom configuration file
- [ ] Conda & Nextflow available in base environment


### Verifying your setup

Before committing to a full run, you can check that the software the pipeline
depends on is actually available on your systema:

```
# Run all available checks
bash tests/run_tests.sh
```

See [`tests/README.md`](tests/README.md) for the full set of options.


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

## Output

The pipeline outputs are organized into sub-directories based on the analysis steps performed.
Below is an overview of the possible directory structure:

```
output_directory/
├── demultiplex/            # Demultiplexed FASTQ files
├── pipeline_info/          # Run execution details (configs, samplesheet, execution trace/report/timeline)
├── fastqc/                 # Quality control reports for raw FASTQ files
├── mapping_STARsolo/       # STARsolo count matrices and stats
├── mapping_alevin/         # (Optional) Alevin-fry count matrices and AlevinQC
├── summary_results/        # MultiQC report, Per-Cell metrics, and mapping_stats.tsv table
├── dashboard.html          # Portable HTML dashboard displaying results
│
├── cellbender/             # (Optional) CellBender filtered matrices
├── saturation/             # (Optional) Sequencing saturation analysis
├── gene_ext/               # (Optional) Extended GTF file and outputs from GeneExt
├── rRNA_mtDNA/             # (Optional) mtDNA and rRNA results from FeatureCounts
├── kraken/                 # (Optional) Kraken2 taxonomic classification of unmapped reads
│
├── CellRanger_pipeline/    # (Optional) External Cell Ranger outputs
├── ParseBio_pipeline/      # (Optional) External Split-pipe outputs
└── BDrhapsody_pipeline/    # (Optional) External BD Rhapsody outputs
```

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [CITATIONS.md](CITATIONS.md) file.

## Contact us

[<img src="img/LOGOs-CRG-ENG_2014_transparent_back.png" width="130" target="_blank" alt="Centre for Genomic Regulation (CRG)"/>][CRG] [<img src="img/EMBL_EBI_Logo_black.svg" width="200" target="_blank" alt="European Bioinformatics Institute (EMBL-EBI)"/>][EBI] [<img src="img/Wellcome_Sanger_Institute_Logo_Landscape_Digital_RGB_Full_Colour.svg" width="200" target="_blank" alt="Wellcome Sanger Institute (Sanger)"/>][Sanger] [<img src="img/Gordon_and_Betty_Moore_Foundation_logo.svg" width="160" target="_blank" alt="Gordon and Betty Moore Foundation"/>][Moore]

[CRG]: https://crg.eu
[EBI]: https://ebi.ac.uk/
[Sanger]: https://sanger.ac.uk/
[Moore]: https://moore.org
