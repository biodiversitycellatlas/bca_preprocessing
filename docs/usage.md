## Setup

### 1. Create a samplesheet

The samplesheet should be a comma-seperated file, specifying the names and locations of the files and details necessary for pipeline execution. Depending on the chosen sequencing technique the order of the FASTQ files is altered, R1 might contain the cDNA while in other cases this might contain the Cell barcode & UMI's, check the available documentation or manually.

When analysing sci-RNA-seq3 data, it is necessary to also provide the index to the p5, p7 and rt's for each sample to analyse this data successfully. These indexes are defined in the [sci-RNA-seq3 barcode whitelist file](assets/sciRNAseq3_bwl.txt), where you can match the sequences to the corresponding barcode, which is used in the samplesheet.

In the case of Parse Biosciences data, the column of p5 should be filled with the group-well definitions, where:

    Wells are specified in blocks, ranges, or individually like this:
        'A1:C6' specifies a block as [top-left]:[bottom-right]; A1-A6, B1-B6, C1-C6.
        'A1-B6' specifies a range as [start]-[end]; A1-A12, B1-6.
        'C4' specifies a single well.
        Multiple selections are joined by commas (no space), e.g. 'A1-A6,B1:D3,C4'

In the table below, the variables are summarized, with an example samplesheet listed [here](conf/example_samplesheet.csv).

| Variable       | Required/Optional | Description                                                                                                                                                             |
| -------------- | ----------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| sample         | **Required**      | Must be unique unless you want the FASTQ files to be merged after demultiplexing.                                                                                       |
| fastq_cDNA     | **Required**      | Path to the FASTQ file containing cDNA                                                                                                                                  |
| fastq_CB_UMI   | **Required**      | Path to the FASTQ file containing the cell barcode & UMI                                                                                                                |
| expected_cells | **Required**      | Number of expected cells                                                                                                                                                |
| p5             | Optional          | Only required for sci-RNA-seq3 (p5-associated barcode from [sci-RNA-seq3 barcode whitelist file](assets/sciRNAseq3_bwl.txt) & Parse Biosciences (group-well definition) |
| p7             | Optional          | Only required for sci-RNA-seq3 (p7-associated barcode from [sci-RNA-seq3 barcode whitelist file](assets/sciRNAseq3_bwl.txt)                                             |
| rt             | Optional          | Only required for sci-RNA-seq3 (rt-associated barcode from [sci-RNA-seq3 barcode whitelist file](assets/sciRNAseq3_bwl.txt)                                             |


### 2. Edit (or create) a custom configuration file

To set the custom parameters for each run, the easiest solution is to fill in the fiels in [`conf/custom_parameters.config`](conf/custom_parameters.config). This custom configuration file extends the general `nextflow.config` file in the base of this repository. `conf/custom_parameters.config` contains the minimum variables in order to run this pipeline, and are described in the table below in more detail.

To customize the run, you can add other (optional) variables to the `conf/custom_parameters.config` file. If you have a multi-species experiment, one configuration file per species must be created in order to analyze the data with the corresponding genome annotation files.

Within each custom configuration file the following variables can be defined:

| Variable                 | Required/Optional | Description                                                                                                                                                                                                                       |
| ------------------------ | ----------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `input`                  | **Required**      | Path to the samplesheet.                                                                                                                                                                                                          |
| `outdir`                 | **Required**      | Path to the results/output directory; must exist before running.                                                                                                                                                                  |
| `protocol`               | **Required**      | Specifies the sequencing technology used (must be one of the following: `"oak_seq"`, `"10xv3"`, `"parse_biosciences_WT_mini"` or `"parse_biosciences_WT"`, `"bd_rhapsody"`, `"sciRNAseq3"` , `"ultima_genomics"` or `"seqspec"`). |
| `ref_fasta`              | **Required**      | Path to the genome FASTA file used for mapping reads.                                                                                                                                                                             |
| `ref_gtf`                | **Required**      | Path to the GTF/GFF file formatted for STARsolo.                                                                                                                                                                                  |
| `ref_gtf_alt`            | Optional          | Path to the GTF/GFF file formatted specifically for analysis with Parse Biosciences / CellRanger pipeline. Defaults to the same path as `ref_gtf`.                                                                                |
| `seqspec_file`           | Optional          | Path to the seqspec file.                                                                                                                                                                                                         |
| `mt_contig`              | Optional          | Name of the mitochondrial contig in the reference annotation, used to calculate mtDNA content. Default set to `"chrM M MT"`.                                                                                                      |
| `grep_rrna`              | Optional          | String used to grep ribosomal RNA (rRNA) reads from annotations. Default set to `"rRNA"`                                                                                                                                          |
| `mapping_software`       | Optional          | Software used to map reads (must be one of the following: `"starsolo"`, `"alevin"` or `"both"`). Default set to `"starsolo"`.                                                                                                     |
| `perform_geneext`        | Optional          | Boolean flag to enable or disable the gene extension step in preprocessing. Default is `false`.                                                                                                                                   |
| `perform_featurecounts`  | Optional          | Boolean flag to enable or disable calculation of mtDNA & rRNA percentages. Default is `true`.                                                                                                                                     |
| `perform_kraken`         | Optional          | Boolean flag to enable or disable Kraken2 classification of unmapped reads. Default is `false`.                                                                                                                                   |
| `perform_cellbender`     | Optional          | Boolean flag to enable or disable removal of ambient RNA using CellBender. Default is `false`.                                                                                                                                    |
| `kraken_db_path`         | Optional          | Path to the Kraken2 database used for taxonomic classification of unmapped reads, if empty, a default database will be installed.                                                                                                 |
| `perform_cellranger`     | Optional          | Boolean flag to enable or disable the CellRanger pipeline. Default is `false`.                                                                                                                                                    |
| `splitpipe_installation` | Optional          | Path to the split-pipe installation folder, that can be used as a control.                                                                                                                                                        |
| `splitpipe_conda_env`    | Optional          | Path to the split-pipe conda environment created by following these [instructions](docs/INSTALLATION_EXTERNAL_PIPELINES.md), required if running the split-pipe pipeline.                                                         |
| `rhapsody_installation`  | Optional          | Path to the BD-Rhapsody pipeline installation folder, that can be used as a control.                                                                                                                                              |

To modify the behaviour of certain processes, additional variables can be added to the configuration file. An overview of the extended custom parameters is listed [here](docs/CONFIGURATION_PARAMETERS.md).


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

nextflow run -profile <institution_config>,conda -c </path/to/custom_parameters.config> -w </path/to/workdir>

####### OR - submitting pipeline through a bash script #######

# Submit pipeline to SLURM queue
sbatch submit_nextflow.sh main.nf
```
