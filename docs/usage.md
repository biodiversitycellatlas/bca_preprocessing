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
| fastq_cDNA     | **Required**      | Path to the FASTQ file containing cDNA. Uncompressed (`.fastq`/`.fq`) and gzipped (`.fastq.gz`/`.fq.gz`) files are both accepted, and may be mixed within a sample.      |
| fastq_CB_UMI   | **Required**      | Path to the FASTQ file containing the cell barcode & UMI. Uncompressed and gzipped files are both accepted.                                                              |
| expected_cells | **Required**      | Number of expected cells                                                                                                                                                |
| manual_cutoff  | Optional          | Per-sample UMI threshold used to call cells. Required for every sample when `cellfilter_method = "manual_cutoff"`, ignored otherwise. See [Cell calling](CONFIGURATION_PARAMETERS.md#cell-calling). |
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
| `protocol`               | **Required**      | Specifies the sequencing technology used (must be one of the following: `"oak_seq"`, `"10xv1"`, `"10xv2"`, `"10xv3"`, `"10xv4"`, `"parse_biosciences_WT_mini"` or `"parse_biosciences_WT"`, `"bd_rhapsody_v1"`, `"bd_rhapsody_enhancedbeads"`, `"sciRNAseq3"` , `"ultima_genomics"`, `"marsseq_v1"`, `"marsseq_v2"` or `"seqspec"`). |
| `ref_fasta`              | **Required**      | Path to the genome FASTA file used for mapping reads.                                                                                                                                                                             |
| `ref_gtf`                | **Required**      | Path to the GTF/GFF file formatted for STARsolo.                                                                                                                                                                                  |
| `ref_gtf_alt`            | Optional          | Path to the GTF/GFF file formatted specifically for analysis with Parse Biosciences / CellRanger pipeline. Defaults to the same path as `ref_gtf`.                                                                                |
| `seqspec_file`           | Optional          | Path to the seqspec file.                                                                                                                                                                                                         |
| `mt_contig`              | Optional          | Name(s) of the mitochondrial contig(s) in the reference annotation, used to calculate mtDNA content. Multiple contigs can be given separated by whitespace, in which case reads on any of them count as mitochondrial. Default set to `"chrM M MT"`. |
| `grep_rrna`              | Optional          | String used to grep ribosomal RNA (rRNA) reads from annotations. Default set to `"rRNA"`                                                                                                                                          |
| `mapping_software`       | Optional          | Software used to map reads (must be one of the following: `"starsolo"`, `"alevin"` or `"both"`). Default set to `"starsolo"`.                                                                                                     |
| `cellfilter_method`      | Optional          | How cells are called: `"star_solocellfilter"` (each mapper's own cell filter), `"second_derivative"` or `"manual_cutoff"` (a per-sample threshold from the samplesheet). Default set to `"second_derivative"`, which re-calls cells for both mappers and recomputes every cell-dependent statistic in the report on the retained cells. See [Cell calling](CONFIGURATION_PARAMETERS.md#cell-calling). |
| `run_method`             | Optional          | One of `"standard"`, `"geneext_only"`, `"external_pipeline_only"` or `"post_mapping"`. Default set to `"standard"`. Use `"post_mapping"` to redo everything after mapping on a finished run — for instance to re-call cells with a different `expected_cells` or a manual cutoff — without mapping again. See [Re-running after mapping](CONFIGURATION_PARAMETERS.md#re-running-after-mapping). |
| `previous_outdir`        | Optional          | Only with `run_method = "post_mapping"`: the finished run's results directory to read the mapping results back from. Defaults to `outdir`.                                                                                        |
| `perform_geneext`        | Optional          | Boolean flag to enable or disable the gene extension step in preprocessing. Default is `false`.                                                                                                                                   |
| `perform_featurecounts`  | Optional          | Boolean flag to enable or disable calculation of mtDNA & rRNA percentages. Default is `false`.                                                                                                                                     |
| `perform_kraken`         | Optional          | Boolean flag to enable or disable Kraken2 classification of unmapped reads. Default is `false`.                                                                                                                                   |
| `ambient_rna_remover`    | Optional          | Software used to remove ambient RNA: `"cellsweep"` or `"none"`. Default is `"cellsweep"`.                                                                                                                                         |
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
