# Usage

### Table of Contents

1. [Before you start](#before-you-start)
2. [Running the pipeline](#running-the-pipeline)
    - [Profiles](#profiles)
    - [Useful options](#useful-options)
    - [Testing your setup](#testing-your-setup)
3. [Run methods](#run-methods)
    - [Full run (`standard`)](#full-run-standard)
    - [Gene extension only (`geneext_only`)](#gene-extension-only-geneext_only)
    - [External pipelines only (`external_pipeline_only`)](#external-pipelines-only-external_pipeline_only)
    - [Re-running after mapping (`post_mapping`)](#re-running-after-mapping-post_mapping)
4. [Filtering UMI matrices (two passes)](#filtering-umi-matrices-two-passes)

---


## Before you start

- [ ] A samplesheet in CSV format, see [`example_samplesheet.csv`](../examples/example_samplesheet.csv) and the [samplesheet configuration per protocol](PROTOCOLS_AND_WHITELISTS.md).
- [ ] A custom configuration file, see [`custom_parameters.config`](../conf/custom_parameters.config) and [all configuration parameters](CONFIGURATION_PARAMETERS.md).
- [ ] Nextflow and one of: Conda, Mamba, Micromamba, Docker, Podman, Singularity or Apptainer.

> [!WARNING]
> When using conda, mamba or micromamba: the pipeline must be run in the base environment, it cannot activate the different environments properly with a prior environment activated.


## Running the pipeline

```
# Submitting the pipeline in a local environment or interactive session on a HPC
nextflow run main.nf
    -profile <institution_config>,<environment_manager>         # e.g. -profile crg,conda
    -c </path/to/custom_parameters.config>                      # e.g. -c conf/custom_parameters.config
    -w </path/to/workdir>                                       # Temporary nextflow files are stored here

# Submitting pipeline to an HPC using SLURM
sbatch submit_nextflow.sh main.nf
```

[`submit_nextflow.sh`](../submit_nextflow.sh) runs `-profile crg,conda -c conf/custom_parameters.config`; edit that line for your own setup. Any extra arguments are passed on to Nextflow, e.g. `sbatch submit_nextflow.sh main.nf -resume`.

Every parameter of the configuration file can also be given on the command line with a double dash (`--run_method post_mapping`), which takes precedence over the file.


### Profiles

Combine one institution profile with one environment manager, separated by a comma.

| Profile | Description |
|---------|-------------|
| `<institution_config>` | Executor and queue settings of your cluster, loaded from [nf-core/configs](https://github.com/nf-core/configs) (e.g. `crg`). Leave it out to run locally. |
| `conda`, `mamba`, `micromamba` | Builds one environment per module from its `environment.yml`. |
| `apptainer`, `singularity` | Pulls a pre-built [Seqera container](https://seqera.io/containers/) (Singularity image) per module on first use. |
| `docker`, `podman` | Pulls the Docker image of the same environment per module. Add `arm` on Apple Silicon to run the amd64 images. |
| `wave` | Builds each module's container from its conda environment with [Wave](https://seqera.io/wave/); use together with a container profile. |
| `test` | Caps every task at 2 CPUs, 4 GB and 1 h, and turns off retries ([`conf/test.config`](../conf/test.config)). |


### Useful options

| Option | Description |
|--------|-------------|
| `-resume` | Re-uses the cached results of a previous run in the same work directory, only re-running what changed. |
| `-preview` | Checks that every process can be launched with the current configuration, without running anything. |
| `--help` | Lists the pipeline parameters (`--help_full` for all of them). |


### Testing your setup

Small subsampled datasets with a matching configuration are included in [`examples/`](../examples/README.md). Run one with the `test` profile before starting on your own data:

```
nextflow run main.nf -profile <institution_config>,conda,test -c examples/10x/conf/spis_10x_parameters.config
```

The environments and container images can also be checked without any data, using the [test suite](../tests/README.md): `bash tests/run_tests.sh`.


## Run methods

The `run_method` parameter decides which part of the pipeline runs.

| `run_method` | Maps reads | Stops after | Use it to |
|--------------|:----------:|-------------|-----------|
| `"standard"` (default) | Yes | The reports | Process a dataset from FASTQ to filtered count matrices. |
| `"geneext_only"` | GeneExt pass only | GeneExt | Build an extended annotation once, to reuse as `ref_gtf`. |
| `"external_pipeline_only"` | No | Vendor pipeline | Run only the vendor pipelines (Cell Ranger, split-pipe, BD Rhapsody). |
| `"post_mapping"` | No | The reports | Redo cell calling and everything after it on a mapped sample. |


### Full run (`standard`)

Pre-processing, mapping, cell calling, filtering and reporting, as in the [pipeline diagram](../img/Preprocs_Pipeline.png). This is the default, so `run_method` can be left out of the configuration file. Use the command under [Running the pipeline](#running-the-pipeline).


### Gene extension only (`geneext_only`)

Maps every sample once with STARsolo, runs [GeneExt](https://github.com/zolotarovgl/GeneExt) on the merged alignments and stops. The result is `<outdir>/gene_ext/geneext.gtf`, to pass back as `ref_gtf` on later runs. `perform_geneext` does not have to be set.

```
nextflow run main.nf -profile <institution_config>,conda -c </path/to/custom_parameters.config> \
    --run_method geneext_only
```

By default the mapping is reduced to what GeneExt needs: a BAM without cell tags, no count matrices and no cell calling. Set `geneext_bam_only = false` to map as a `standard` run would and still stop after GeneExt. See [The `geneext_only` run method](CONFIGURATION_PARAMETERS.md#the-geneext_only-run-method).


### External pipelines only (`external_pipeline_only`)

Runs the protocol's pre-processing and the enabled vendor pipelines, then stops. Nothing is mapped with STARsolo or alevin-fry, and no reports are made. A vendor pipeline only runs when it is enabled:

| Pipeline | Protocols | Enabled by |
|----------|-----------|------------|
| Cell Ranger | 10x Genomics, OAK | `perform_cellranger = true` |
| split-pipe | Parse Biosciences | `splitpipe_installation` and `splitpipe_conda_env` |
| BD Rhapsody | BD Rhapsody | `rhapsody_installation` |

```
nextflow run main.nf -profile <institution_config>,conda -c </path/to/custom_parameters.config> \
    --run_method external_pipeline_only
```

For installing split-pipe and BD Rhapsody, see [supported external pipelines](INSTALLATION_EXTERNAL_PIPELINES.md).


### Re-running after mapping (`post_mapping`)

Picks a finished run up where mapping ended. The previous run's STARsolo and alevin-fry results are read back in, cells are called again from the current samplesheet, and every step after that is redone: matrix filtering, saturation, read and per-cell metrics, ambient RNA removal, doublet detection, metacells, the dashboard and MultiQC. Demultiplexing, trimming, FastQC, indexing, mapping and GeneExt are not repeated.

Use it to correct the cell calling without mapping again:

- **A different number of expected cells:** edit `expected_cells` in the samplesheet. It bounds the window the second-derivative search looks for the knee in.
- **A cutoff chosen by hand:** set `cellfilter_method = "manual_cutoff"` and give every sample a `manual_cutoff` in the samplesheet.

See [Cell calling](CONFIGURATION_PARAMETERS.md#cell-calling) for both methods.

```
nextflow run main.nf -profile <institution_config>,conda -c </path/to/custom_parameters.config> \
    --run_method post_mapping \
    --previous_outdir </path/to/finished/results> \
    --cellfilter_method manual_cutoff
```

| Variable | Notes |
|----------|-------|
| `previous_outdir` | Results directory of the finished run. Defaults to `outdir`, which overwrites the cell-called outputs in place; set a different `outdir` to keep the original results. |
| `input` | The samplesheet, read again. Sample names must match the previous run; the FASTQs do not have to exist anymore. |
| `mapping_software`, `perform_geneext`, `geneext_downstream` | Must match the previous run, as they decide which analytical runs are picked up. They are recorded in its `pipeline_info/run_config_*.txt`. |
| `ref_gtf`, `ref_fasta`, `protocol` | Still required, for the read metrics. |
| `geneext_gtf` | Only with `perform_geneext`. Defaults to `<previous_outdir>/gene_ext/geneext.gtf`. |

Notes:

- The saturation, read metrics and Kraken steps read the previous run's BAMs, so they need a run made with `star_generateBAM = true`.
- An analytical run with missing or incomplete outputs is named in a warning and skipped; the other samples are still processed.
- With `perform_geneext`, the GeneExt runs are picked up (and the standard ones too with `geneext_downstream = "both"`). Without it, only the standard runs are processed, which needs a previous run that mapped them.


## Filtering UMI matrices (two passes)

Available with `run_method` `"standard"` or `"post_mapping"`.

1. Run with `perform_metacells = true`. The run ends with `<outdir>/filtering_report.html`.
2. Open the report and set the cell-, metacell- and gene-level filters for each sample. *Apply settings to all samples* copies them to the rest.
3. Download `metacell_selection.json` from the *Export selection* tab. It holds the samples whose settings were applied at all three levels.
4. Rerun the same command with `-resume` and the selection:

```
nextflow run main.nf -profile <institution_config>,conda -c </path/to/custom_parameters.config> -w </path/to/workdir> \
    -resume --metacell_selection metacell_selection.json
```

Only the filtering step runs. The filtered matrices are written to `<outdir>/metacell_filtering/<sample>/<sample>_final/`. See [Metacells and the filtering report](CONFIGURATION_PARAMETERS.md#metacells-and-the-filtering-report) for the report's settings.
