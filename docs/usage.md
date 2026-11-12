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
5. [Combining samples](#combining-samples)

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

Samples in a `filtering_report` group are filtered together, as one entry of the report, see [Combining samples](#combining-samples).


## Combining samples

Samplesheet rows can be combined at three points. Which one is right depends on whether the rows hold the same cells.

| Mechanism | Combined | Use it when the rows hold |
|-----------|----------|---------------------------|
| The same `sample` name on several rows | The FASTQs, before mapping | **The same cells**: several lanes, or a top-up sequencing of one library. |
| `qc_report` column | Only the numbers in `dashboard.html` | **Different cells** whose quality you want to judge together. |
| `filtering_report` column | The cell-called objects, just before Metacell2 | **Different cells** you want to filter and analyse as one set of metacells. |

Rows that hold different cells must not share a `sample` name. Their barcodes overlap, so merging their reads before mapping turns two cells with the same barcode into one, and mixes their empty droplets. This happens with:
- one sample loaded on several 10x or BD channels;
- the sublibraries of a Parse, ScaleBio or sci-RNA-seq3 experiment, which reuse the same barcode combinations;
- the plates of a ScaleBio Extended Throughput run;
- the aliquots of an overloaded 10x channel (OAK).

Give such rows their own `sample` names and combine them with the two columns instead. Both are optional. An empty value leaves a sample on its own, and the two columns are independent of each other:

```
sample,fastq_cDNA,fastq_BC_UMI,expected_cells,rt,qc_report,filtering_report
nvec_well1_larva_rep1,...,...,2500,A1-A6,nvec_well1,nvec_larva
nvec_well1_larva_rep2,...,...,2500,A1-A6,nvec_well1,nvec_larva
nvec_well2_larva,...,...,2500,B1-B6,nvec_well2,nvec_larva
nvec_well3_adult,...,...,2500,C1-C6,nvec_well3,nvec_adult
nvec_well4_adult,...,...,2500,D1-D6,nvec_well4,nvec_adult
```

Here the two replicates of well 1 are shown together in the dashboard. The metacells and the filtering are done per stage: wells 1 and 2 as `nvec_larva`, wells 3 and 4 as `nvec_adult`.

### Which mechanism to use

| Situation | `qc_report` | `filtering_report` |
|-----------|-------------|--------------------|
| One sample on several 10x or BD channels | the sample | the sample, or its condition |
| Sublibraries of one Parse, ScaleBio or sci-RNA-seq3 experiment, one row per sublibrary | the sample | the sample |
| ScaleBio Extended Throughput plates (`sampleX_AP1`, `sampleX_AP2`) | `sampleX` | `sampleX` |
| Aliquots of one overloaded 10x channel (OAK) | the channel | the sample or condition |
| The wells of a sample split into rows, to see the quality of each well | the well or replicate | the condition |
| Biological replicates of one condition | the condition, to compare the replicates | the condition, for more cells per Metacell2 run and one filtering applied to all |
| Stages or tissues meant to share metacells (an atlas) | - | one group for all of them |

Do not group:
- **Two sequencing runs of the same library under `filtering_report`.** The same cells would enter Metacell2 twice. Give the runs one `sample` name instead. The pooling step warns when two samples share more than half of their cell barcodes. Unrelated captures share far fewer, a few percent at most.
- **Samples that need very different filtering.** The report's thresholds act on the pooled metacells, so one setting applies to the whole group.
- **Samples of different species or references.** Each needs its own samplesheet anyway, and the pooling refuses samples counted against different genes.

Pooled metacells can be made of one sample only, from a batch effect or from a cell type only that sample has. The report's *Main sample %* column shows this per metacell.

Rules, checked when the run starts:
- every row of a sample gives the same groups;
- a group cannot have the name of a sample, or end in an analytical-run suffix such as `_starsolo`;
- the name of a sample in a `filtering_report` group cannot contain `:`.

A group of one sample is allowed, and runs under the group's name.

### `qc_report`: pooled in the dashboard

Every sample is mapped, its cells called, and its ambient RNA and doublets assessed on its own. Only `dashboard.html` brings a group together, once per analytical run (a group's `_starsolo` runs with each other, never with its `_alevinfry` runs):
- a **Sample groups** box, with the members' numbers and a pooled row, every member's barcode-rank curve on one plot and the mapped reads per called cell of every member side by side;
- a *Group* column and the pooled rows in the summary table;
- the members listed under their group in the sample selector.

How the pooled row is built:

| Value | Pooled as |
|-------|-----------|
| Input reads, called cells, reads needed for the target saturation | summed |
| % mapped reads, mean reads per cell | recomputed from the sums |
| Noise, % intronic reads | members' values weighted by their uniquely mapped reads |
| Median UMIs per cell | the median over every member's called cells |
| Saturation | approximate: weighted by input reads |
| rRNA and mtDNA % | approximate: weighted by uniquely mapped reads |
| Median genes per cell | approximate: members' medians weighted by their cells |

Approximate values are marked ≈ in the table. Saturation and barcode-rank curves cannot be pooled, so they are shown side by side. MultiQC and `summary_results/mapping_stats.tsv` stay per sample.

### `filtering_report`: pooled for metacells and filtering

Every sample goes through cell calling, ambient RNA removal (CellSweep) and doublet detection on its own. These belong to one capture:
- the ambient RNA profile is that of one channel or reaction;
- doublets form only within one capture;
- scDblFinder expects a doublet rate that grows with the number of cells it is given, so three pooled captures of 5,000 cells would be modelled at about 15 % doublets instead of about 5 %.

Just before Metacell2, `MERGE_GROUP_H5AD` pools the samples' cell-called objects, per analytical run. The doublet calls and CellSweep annotations travel with each cell. From there on, the group takes its samples' place:
- one Metacell2 run, under `metacells/<group><suffix>/`;
- one entry in `filtering_report.html`, listing its samples. Each cell's reads come from its own sample's `CellReads.stats`, and the read metrics are the samples' values weighted by their cells;
- one entry in `metacell_selection.json`, keyed by the group;
- after the selection is applied, the group's filtered matrices in `metacell_filtering/<group><suffix>/`, and each sample's kept cells under `members/<sample>/`.

In the pooled object, the cells are named `<sample>:<barcode>`, so that cells of different samples with the same barcode stay apart. `obs["sample"]` gives each cell's sample. The files under `members/` use the original barcodes again. See [Metacells and filtering](output.md#metacells-and-filtering) for the files.

### Changing the groups of a finished run

The groups are read from the samplesheet when the run starts, and are not part of any step before Metacell2. So the groups can be changed after a run, for example after looking at the per-sample dashboard. Edit the two columns and rerun with `-resume`. Only these steps run again:
- the pooling;
- Metacell2 for the groups that changed, and for samples taken out of a group;
- the two reports.

A `post_mapping` run reads the groups from the current samplesheet as well.

A selection made before the groups changed no longer applies to the changed groups. Export a new one from the new report.
