# Test data

This folder contains small, subsampled datasets with a matching parameter config and samplesheet. Use them to check that the pipeline and your compute setup work before you run it on your own data.

```
examples/
├── 10x/
│   ├── conf/        # parameter config + samplesheet
│   └── test_data/   # subsampled FASTQs + reference
└── mars-seq/
    ├── conf/        # parameter config + samplesheet
    └── test_data/   # subsampled FASTQs + reference
```

## Running a test

Run the commands from the repository root, because the paths in the configs and samplesheets are relative to it.

Add the `test` profile to your usual profiles. It caps every process at 2 CPUs, 4 GB of memory and 1 hour ([conf/test.config](../conf/test.config)), and it turns off retries, so a failure on the toy data shows up straight away.

**10x genomics data:**

```bash
nextflow run main.nf \
    -profile crg,conda,test \
    -c examples/10x/conf/spis_10x_parameters.config \
    -ansi-log false
```

Output is written to `examples/10x/output/`.

Replace `crg` with the profile for your own institution, and `conda` with another software profile if you need one (for example `singularity`).

### Or: on a SLURM cluster via `submit_nextflow.sh`

[submit_nextflow.sh](../submit_nextflow.sh) already contains the test command. Comment out the default `nextflow run ...` line, and uncomment the test line under *Test the pipeline*:

```bash
nextflow run -profile crg,conda,test -c examples/10x/conf/spis_10x_parameters.config -ansi-log false "$@" & pid=$!
```

Then submit it:

```bash
sbatch submit_nextflow.sh main.nf
```

## Data sources

Both test datasets were **subsampled** from the public releases below. They are only meant for testing the pipeline and are not suitable for biological analysis.

### 10x Genomics (*Stylophora pistillata*)

| | |
|---|---|
| Publication | <https://www.nature.com/articles/s41586-025-09623-6> |
| SRA | [SRR32332771](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR32332771&display=metadata) |

### MARS-seq (*Nematostella vectensis*)

| | |
|---|---|
| Publication | <https://www.cell.com/cell/fulltext/S0092-8674(18)30596-8> |
| SRA | [SRR6502902](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR6502902&display=download) |
