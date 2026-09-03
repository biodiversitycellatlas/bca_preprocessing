# Resource Tuning

Find out how much CPU, memory and time the pipeline runs actually used, and get a
ready-to-use config that fixes over- or under-sized requests for your next run.

---

### **Table of Content**:
- [Quick start](#quick-start)
- [Reading the report](#reading-the-report)
- [Applying the recommendations](#applying-the-recommendations)
- [Dynamic memory](#dynamic-memory)
- [Tips](#tips)

---

## Quick start

1. **Run the pipeline as normal.** 
  Each run automatically creates the files required into `<outdir>/pipeline_info/`, named with the run's
  timestamp, so repeated runs into the same directory accumulate instead of overwriting. 


2. **Run the tool** against your results directory:

   ```bash
   bin/resource_efficiency.py --results /path/to/outdir
   ```

   This prints a summary table and writes an HTML report you can open in a browser.
   To skip the report and just get the table (no extra setup needed):

   ```bash
   bin/resource_efficiency.py --results /path/to/outdir --no-plots
   ```

3. **Got more than one run?** Point at all of them for more accurate recommendations,
   especially if the datasets were different sizes:

   ```bash
   bin/resource_efficiency.py -r /path/run_A -r /path/run_B -r /path/run_C
   # or, if they all live under one parent folder:
   bin/resource_efficiency.py --results-root /path/to/runs
   ```

## Reading the report

The table has one row per pipeline step, with the biggest resource users at the top.

| Column | What it tells you |
| --- | --- |
| **mem eff** | How much of the requested memory was actually used. Close to 100% means the next run might run out of memory. Very low means memory is being wasted. |
| **cpu eff** | How much of the requested CPUs were kept busy. Low means fewer CPUs would do just as well. |
| **retr** | Tasks that failed and had to retry. Any number above 0 is a sign the request was too small. |
| **-> recommendation** | The suggested new request for cpus / memory / time. |

The HTML report additionally shows a chart per step: how much it used, plotted against
the size of the data it processed, so you can see at a glance whether a step scales with
dataset size or has a mostly fixed cost.

## Applying the recommendations

```bash
bin/resource_efficiency.py --results /path/to/outdir --apply
```

This will generate a new configuration file `conf/resources_tuned.config`.
To apply these new settings to your next run, add it as a new config file to your nextflow command:

```bash
nextflow run
    -profile <profile>
    -c </path/to/custom_parameters.config>
    -c conf/resources_tuned.config
    ...
```

## Dynamic memory

A label like `process_low` gives every step the same memory, and it has to be big
enough for the largest dataset you will ever run. On everything smaller, most of that
is reserved and left idle.

For eight steps we noticed that memory tracks the size of their input, so
they scale with it instead of using a fixed amount of resources:

| Step | Dependent on | Scale |
| --- | --- | --- |
| `STARSOLO_ALIGN` | FASTQ size | ~5 GB small, ~160 GB very large |
| `STARSOLO_INDEX` | Reference genome size |  ~4 GB small, ~34 GB large |
| `MTX_TO_H5AD` | Matrix size, plus CellSweep's denoised matrix where it ran | ~3 GB small, up to 96 GB |
| `MTX_TO_10X` | Matrix size | ~3 GB small, up to 48 GB |
| `DOUBLET_FILTER` | Matrix size | ~3 GB small, up to 48 GB |
| `SATURATION_TABLE` | Filtered BAM size | ~13 GB small, up to 128 GB |
| `SUBSET_VELOCYTO_MATRICES` | Velocyto matrix size | ~3 GB small, up to 48 GB |
| `VELOCITY_H5AD` | Velocyto or USA matrix size | ~3 GB small, up to 48 GB |

`DOUBLET_FILTER` only runs with `perform_doublet_filtering = true`, and the last two only
with `perform_velocity = true`. `MTX_TO_H5AD` reads two matrices where CellSweep ran — the
counts and the denoised counts it merges in as a layer — which is why its ceiling is the
higher one; its coefficients were measured before that and are worth re-measuring.

The numbers live in `dynamic_memory` in [`nextflow.config`](../nextflow.config), and you can refresh them
based on your own runs:

```bash
bin/resource_efficiency.py --results /path/to/outdir --dynamic
```

That writes a `resources_dynamic_*.config` you can pass with `-c`, or copy into
`nextflow.config`.

`perform_velocity = true` is worth a note of its own. It adds a third feature structure
to STARsolo's solo post-map step, so `STARSOLO_ALIGN` needs more memory — but its entry
scales on **FASTQ size**, which does not change when a feature is added. The anchor
(`mem_gb` at `ref_gb`) therefore under-provisions with velocity on. The retry ladder
absorbs a single overrun, but if `STARSOLO_ALIGN` starts retrying on memory once you
enable velocity, raise `dynamic_memory.STARSOLO_ALIGN.mem_gb` (and `cap_gb` if you are
hitting it) rather than assuming the fit is wrong.

One thing to watch when copying an entry across: `--dynamic` calibrates `ref_gb` from
the trace's `rchar`, the total bytes a task read, while the module measures the file
as it was staged. For a step that reads its input more than once the two differ —
`SATURATION_TABLE` reads its BAM roughly seven times — so a generated `ref_gb` has to
be converted back to staged size first, or the anchor sits far too high and every
request collapses to the floor.


## Tips

- **The more runs, the better.** A single run still gives useful numbers, but seeing a
  few runs on different-sized datasets lets the tool predict resource needs for future
  larger or smaller datasets too.
- **Don't delete `pipeline_info/`** in your results folders, that's where the usage
  data comes from.
- **A step marked "insufficient data"** hasn't run enough times yet for a reliable
  recommendation; its current setting is left unchanged.
