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
  Each run writes its execution trace into `execution_reports/` in the directory you
  launched it from, and its run config, samplesheet and parameters into
  `<outdir>/pipeline_info/`. All are named with the run's timestamp, so repeated runs
  accumulate instead of overwriting.


2. **Run the tool** against your results directory:

   ```bash
   bin/resource_efficiency.py --results /path/to/outdir
   ```

   The tool looks for the traces in `execution_reports/` in the pipeline folder (where
   `submit_nextflow.sh` launches from) and in the current folder. It only takes the
   traces of runs recorded in that outdir's `pipeline_info/`, so runs with other outdirs
   launched from the same folder are left out. If you launched from somewhere else, say
   where:

   ```bash
   bin/resource_efficiency.py --results /path/to/outdir --traces /launch/dir/execution_reports
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

## Dynamic time and memory

A label like `process_low` gives every step the same memory, and it has to be big
enough for the largest dataset you will ever run. On everything smaller, most of that
is reserved and left idle.

For these steps memory tracks the size of their input, so they scale with it instead of
using a fixed amount of resources. The numbers live in `dynamic_memory` in 
[`nextflow.config`](../nextflow.config), and you can refresh them based on your own runs:

```bash
bin/resource_efficiency.py --results /path/to/outdir --dynamic
```

That writes a `resources_dynamic_*.config` you can pass with `-c`, or copy into
`nextflow.config`. It holds the `dynamic_time` and `dynamic_memory`.

`dynamic_cpus` steps a CPU count up above an input size. It has one entry,
`STARSOLO_ALIGN: [ above_gb: null, cpus: 16 ]`, which does nothing while `above_gb` is
`null`: the size at which 16 threads pay off has to come from a thread sweep on the
cluster. With a threshold set, a larger input gets 16 cores, and shorter runs.


## Tips

- **The more runs, the better.** A single run still gives useful numbers, but seeing a
  few runs on different-sized datasets lets the tool predict resource needs for future
  larger or smaller datasets too.
- **Don't delete `execution_reports/`** in your launch folder, or `pipeline_info/` in
  your results folders. The first holds the usage data, the second tells the tool which
  run each trace belongs to and what data it processed.
- **A step marked "insufficient data"** hasn't run enough times yet for a reliable
  recommendation; its current setting is left unchanged.
