# Resource Tuning

Find out how much CPU, memory and time the pipeline runs actually used, and get a
ready-to-use config that fixes over- or under-sized requests for your next run.

---

### **Table of Content**:
- [Quick start](#quick-start)
- [Reading the report](#reading-the-report)
- [Applying the recommendations](#applying-the-recommendations)
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

## Tips

- **The more runs, the better.** A single run still gives useful numbers, but seeing a
  few runs on different-sized datasets lets the tool predict resource needs for future
  larger or smaller datasets too.
- **Don't delete `pipeline_info/`** in your results folders, that's where the usage
  data comes from.
- **A step marked "insufficient data"** hasn't run enough times yet for a reliable
  recommendation; its current setting is left unchanged.
