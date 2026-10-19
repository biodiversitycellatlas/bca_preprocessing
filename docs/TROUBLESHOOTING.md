# Common Issues & Troubleshooting

This page describes common issues encountered when running the Nextflow pipeline, along with their causes and fixes.

If you fix an issue, you can usually continue the pipeline using:

```bash
nextflow run <pipeline> -resume
```

Always inspect `.nextflow.log` for additional debugging details.


---

## Exit Statuses

### Exit status 126 – Permission denied on script in `bin/`

The pipeline fails with an error similar to:

```
Command exit status:
  126

Command error:
  .command.sh: line 14: \
  /path/to/bca_preprocessing/bin/salmon_create_splici_ref.R: Permission denied
```

Exit status **126** indicates that a file exists but **cannot be executed**.

Nextflow automatically prepends the pipeline’s `bin/` directory to the `PATH` and assumes that all scripts inside it are executable. In this case, the script `bin/salmon_create_splici_ref.R` is missing execute permissions.

Make the scripts within the bin/ folder executable from the root of the pipeline repository:

```bash
chmod -R +x bin/
```

## Tool-internal memory ceilings

A scheduler allocation is not the same thing as the ceiling a tool applies to itself.
STAR has one: `--limitBAMsortRAM`, the buffer it sorts the BAM in. Left at `0`, STAR
sets it from the size of the genome index, which has nothing to do with the memory
the job was actually given — on a generous allocation it sorts in far less RAM than
it could, and on a tight one it reserves memory the scheduler never granted and the
task is killed.

The pipeline therefore derives it per task, in
[`lib/BcaResources.groovy`](../lib/BcaResources.groovy): the allocation, minus the
genome index STAR keeps resident for the whole run, minus a slack allowance for its
alignment buffers and the solo structures. The task log records the resulting
arithmetic on one line:

```
limitBAMsortRAM: 91214942208 = allocation 137438953472 - genome index 41929244672 - slack 4294967296
```

Set [`star_limitBAMsortRAM`](CONFIGURATION_PARAMETERS.md) to a non-zero number of
bytes to pin it instead — which is how to answer a STAR error that names the exact
figure it wants.