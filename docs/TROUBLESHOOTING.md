# Common Issues & Troubleshooting

This page documents common issues encountered when running the Nextflow pipeline, along with their causes and fixes.

If you fix an issue, you can usually continue the pipeline using:

```bash
nextflow run <pipeline> -resume
```

Always inspect `.nextflow.log` for additional debugging details.


---

## Exit status 126 – Permission denied on script in `bin/`

### Error

The pipeline fails with an error similar to:

```
Command exit status:
  126

Command error:
  .command.sh: line 14: \
  /path/to/bca_preprocessing/bin/salmon_create_splici_ref.R: Permission denied
```

Execution stops while trying to run an R script located in the pipeline `bin/` directory.

Inspecting the `bin/` directory shows:

```
-rw-rw---- salmon_create_splici_ref.R
```

---

### Solution

Exit status **126** indicates that a file exists but **cannot be executed**.

Nextflow automatically prepends the pipeline’s `bin/` directory to the `PATH` and assumes that all scripts inside it are executable. In this case, the script `bin/salmon_create_splici_ref.R` is missing execute permissions

Make the scripts within the bin/ folder executable from the root of the pipeline repository:

```bash
chmod -R +x bin/
```

Verify the permissions:

```bash
ls -lah bin/
```

Expected output is similar to:

```
-rwxrwx--- salmon_create_splici_ref.R
...
```
