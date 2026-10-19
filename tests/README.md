# Test suite

Checks that validate the pipeline's environments and container images without
needing sequencing data. Each check is skipped rather than failed when the tool it needs is missing, so the suite is safe to run anywhere.  

Run everything with:

```bash
bash tests/run_tests.sh
```

List all tests that are available with:

```bash
bash tests/run_tests.sh --list
```

## Layout

```
tests/
├── run_tests.sh                # entry point; discovers and runs the checks
├── checks/                     # one file per check, discovered automatically
│   ├── alevin_usa.sh
│   ├── annotated_h5ad.sh
│   ├── conda_envs.sh
│   ├── containers.sh
│   ├── dashboard_geneext.sh
│   ├── dynamic_resources.sh
│   ├── metacell_filtering.sh
│   ├── mt_rrna_metrics.sh
│   ├── nextflow_versions.sh
│   ├── resource_efficiency.sh
│   └── velocity_matrix.sh
├── fixtures/
│   └── execution_trace_formats.txt  # hand-written trace covering awkward formats
└── lib/
    ├── common.sh                  # logging, result recording, summary
    ├── extract_containers.awk     # pulls image refs out of module definitions
    ├── make_annotation_fixture.py # generates both mappers' triplets plus doublet
    │                              # and CellSweep annotations to merge
    ├── make_geneext_fixture.py    # generates a synthetic GeneExt report and log
    ├── make_metacell_fixture.py   # generates a sample with a planted metacell assignment
    ├── report_harness.js          # runs filtering_report.html's script to export a selection
    ├── make_trace_fixture.py      # generates synthetic runs with a known exponent
    └── make_velocyto_fixture.py   # generates synthetic Velocyto and USA matrices
```

Logs land in `tests/.logs/<timestamp>/`, one file per test case. Those, along
with `tests/.cache/` and `tests/.envs/`, are gitignored.

## Running

```bash
bash tests/run_tests.sh                              # every check
bash tests/run_tests.sh conda_envs                   # a single check
bash tests/run_tests.sh conda_envs containers        # several

# Everything after `--` is passed through to the check
bash tests/run_tests.sh conda_envs -- --solve-only
bash tests/run_tests.sh containers -- --depth manifest

# Each check also runs standalone and documents its own options
bash tests/checks/containers.sh --help
```

The runner exits non-zero if any case failed, and prints the failures with a
pointer to the relevant log.

## `conda_envs`

Creates every `modules/**/environment.yml` with the solver behind each
conda-based profile in `nextflow.config` (`conda`, `mamba`, `micromamba`), so
unsolvable pins and packages dropped from a channel surface here instead of
halfway through a run.

A full build of all 39 environments is slow and disk-hungry. Useful variants:

```bash
# Fast: resolve dependencies without installing anything
bash tests/run_tests.sh conda_envs -- --solve-only

# Only the solver you actually use, only some modules
bash tests/run_tests.sh conda_envs -- --profiles mamba --only samtools

# Build several at once
bash tests/run_tests.sh conda_envs -- --jobs 4
```

Environments are created under `tests/.envs/` and removed once verified, unless
`--keep` is given. `--solve-only` needs a solver new enough to support
`--dry-run`; if yours is older, the log will show an unrecognised-argument
error.

Note the recipes are built directly from each `environment.yml`, so the channel
list comes from the file itself rather than from `conda.channels` in
`nextflow.config`.

## `containers`

Extracts the image references declared by the modules and checks that each one
resolves for the container profiles in `nextflow.config` (`docker`, `podman`,
`singularity`, `apptainer`, `wave`).

```bash
# Registry lookups only, no layers downloaded
bash tests/run_tests.sh containers -- --depth manifest

# Actually download the singularity images
bash tests/run_tests.sh containers -- --profiles singularity --depth pull
```

Notes:

- `--depth manifest` cannot check `oras://` references, since singularity has no
  way to inspect a remote image without fetching it. Use `--depth pull` for
  those. For ordinary registry images, a manifest check borrows `docker` when it
  is available.
- The `wave` profile builds images on demand through Seqera, so it needs the
  `wave` CLI and `TOWER_ACCESS_TOKEN`; without them it is skipped.

Pulled images are cached in `tests/.cache/containers/` and reused on later runs
unless `--force` is passed.

## `nextflow_versions`

Takes the pipeline through a dry run with every Nextflow version installed on
the machine, which is how you find out whether a release still works before
upgrading, and whether older ones still work after a change.

Versions come from `$NXF_HOME/framework/`, where Nextflow caches every version
it has been asked to run, plus whatever `nextflow -version` reports. To test
against another one, install it and it is picked up from then on:

```bash
NXF_VER=24.10.5 nextflow -version      # downloads and caches that version
bash tests/run_tests.sh nextflow_versions --list-versions
```

Each version goes through a ladder of increasingly demanding steps, selected
with `--depth`:

| depth | what it runs | what it catches |
| --- | --- | --- |
| `launch` | `nextflow -version` | the version is usable at all |
| `config` | `nextflow config` | config syntax, profiles, included configs |
| `parse` (default) | `nextflow run --help` | plugin compatibility (`nf-schema`), entry script compilation, parameter schema |
| `preview` | `nextflow run -preview` | the whole include graph and every channel |

```bash
# Quick: does each version still read the config?
bash tests/run_tests.sh nextflow_versions -- --depth config

# Full dry run of two versions, downloading them if needed
bash tests/run_tests.sh nextflow_versions -- \
    --versions 24.10.5,25.04.6 --install \
    --depth preview --config tests/test_parsebio.config
```

Notes:

- Only `preview` resolves the whole workflow, and it is the only depth that
  needs data: it runs the same dry-run as [`submit_nextflow.sh`](../submit_nextflow.sh)
  and so needs a config with paths you have access to. Without `--config` it
  uses the first `tests/*.config` it finds, and skips if there is none.
- Versions below the minimum in `manifest.nextflowVersion` are skipped, since
  Nextflow refuses to run the pipeline at all with those. Pass
  `--ignore-manifest` to try them anyway, which is how you check whether the
  declared minimum is still honest.
- `--profile` takes the same string as a real run, so use `--profile crg,conda`
  if that is what you launch with. The default is `conda`.
- Each version runs in its own directory under `tests/.logs/`, so `.nextflow.log`
  and `work/` never land in the repository. They are removed afterwards unless
  `--keep` is given.

## `resource_efficiency`

Validates [`bin/resource_efficiency.py`](../bin/resource_efficiency.py), which turns
execution traces into scheduled-resource recommendations. A parsing error there does
not raise an exception — it produces a plausible-looking number that is wrong, and the
next run gets sized from it. So the check asserts on the numbers rather than on the
exit status.

There is no real trace to check against on a development machine, so it builds two
kinds of input:

- `tests/fixtures/execution_trace_formats.txt`, written by hand, covering the awkward
  things Nextflow emits: `1.2 GB` and `512 MB`, `2h 3m 4s` and `250ms` and `1d 2h`,
  `540.5%` (over 100, because `%cpu` is summed across cores), `-` for unmeasured
  fields, a `CACHED` row of nothing but dashes, a task killed with exit 137 followed
  by its successful retry, an aliased process name, an unlabelled process and a
  process no module declares.
- `tests/lib/make_trace_fixture.py`, which generates whole runs whose memory follows a
  **known** power law, so the recovered scaling exponent can be compared against the
  one that was asked for. Deterministic under `--seed`.

```bash
bash tests/run_tests.sh resource_efficiency

# Use a specific interpreter, and embed a different exponent
bash tests/run_tests.sh resource_efficiency -- \
    --python "$(command -v python3)" --exponent 1.0 --tolerance 0.2

# Keep the generated fixtures to look at them
bash tests/run_tests.sh resource_efficiency -- --keep
```

Notes:

- Needs only a `python3`; the tool's analysis path is standard library only and the
  check runs it with `--no-plots`. It skips cleanly if no interpreter is found.
- `default_fields` covers a trace written with Nextflow's *default* `trace.fields`,
  which records no `cpus`/`memory`/`time`. Without a denominator every efficiency is
  blank and both overview charts come out empty, so the tool substitutes the declared
  `conf/base.config` values; this case asserts it does, and that it still tells the
  reader the values were inferred. It also covers pairing the run config by nearest
  timestamp, which is what a `-resume` leaves behind.
- `label_coverage` asserts that every process under `modules/` maps to a known tier
  and that the four aliased inclusions resolve. This is the case that fails when a new
  module is added with a label the mapper cannot see.
- `base_config` asserts that each tier in `conf/base.config` parses to a real
  cpus/memory/time triple — it catches the parser silently returning `None` after a
  formatting change to that file.
- `coverage_guard` asserts that a tier only partly exercised by the runs is never
  *lowered*, since the members that did not run would be under-provisioned.
- `emitted_config` runs `nextflow config` over the generated file and skips, rather
  than fails, when Nextflow is not installed.

## `dynamic_resources`

Unit-tests [`lib/BcaResources.groovy`](../lib/BcaResources.groovy), which turns the
size of a process's input into its memory request, and an allocation into STAR's
`--limitBAMsortRAM`. That code decides what tasks ask the scheduler for, and a
mistake in it does not throw — it produces a plausible number that is too small, and
the task is killed part-way through a long run.

The cases pin the properties that have to hold regardless of the coefficients: a
process with no entry in `params.dynamic_memory` gets *exactly* its label value, a
malformed or absent entry falls back rather than throwing (Nextflow's config returns
an empty `ConfigObject` for a missing key, not `null`, which is easy to get wrong), a
large resident reference raises the floor so a small input can never be handed less
memory than the reference alone needs, and `task.attempt` still escalates. It also
asserts that every entry actually shipped in `nextflow.config` is complete, since an
entry missing `ref_gb` or `mem_gb` silently does nothing.

Two cases work on a `TaskPath`, which is what an input really is once it reaches a
process: a path named after the staged file, whose file system provider is not the
one that owns it, so every `java.nio.file.Files` call on it throws. Measured
carelessly through one of those, *every* input reads as 0 bytes and every request
falls back to its label without a word in the log — so the cases construct one and
assert the size comes back the same as for the path it stands for.

The `bamsort_*` cases cover the STAR sort budget: the allocation less the resident
genome index and STAR's own working memory, a pinned `params.star_limitBAMsortRAM`
used verbatim, and the two degenerate cases — no allocation to derive from, and an
index too large for the allocation to hold.

```bash
bash tests/run_tests.sh dynamic_resources
bash tests/run_tests.sh dynamic_resources -- --java /path/to/java --jar /path/to/nextflow-one.jar
```

Groovy comes from the Nextflow fat jar that Nextflow caches under
`$NXF_HOME/framework/`, so no extra toolchain is needed — but the check skips
cleanly when neither a JDK nor a cached jar is present.

## `dashboard_geneext`

Checks the gene-extension tab of `dashboard.html`, which
[`bin/generate_dashboard.py`](../bin/generate_dashboard.py) fills by re-reading the
JSON payload GeneExt embeds in its own HTML report, falling back to GeneExt's
plain-text log when that report is absent. Neither failure mode is loud: a payload
that cannot be located yields an empty object, which just hides the tab, and a
fallback that matches nothing yields a tab full of blanks.

GeneExt itself needs a genome, an annotation and an aligned BAM, so
[`lib/make_geneext_fixture.py`](lib/make_geneext_fixture.py) reconstructs the two
files the dashboard reads with distinctive, known numbers — including the one-gene
discrepancy GeneExt has between its report and its log, which is what tells the two
code paths apart in the assertions.

The cases also pin the two reductions the dashboard makes deliberately: GeneExt's
per-gene extension table and its orphan-peak BED are *not* carried over, since they
dominate the report's size and are one click away in the report itself.

```bash
bash tests/run_tests.sh dashboard_geneext
bash tests/run_tests.sh dashboard_geneext -- --keep    # keep the generated dashboards
```

## `metacell_filtering`

Checks the metacell filtering loop end to end:

- [`bin/run_metacells.py`](../bin/run_metacells.py) summarises each metacell.
- [`bin/filtering_report.html`](../bin/filtering_report.html) turns the user's choices into a
  selection.
- [`bin/apply_metacell_filter.py`](../bin/apply_metacell_filter.py) writes the matrices that
  selection asks for.

Every link can fail silently:

- A mitochondrial percentage computed after the mitochondrial genes were excluded from the
  grouping reads 0 everywhere.
- A CellSweep annotation lifted by position lands on the wrong cells.
- A selection applied to a regrouped sample keeps a plausible set of the wrong cells.

[`lib/make_metacell_fixture.py`](lib/make_metacell_fixture.py) therefore plants the assignment,
so Metacell2 is not needed for most cases. It also plants every number the report shows:

- the ambient object's barcodes are shuffled among empty droplets;
- the domain annotation is keyed by versioned transcript and protein ids;
- one gene is named `</script><script>alert(1)` and another like the report's placeholder.

Where a JavaScript runtime is found, [`lib/report_harness.js`](lib/report_harness.js) runs the
report's own script in a stubbed DOM. It applies a set of choices and exports the selection,
which the apply cases then check. So what is tested is what the page really downloads,
including a download → load round trip. Without a runtime, those cases are skipped and an
equivalent selection is built in Python.

`mc2_smoke` runs the real Metacell2 grouping twice with the same seed and expects the same
fingerprint. It is skipped where `metacells` cannot be imported. metacells publishes no
Windows wheels, so on Windows run it from WSL or a Linux node.

```bash
bash tests/run_tests.sh metacell_filtering
bash tests/run_tests.sh metacell_filtering -- --python /path/to/env/bin/python --node node --keep
```

## `mt_rrna_metrics`

Checks the read-level mtDNA and rRNA metrics of
[`bin/calculate_read_metrics.py`](../bin/calculate_read_metrics.py) (CALC_READ_METRICS) and
[`bin/per-cell_images.py`](../bin/per-cell_images.py) (PERCELL_METRICS).

A BAM holds one record per alignment. A count of records therefore weights every
multimapper by its `NH`, and the percentages stay plausible while measuring the wrong thing.
The check writes a small SAM by hand, so that every expected number can be counted on paper.
It contains:

- a nuclear multimapper with a secondary alignment on chrM;
- an mtDNA multimapper on mitochondrial rRNA;
- reads sense, antisense and outside the mitochondrial genes;
- an rRNA read on an added rRNA contig;
- reads from a called and an uncalled barcode.

The cases assert exact values for:

- the library-level rows;
- the split of mtDNA reads by strand;
- the called-cell rows;
- the antisense rows, from a hand-written `CellReads.stats`, and their absence for an
  unstranded run;
- the per-barcode reads (`barcode_reads.tsv.gz`), and that the called barcodes add up to the
  called-cell rows;
- the per-cell table: mapped reads, mtDNA %, rRNA %, intronic % from `CellReads.stats`,
  and unspliced % from Velocyto. `per-cell_images.py` reads the per-barcode reads of the
  stranded run, as `PERCELL_METRICS` reads `CALC_READ_METRICS`', so the case also checks
  that the two scripts fit together.

The check needs samtools, to build the fixture BAM, and a Python with pysam, which has no
Windows build, so on Windows run it from WSL. The per-cell case also needs pandas, scipy and
matplotlib; it no longer needs pysam itself. Each case is skipped when its tools are missing.

```bash
bash tests/run_tests.sh mt_rrna_metrics
bash tests/run_tests.sh mt_rrna_metrics -- --python /path/to/env/bin/python --keep
```

## `velocity_matrix`

Checks the intronic / RNA-velocity outputs produced under `perform_velocity = true`:
[`bin/subset_matrices_to_cells.py`](../bin/subset_matrices_to_cells.py),
[`bin/velocity_matrices_to_h5ad.py`](../bin/velocity_matrices_to_h5ad.py) and the
`--counts U` path through
[`bin/collapse_alevin_usa.py`](../bin/collapse_alevin_usa.py).

Every failure mode here is silent. Subsetting the velocity matrices on a UMI cutoff of
their own rather than the `GeneFull_Ex50pAS` cell call gives a perfectly plausible matrix
describing the wrong cells. Transposing one layer and not the others, or mapping
alevin-fry's spliced block onto the unspliced layer, gives an object of exactly the right
shape carrying the wrong numbers. None of it raises.

[`lib/make_velocyto_fixture.py`](lib/make_velocyto_fixture.py) therefore writes counts
that identify their own layer, gene and cell on sight — `spliced` in the hundreds,
`unspliced` in the tens, `ambiguous` in the units — so a swapped layer cannot pass. It
emits both mappers' layouts: STARsolo's genes × cells matrices and alevin-fry's cells ×
columns USA matrix in alevin-fry's real column layout (bare gene IDs, then `-U`, then
`-A`). One gene ID ends in `-A`, which a split by suffix would mistake for another gene's
ambiguous column. The alevin-fry h5ad is checked against alevin-fry's velocity convention
(`spliced` = S + A, `unspliced` = U).

The `h5ad_layers` cases are skipped where `anndata` is not importable; the rest need only
numpy, pandas and scipy.

```bash
bash tests/run_tests.sh velocity_matrix
bash tests/run_tests.sh velocity_matrix -- --keep    # keep the generated matrices
```

## `alevin_usa`

Checks every script that reads alevin-fry's USA matrix through
[`bin/alevin_usa.py`](../bin/alevin_usa.py): the collapse
([`collapse_alevin_usa.py`](../bin/collapse_alevin_usa.py)) for every `alevin_usa_counts`
value, the cell-calling helper
([`secondderiv_alevin.py`](../bin/secondderiv_alevin.py)) `umis` and `filter`, and the gene
count in both dashboards. It reuses the velocity fixture, so every expected count is the
sum of known S, U and A values.

The failures it guards against were silent: a collapse that did not recognise the layout
copied the matrix through with three columns per gene, and cell calling then counted every
column, so `alevin_usa_counts` had no effect and genes were counted up to three times in
the statistics. A matrix that is not in the USA layout (plain gene columns, or `-S`/`-U`/`-A`
names) must now fail without writing output. Needs only numpy, pandas and scipy.

```bash
bash tests/run_tests.sh alevin_usa
```

## Adding a check

Drop an executable `*.sh` file into `tests/checks/`. The runner picks it up by
filename, so no registration is needed. Follow the shape of the existing ones:

```bash
#!/usr/bin/env bash
# description: One line shown by --list

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/common.sh"

record PASS "some case"
record FAIL "another case" "why it failed"
record SKIP "third case" "why it was skipped"

finish_check
```

`record` writes results to a shared file, which is what lets the runner
summarise across checks and lets a check run its cases in parallel. Report a
missing prerequisite as `SKIP` when it was auto-detected, and as `FAIL` only
when the user explicitly asked for it.
