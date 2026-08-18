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
│   ├── conda_envs.sh
│   ├── containers.sh
│   ├── nextflow_versions.sh
│   └── resource_efficiency.sh
├── fixtures/
│   └── execution_trace_formats.txt  # hand-written trace covering awkward formats
└── lib/
    ├── common.sh               # logging, result recording, summary
    ├── extract_containers.awk  # pulls image refs out of module definitions
    └── make_trace_fixture.py   # generates synthetic runs with a known exponent
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
