# `biodiversitycellatlas/bca_preprocessing`: Contributing Guidelines

Hi there!
Many thanks for taking an interest in improving biodiversitycellatlas/bca_preprocessing.

We try to manage the required tasks for biodiversitycellatlas/bca_preprocessing using GitHub issues, you probably came to this page when creating one.
Please use the pre-filled template to save time.

However, don't be put off by this template - other more general issues and suggestions are welcome!
Contributions to the code are even more welcome ;)

## Contribution workflow

If you'd like to write some code for biodiversitycellatlas/bca_preprocessing, the standard workflow is as follows:

1. Check that there isn't already an issue about your idea in the [biodiversitycellatlas/bca_preprocessing issues](https://github.com/biodiversitycellatlas/bca_preprocessing/issues) to avoid duplicating work. If there isn't one already, please create one so that others know you're working on this
2. [Fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the [biodiversitycellatlas/bca_preprocessing repository](https://github.com/biodiversitycellatlas/bca_preprocessing) to your GitHub account
3. Make the necessary changes / additions within your forked repository following [Pipeline conventions](#pipeline-contribution-conventions)
4. Use `nf-core pipelines schema build` and add any new parameters to the pipeline JSON schema (requires [nf-core tools](https://github.com/nf-core/tools) >= 1.10).
5. Submit a Pull Request against the `dev` branch and wait for the code to be reviewed and merged

If you're not used to this workflow with git, you can start with some [docs from GitHub](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests) or even their [excellent `git` resources](https://try.github.io/).

## Testing

The checks in [`tests/`](../tests/README.md) validate the pipeline's environments
and container images and need no real data:

```bash
bash tests/run_tests.sh           # run everything
bash tests/run_tests.sh --list    # see the available checks
```

- `conda_envs` builds every `modules/**/environment.yml` with each conda-based
  profile (`conda`, `mamba`, `micromamba`), catching dependency specs that no
  longer solve.
- `containers` checks that every container image declared by the modules still
  resolves for each container profile.
- `nextflow_versions` dry-runs the pipeline with each Nextflow version installed
  on your machine, so a release that breaks the config, the plugins or the
  workflow is caught before you upgrade.

They are slow at full depth, so run the fast variants while iterating:

```bash
bash tests/run_tests.sh conda_envs -- --solve-only
bash tests/run_tests.sh containers -- --depth manifest
bash tests/run_tests.sh nextflow_versions -- --depth config
```

A check is skipped, not failed, when the tool it needs is not installed, so the
suite is safe to run on any machine. If you add a new step to the pipeline,
please make sure it still passes. See [`tests/README.md`](../tests/README.md) for
how to add a check of your own.

## Patch

:warning: Only in the unlikely and regretful event of a release happening with a bug.

- On your own fork, make a new branch `patch` based on `upstream/main` or `upstream/dev`.
- Fix the bug, and bump version (X.Y.Z+1).
- Open a pull-request from `patch` to `main`/`dev` with the changes.

## Pipeline contribution conventions

To make the `biodiversitycellatlas/bca_preprocessing` code and processing logic more understandable for new contributors and to ensure quality, we semi-standardise the way the code and other contributions are written.

### Adding a new step

If you wish to contribute a new step, please use the following coding standards:

1. Define the corresponding input channel into your new process from the expected previous process channel.
2. Write the process block (see below).
3. Define the output channel if needed (see below).
4. Add any new parameters to `nextflow.config` with a default (see below).
5. Add any new parameters to `nextflow_schema.json` with help text (via the `nf-core pipelines schema build` tool).
6. Add sanity checks and validation for all relevant parameters.
7. Perform local tests to validate that the new code works as expected.
8. If applicable, add a new test in the `tests` directory.
9. Update MultiQC config `assets/multiqc_config.yml` so relevant suffixes, file name clean up and module plots are in the appropriate order. If applicable, add a [MultiQC](https://https://multiqc.info/) module.
10. Add a description of the output files and if relevant any appropriate images from the MultiQC report to `docs/output.md`.

### Default values

Parameters should be initialised / defined with default values within the `params` scope in `nextflow.config`.

Once there, use `nf-core pipelines schema build` to add to `nextflow_schema.json`.

### Default processes resource requirements

Sensible defaults for process resource requirements (CPUs / memory / time) for a process should be defined in `conf/base.config`. These should generally be specified generic with `withLabel:` selectors so they can be shared across multiple processes/steps of the pipeline. A nf-core standard set of labels that should be followed where possible can be seen in the [nf-core pipeline template](https://github.com/nf-core/tools/blob/main/nf_core/pipeline-template/conf/base.config), which has the default process as a single core-process, and then different levels of multi-core configurations for increasingly large memory requirements defined with standardised labels.

The process resources can be passed on to the tool dynamically within the process with the `${task.cpus}` and `${task.memory}` variables in the `script:` block.

### Naming schemes

Please use the following naming schemes, to make it easy to understand what is going where.

- initial process channel: `ch_output_from_<process>`
- intermediate and terminal channels: `ch_<previousprocess>_for_<nextprocess>`

### Nextflow version bumping

If you are using a new feature from core Nextflow, you may bump the minimum required version of nextflow in the pipeline with: `nf-core pipelines bump-version --nextflow . [min-nf-version]`

### Images and figures

For overview images and other documents we follow the nf-core [style guidelines and examples](https://nf-co.re/developers/design_guidelines).


