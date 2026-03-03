# biodiversitycellatlas/bca_preprocessing: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-02-24

- Added pavianCore as submodule and increased resources by @bonitavw in #41
- Set limitGenomeGenerateRAM and limitBAMsortRAM as custom config variables by @bonitavw in #42
- Converted to Seqera containers & Fixed permissions by @bonitavw in #43
- MAJOR refractoring of the workflows and subworkflows & added PavianCore by @bonitavw in #44
- Added possibility of specifying kit for split-pipe by @bonitavw in #45
- Remove strict conda dependency by @bonitavw in #46
- Increase resources withLabel process_high_memory by @bonitavw in #47
- Configuration updates and improvements by @bonitavw in #48
- Major refactoring of input and emit channels and fixed env GeneExt by @bonitavw in #49
- Added html dashboard generation and reporting workflow by @bonitavw in #50
- Modify emit channels preprocs subworkflows by @bonitavw in #51

## [0.1.2] - 2025-12-18

- Fix boolean expression to disable generation of BAM during STARsolo by @bonitavw in #36
- Fix STARsolo silent error because of missing pigz by @bonitavw in #37
- Fix empty MAPREADS variable by writing to file + fix seqera container by @bonitavw in #38
- Swithed to cellbender container exclusively by @bonitavw in #39
- Removed conda dependencies from geneext module by @bonitavw in #40

## [0.1.1] - 2025-12-16

- docs: Add troubleshooting guide by @bonitavw in #28
- Allow disabling BAM generation via star_generateBAM config by @bonitavw in #29
- Expanded well-range support (e.g. WT kit) for parse demultiplexing by @bonitavw in #30
- Update README.md & bump manifest version by @bonitavw in #31
- Added Seqera containers to custom modules by @bonitavw in #32
- Split 10x_saturate saturation_table process into two modules by @bonitavw in #33
- Added 10x v4 to seqtech parameters by @bonitavw in #34
- Missing out specification in subworkflow by @bonitavw in #35

## [0.1.0] - 2025-12-11

- Initial release.

[unreleased]: https://github.com/biodiversitycellatlas/bca_preprocessing/v0.2.0...HEAD
[0.2.0]: https://github.com/biodiversitycellatlas/bca_preprocessing/compare/v0.1.2...v0.2.0
[0.1.2]: https://github.com/biodiversitycellatlas/bca_preprocessing/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/biodiversitycellatlas/bca_preprocessing/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/biodiversitycellatlas/bca_preprocessing/releases/tag/v0.1.0
