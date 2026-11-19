# biodiversitycellatlas/bca_preprocessing: Documentation

The biodiversitycellatlas/bca_preprocessing documentation is split into the following pages:

- [Installation](../README.md#installation)
  - Cloning the pipeline and checking its dependencies.

- [Samplesheet](SAMPLESHEET.md)
  - Every samplesheet column, the columns per protocol, CRAM input, and how to combine samples.

- [Protocols & barcode whitelists](PROTOCOLS_AND_WHITELISTS.md)
  - What the pipeline does with each protocol's reads, and which barcode whitelists it needs.

- [Configuration parameters](CONFIGURATION_PARAMETERS.md)
  - All parameters of the configuration file.

- [Usage](usage.md)
  - How to run the pipeline, the profiles, the run methods and the two-pass filtering.

- [Output](output.md)
  - An overview of the various files & results produced by the pipeline and how to interpret them.

- [External pipelines](INSTALLATION_EXTERNAL_PIPELINES.md)
  - Installing Cell Ranger, split-pipe and BD Rhapsody to run alongside the pipeline.

- [Resource Tuning](RESOURCE_TUNING.md)
  - How to see what the cluster actually used per process, and how to adjust the scheduled resources for the next run.

- [Troubleshooting](TROUBLESHOOTING.md)
