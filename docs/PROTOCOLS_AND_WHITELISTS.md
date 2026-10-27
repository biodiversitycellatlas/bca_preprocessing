# Protocols & barcode whitelists

The `protocol` parameter selects the pre-processing steps and the read layout used for mapping. The layouts (`star_soloTypestring`, `alevin_*_geometry`) and the protocol-specific mapping settings live in [`conf/seqtech_parameters.config`](../conf/seqtech_parameters.config).



### Table of Contents

1. [Overview](#overview)
2. [Barcode whitelists](#barcode-whitelists)
3. [10x Genomics](#10x-genomics)
4. [OAK](#oak)
5. [Parse Biosciences](#parse-biosciences)
6. [BD Rhapsody](#bd-rhapsody)
7. [sci-RNA-seq3](#sci-rna-seq3)
8. [MARS-seq](#mars-seq)
9. [seqspec](#seqspec)

---

## Overview

| Platform            | `protocol`                                             | Pre-processing before mapping                   | Extra samplesheet columns            |
| ------------------- | ------------------------------------------------------ | ----------------------------------------------- | ------------------------------------ |
| 10x Genomics        | `10xv1`, `10xv2`, `10xv3`, `10xv4`                     | none                                            | -                                    |
| OAK                 | `oak_v1`                                               | split aliquots on their i5/i7 index             | `fastq_indices`, `p5`, `p7` (optional) |
| Parse Biosciences   | `parse_biosciences_WT_mini`, `parse_biosciences_WT`    | split samples on their RT wells                 | `rt`                                 |
| BD Rhapsody         | `bd_rhapsody_v1`, `bd_rhapsody_enhancedbeads`          | remove variable bases (Enhanced Beads only)     | -                                    |
| sci-RNA-seq3        | `sciRNAseq3`                                           | demultiplex with sci-rocket, trim with fastp    | `p5`, `p7`, `rt`                     |
| MARS-seq            | `marsseq_v1`, `marsseq_v2`                             | rebuild the reads                               | -                                    |
| Other               | `seqspec`                                              | trim with fastp                                 | -                                    |

Pre-processing steps that split or rebuild reads can be skipped with `perform_demultiplexing = false`, and fastp trimming with `perform_trimming = false`. Results of each step are written to `demultiplex/`, `fastp/` or `marsseq_reformat/`.

The commercial pipelines (Cell Ranger, split-pipe, BD Rhapsody) can run alongside for comparison, see [supported external pipelines](INSTALLATION_EXTERNAL_PIPELINES.md).


## Barcode whitelists

`bc_whitelist` takes one or more paths **or** URLs (not mixed), separated by whitespace. URLs are downloaded and unzipped automatically.

- Each file is a plain list with one barcode per line.
- Protocols with a split cell barcode (Parse, BD Rhapsody) take one file per barcode segment, **in the order listed below**.

| `protocol`                    | `bc_whitelist`                                                                 |
| ----------------------------- | ------------------------------------------------------------------------------ |
| `10xv1` - `10xv4`, `oak_v1`   | the 10x whitelist matching the chemistry                                       |
| `parse_biosciences_WT_mini`   | `bc_data_n26_R1_v3_4` `bc_data_v1` `bc_data_R3_v4`                             |
| `parse_biosciences_WT`        | `bc_data_n107_R1_v3_4` `bc_data_v1` `bc_data_R3_v4`                            |
| `bd_rhapsody_v1`              | `A96_cell_key1` `A96_cell_key2` `A96_cell_key3`                                |
| `bd_rhapsody_enhancedbeads`   | `A96_cell_key1` `A96_cell_key2` `A96_cell_key3` (V1 beads), or `B384_cell_key1` `B384_cell_key2` `B384_cell_key3` (V2/V3 beads) |
| `sciRNAseq3`                  | the sci-rocket barcode definition file, see [sci-RNA-seq3](#sci-rna-seq3)      |
| `marsseq_v1`, `marsseq_v2`    | none, leave empty                                                              |
| `seqspec`                     | the whitelist(s) belonging to the seqspec layout                               |


## 10x Genomics

Put read 2 in `fastq_cDNA` and read 1 (cell barcode + UMI) in `fastq_BC_UMI`. No pre-processing is needed. Set `perform_cellranger = true` to also run Cell Ranger.


## OAK

OAK splits an overloaded 10x 3' run into aliquots, each with its own i5 index (and optionally its own i7 index). Because the same cell barcode can occur in several aliquots, every aliquot is demultiplexed and processed as a separate sample.

**Samplesheet.** Add one row per aliquot, with a unique `sample` name. All rows point to the same FASTQs:

| Column          | Content                                                 |
| --------------- | ------------------------------------------------------- |
| `fastq_cDNA`    | read 2                                                  |
| `fastq_BC_UMI`  | read 1                                                  |
| `fastq_indices` | both index reads, e.g. `/path/Undetermined_S0_I*`       |
| `p5`            | i5 index sequence of the aliquot (**required**)         |
| `p7`            | i7 index sequence of the aliquot (optional)             |

`OAK_DEMUX` selects the reads of each aliquot on its i7 (if given) and i5 index, allowing `oak_demux_max_mismatches` mismatches (default `1`). The read counts after each step are written to `demultiplex/<sample>/<sample>_oak_demux.log`. Each aliquot is then mapped as a 10x 3' v3.1 sample.


## Parse Biosciences

The two protocols share the same chemistry and differ only in kit size, and therefore in the first-round (RT) whitelist, see [Barcode whitelists](#barcode-whitelists).

Parse pools several samples into one library, which are told apart by the RT wells they were loaded into. The pipeline therefore splits the FASTQs per sample before mapping.

**Samplesheet.** Put read 1 in `fastq_cDNA` and read 2 (barcodes + UMI) in `fastq_BC_UMI`. Fill the `rt` column with the wells of the sample:

- `A1:C6` is a block, from top-left to bottom-right: A1-A6, B1-B6, C1-C6;
- `A1-B6` is a range, from start to end: A1-A12, B1-B6;
- `C4` is a single well;
- selections are joined by commas without spaces, e.g. `A1-A6,B1:D3,C4`.

**Demultiplexing.** By default, `PARSEBIO_CUSTOM_DEMUX` matches the RT barcode against the well whitelist given in `bc_whitelist_parse_splitwells` (a CSV with the columns `bci,sequence,uid,well,stype`, as shipped with split-pipe), allowing up to 2 mismatches. Set `splitpipe_demultiplex_script` to use split-pipe's own demultiplexing script instead.


## BD Rhapsody

**Samplesheet.** Put read 2 in `fastq_cDNA` and read 1 (cell label + UMI) in `fastq_BC_UMI`.

Enhanced Beads reads start with 0-3 variable bases before the cell label. For `bd_rhapsody_enhancedbeads`, `RM_VARBASES` removes these with cutadapt, so the barcode positions are fixed during mapping. `bd_rhapsody_v1` reads are used as they are.

The BD Rhapsody pipeline, if enabled, receives the untrimmed reads, as it handles the variable bases itself.


## sci-RNA-seq3

**Samplesheet.** Put read 2 in `fastq_cDNA` and read 1 in `fastq_BC_UMI`. Fill the `p5`, `p7` and `rt` columns with the barcode **names** (not sequences) from the sci-rocket barcode definition file, e.g. `p5 = A02`, `p7 = D09`, `rt = P02-A01:P02-H07`. The p5 and p7 sequences must be present in the FASTQ headers, see [this bcl2fastq script](../assets/bcl2fastq2.sh).

**Whitelist.** `bc_whitelist` points to the sci-rocket barcode definition file, not to a list of sequences. The per-segment whitelists used for mapping are generated during demultiplexing.

**Pre-processing.** `SCIROCKET_DEMUX` keeps only reads whose p5, p7 and RT barcodes match a sample in the samplesheet, and writes a new 48 nt barcode read: p7, p5, ligation and RT barcode (10 nt each, 9 nt ligation barcodes are padded with a `G`), followed by the 8 nt UMI. The mapping layout in `conf/seqtech_parameters.config` refers to this new read. `FASTP` then trims adapters and low-quality bases.

> [!NOTE]
> With `perform_demultiplexing = false`, no whitelists are generated for mapping.


## MARS-seq

**Samplesheet.** Put read 1 in `fastq_cDNA` and read 2 in `fastq_BC_UMI`. No whitelist is used: cells are called from the raw matrix by the pipeline's cell-calling step (`cellfilter_method`).

**Pre-processing.** Read 1 carries the batch barcode and the cDNA, read 2 the cell barcode and UMI. `MARSSEQ_BUILD_READS` rebuilds them before mapping:

- **cDNA read**: the 45 nt mapping region of read 1;
- **barcode read**: batch barcode + cell barcode + UMI. The batch barcode is part of the cell barcode, so cells from different batches in the same FASTQ stay distinct.

| `protocol`   | Read 1 design | Read 2 design | Cell barcode | UMI   |
| ------------ | ------------- | ------------- | ------------ | ----- |
| `marsseq_v1` | `5I4B45M5I`   | `6W4R5I`      | 10 nt        | 4 nt  |
| `marsseq_v2` | `5I4B45M5I`   | `7W8R5I`      | 11 nt        | 8 nt  |

`B` = batch barcode, `W` = cell barcode, `R` = UMI, `M` = mapping region, `I` = ignored.



## seqspec

For protocols not listed above, set `protocol = "seqspec"` and point `seqspec_file` to a [seqspec](https://github.com/pachterlab/seqspec) file describing the reads. STARsolo reads the barcode and UMI positions from this file. Reads are trimmed with fastp before mapping.
