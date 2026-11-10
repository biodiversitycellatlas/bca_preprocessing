# Protocols & barcode whitelists

The `protocol` parameter selects the pre-processing steps and the read layout used for mapping. The layouts (`star_soloTypestring`, `alevin_*_geometry`) and the protocol-specific mapping settings live in [`conf/seqtech_parameters.config`](../conf/seqtech_parameters.config).



### Table of Contents

1. [Overview](#overview)
2. [Barcode whitelists](#barcode-whitelists)
   - [CRAM input](#cram-input)
3. [10x Genomics](#10x-genomics)
4. [OAK](#oak)
5. [Parse Biosciences](#parse-biosciences)
6. [BD Rhapsody](#bd-rhapsody)
7. [sci-RNA-seq3](#sci-rna-seq3)
8. [MARS-seq](#mars-seq)
9. [ScaleBio](#scalebio)
10. [seqspec](#seqspec)

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
| ScaleBio            | `scalebio_quantum_rna`, `scalebio_rna_v1_1`            | correct barcodes, split samples on RT (+PCR) wells, trim | `fastq_indices`, `rt`, `p5` (optional) |
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
| `scalebio_quantum_rna`        | bead 1, bead 2, bead 3, RT, PCR lists, see [ScaleBio](#scalebio)               |
| `scalebio_rna_v1_1`           | ligation, RT, PCR (i5) lists, see [ScaleBio](#scalebio)                        |
| `seqspec`                     | the whitelist(s) belonging to the seqspec layout                               |


### CRAM input

Instead of FASTQs, a samplesheet row can give its reads as unaligned CRAM (or BAM) files in a `cram` column, in place of `fastq_cDNA`, `fastq_BC_UMI` and `fastq_indices`. A glob is allowed, for example one CRAM per lane. A sample's rows must all use FASTQs or all use CRAMs.

```csv
sample,cram,expected_cells,p5,p7,rt
sampleA,run1/lib_S1_L00*.cram,5000,,,1A-2C
```

`cram_type` says what the files hold:

- **`"illumina"`** (the default, any protocol): both reads of every pair, with the index reads in the `BC` (sequence) and `QT` (quality) tags. Files written by `samtools import --i1 --i2`, or by Picard `FastqToSam`/`IlluminaBasecallsToSam`, are like this.
  - `CRAM_TO_FASTQ` converts each file with `samtools fastq` into bcl-convert-named FASTQs (`<file>_R1_001.fastq.gz`, `_R2_`, `_I1_`, `_I2_`). The protocol's `cdna_read` (`conf/seqtech_parameters.config`) says which read holds the cDNA.
  - The number of index reads is counted from the separator in `BC`, so `ACGT-TTGG` is two. Two indices joined without a separator need `cram_index_format`, for example `"i8i8"` (see `samtools fastq --index-format`).
  - A file with single-end reads stops the run, so that no reads are lost.
- **`"ultima"`** (ScaleBio only): Ultima Genomics reads, whose barcodes the Ultima trimmer moved into SAM tags. See [ScaleBio](#scalebio).

Unaligned CRAMs need no reference. For reference-compressed CRAMs, set `cram_reference` to the FASTA they were written against. `seqspec` does not take CRAMs.


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

**Demultiplexing.** By default, `PARSEBIO_CUSTOM_DEMUX` matches the RT barcode against the well whitelist given in `bc_whitelist_parse_splitwells` (a CSV with the columns `bci,sequence,uid,well,stype`, as shipped with split-pipe), at the smallest edit distance, at most 2.

Set `splitpipe_demultiplex_script` to use split-pipe's own demultiplexing script instead.

**Round-1 barcodes.** Each round-1 well holds an oligo-dT (`T`) and a random-hexamer (`R`) RT primer, and in some kits targeted (`X`) primers, each with its own barcode (the `stype` column of the CSV). split-pipe names cells by well, so these reads form one cell. A mapper given every round-1 barcode would split each cell into a T and an R barcode (in a typical run, about 40% and 60% of its reads). With `parse_collapse_bc1 = true` (the default), `PARSEBIO_COLLAPSE_BC1` therefore rewrites the round-1 barcode in the barcode read to its well's T barcode before mapping:

- the barcode is matched against `bc_whitelist_parse_splitwells` like in the demultiplexing (smallest edit distance, at most 2; no match with an `N`), so it is also corrected;
- reads whose best matches lie in different wells, unmatched reads and reads too short to hold the barcode keep their barcode;
- the mappers get the T barcodes of the CSV as round-1 whitelist instead of the first `bc_whitelist` file; the round-2 and round-3 files are used as given;
- split-pipe, if enabled, still receives the original reads.

The rewrite runs after the demultiplexing, whichever is used, and also with `perform_demultiplexing = false`.


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


## ScaleBio

Supports the QuantumScale RNA kit (`scalebio_quantum_rna`) and the 3-level ScaleRNA v1.1 / Extended Throughput kits (`scalebio_rna_v1_1`).

**Read layout.** A cell is identified by barcodes spread over several reads. Each one is corrected to its list on its own:

| | `scalebio_quantum_rna` | `scalebio_rna_v1_1` |
| --- | --- | --- |
| cDNA | read 1 | read 2 |
| UMI | read 2, nt 1-8 | read 1, 8 nt after the splint |
| RT barcode (sample) | read 2, nt 9-16 | read 1, 10 nt after the UMI |
| other cell barcodes | index 1: bead barcodes at nt 1-8, 13-20, 25-32 | read 1: ligation barcode (8-9 nt) before the splint `TCAGAGC` |
| PCR barcode | index 2 (i5), 8 nt; a well may hold a pool of sequences | index 2 (i5), 10 nt |
| cell barcode | bead 1 + bead 2 + bead 3 + RT + PCR | ligation + RT + PCR |

**Samplesheet.** Put the cDNA read in `fastq_cDNA`, the other read in `fastq_BC_UMI`, and the index reads in `fastq_indices`. QuantumScale needs I1 and I2, v1.1 needs I2. Index reads are recognised by `_I1_`/`_I2_` in their names, as written by bcl-convert with `CreateFastqForIndexReads,1`. The lanes of a read are paired by name, so name them the same apart from the read token, for example `_L001_R1_`, `_L001_I1_`.

A library pools many samples. Give every sample its own row with the **same** FASTQs. The pipeline reads a library once per lane, not once per sample.

- `rt`: the sample's RT wells, for example `1A-2H`, `3A;3B` or `1A;5A-5H`. Ranges of plate wells run column by column, so `1A-2C` is all of column 1 plus `2A-2C`. Separate the parts with `;`.
- `p5` (optional): the sample's PCR wells, as named in the PCR list, for example `QSR-1;QSR-2` or a range `QSR-1-QSR-6`. Left empty, the sample takes every PCR well.
- No RT × PCR well combination may belong to two samples.
- Extended Throughput: the i7 pool (plate) is not part of the cell barcode, so give each plate's FASTQs their own sample names, for example `sampleX_AP1` and `sampleX_AP2`.
- v1.1 reads only the i5 index. A single index read, given as `_I1_` or coming from a single-index CRAM, is taken as the i5.

**CRAM input.** Illumina-style CRAMs are converted to FASTQ first, see [CRAM input](#cram-input).

Ultima CRAMs (`cram_type = "ultima"`) are read as they are. Each file is one task, and a library's samples again share its files. The read is the cDNA, and the barcodes come from SAM tags:
- `scalebio_ultima_tags` names the tag of each barcode, as `segment=TAG[:START],...`. START is where the barcode begins in the tag's value, for tags that hold a constant linker first. QuantumScale defaults to the Ultima trimmer's "Scale RNA" format, `pcr=RG,bead1=b1,bead2=b2:4,bead3=b3:4,rt=b4`: `b2` and `b3` hold a constant 4 nt before the 8-nt bead barcode. There is no default for v1.1.
- A tag may hold:
  - a barcode sequence, corrected like any other;
  - a sequence longer than the barcode with no START given: the barcode is looked for at its end, then at its start;
  - a barcode id of the list, for example `QSR-1`, which is what `RG` holds;
  - an id or sequence within a longer value, for example a read group `run7-QSR-1`.
- `scalebio_ultima_umi` says where the UMI is. QuantumScale defaults to `tag:UM`, as in the "Scale RNA" format. It takes either:
  - `tag:XX`, the SAM tag that holds the UMI; or
  - `seq:START:LENGTH`, its place in the read, 0-based. It is cut from the cDNA.

  The CRAM header's `@CO` lines describe each tag. Check a read with `samtools view file.cram | head -1`.

**Barcode lists.** `bc_whitelist` takes the lists of the kit, in this order:
- QuantumScale: bead 1, bead 2, bead 3, RT, PCR.
- v1.1: ligation, RT, PCR.

Each list is a headerless text file with one barcode per line:
- an id, then its sequence, separated by tabs or spaces, for example `1A  ACGTACGT`;
- a pooled well lists all its sequences after the id;
- a list of bare sequences is also read, numbered by line.

Get the sequences from your kit's documentation or from ScaleBio. Check with ScaleBio whether your licence lets you reuse the lists that come with their own software.

**Demultiplexing.** First, `SCALEBIO_PREPARE_BARCODES` checks the lists: no sequence may appear twice, and every barcode must fit its segment. Then `SCALEBIO_DEMUX` processes each lane (or Ultima CRAM) of each library. For every read it:
1. finds the splint (v1.1 only; at most 1 mismatch);
2. corrects every barcode to its list. At most `scalebio_max_mismatches` substitutions are corrected (default 1). An `N` counts as a substitution, and a barcode equally close to two entries is dropped;
3. drops reads whose UMI contains an `N` or is a single repeated base;
4. trims the cDNA at a poly-A run (8 nt) and a poly-G run (12 nt, QuantumScale) or poly-T run (8 nt, v1.1), plus 3 or more of these bases left at its end, and drops reads shorter than 16 nt;
5. assigns the read to a sample by its RT and PCR wells.

The PCR index orientation on I2 is detected per lane from its first 200,000 reads (`scalebio_i5_revcomp = "auto"`), or set it with `"true"` or `"false"`.

Each sample's reads are written with a new barcode read, which the mapping layout in `conf/seqtech_parameters.config` refers to. It holds the cell barcode, then the UMI:
- **v1.1** (37 nt): the corrected barcodes in list order, with pooled PCR wells as their first sequence and ligation barcodes padded to 9 nt.
- **QuantumScale** (24 nt): its five barcodes together are 40 nt, longer than STARsolo's 31 nt maximum. They also make about 10<sup>9</sup> combinations, too many for its whitelists. So the cell barcode is the cell's combination of list entries, written as one number in 16 bases. It is lossless: `scalebio_demux.py decode --barcodes demultiplex/scalebio/barcodes/barcodes.json --cells barcodes.tsv` gives the bead, RT and PCR ids of each cell barcode.

`SCALEBIO_GATHER` joins a sample's lanes. Demultiplexing cannot be skipped: `perform_demultiplexing` has no effect for ScaleBio.

**Metrics.** Per lane, `demultiplex/scalebio/<library>/<chunk>.metrics.json` gives:
- reads per category: `Pass`, `LinkerError`, `BarcodeError`, `UmiFiltered`, `TooShort`, `Unassigned` (valid barcodes, but no sample);
- per barcode: `Exact`, `Corrected`, `Ambiguous`, `NoMatch`, `Error`;
- reads per sample and per RT well.

`demultiplex/scalebio/<sample>/<sample>.scalebio_demux.json` sums them for each sample. Set `scalebio_write_unassigned = true` to also keep the unassigned reads.

> [!NOTE]
> The mappers get no whitelist for QuantumScale. The barcodes are already corrected, so STARsolo takes the 16 nt cell barcode as it is (`--soloCBwhitelist None`), and alevin-fry uses its knee permit list. v1.1 passes the per-segment lists to both mappers.



## seqspec

For protocols not listed above, set `protocol = "seqspec"` and point `seqspec_file` to a [seqspec](https://github.com/pachterlab/seqspec) file describing the reads. STARsolo reads the barcode and UMI positions from this file. Reads are trimmed with fastp before mapping.
