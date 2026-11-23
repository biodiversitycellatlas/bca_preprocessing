# Protocols & barcode whitelists

The `protocol` parameter selects the pre-processing steps and the read layout used for mapping. The layouts (`star_soloTypestring`, `alevin_*_geometry`) and the protocol-specific mapping settings live in [`conf/seqtech_parameters.config`](../conf/seqtech_parameters.config).

Which read goes in which samplesheet column, and the format of the `p5`, `p7` and `rt` columns, are described in [Samplesheet](SAMPLESHEET.md#columns-per-protocol).

### Table of Contents

1. [Overview](#overview)
2. [Barcode whitelists](#barcode-whitelists)
3. [10x Genomics](#10x-genomics)
4. [OAK](#oak)
5. [Parse Biosciences](#parse-biosciences)
6. [BD Rhapsody](#bd-rhapsody)
7. [sci-RNA-seq3](#sci-rna-seq3)
8. [MARS-seq](#mars-seq)
9. [ScaleBio](#scalebio)
10. [seqspec](#seqspec)
11. [Known limitations](#known-limitations)

---

## Overview

| Platform            | `protocol`                                             | Pre-processing before mapping                   | [Samplesheet columns](SAMPLESHEET.md#columns-per-protocol) besides the FASTQs |
| ------------------- | ------------------------------------------------------ | ----------------------------------------------- | ------------------------------------ |
| 10x Genomics        | `10xv2`, `10xv3`, `10xv4`                              | none                                            | -                                    |
| OAK                 | `oak_v1`                                               | split aliquots on their i5/i7 index             | `fastq_indices`, `p5`, `p7` (optional) |
| Parse Biosciences   | `parse_biosciences_WT_mini`, `parse_biosciences_WT`    | split samples on their RT wells, collapse round-1 barcodes | `rt`                       |
| BD Rhapsody         | `bd_rhapsody_v1`, `bd_rhapsody_enhancedbeads`          | remove variable bases (Enhanced Beads only)     | -                                    |
| sci-RNA-seq3        | `sciRNAseq3`                                           | demultiplex with sci-rocket, trim with fastp    | `p5`, `p7`, `rt`                     |
| MARS-seq            | `marsseq_v1`, `marsseq_v2`                             | rebuild the reads                               | -                                    |
| ScaleBio            | `scalebio_quantum_rna`, `scalebio_rna_v1_1`            | correct barcodes, split samples on RT (+PCR) wells, trim | `fastq_indices`, `rt`, `p5` (optional) |
| Other               | `seqspec`                                              | none; fastp with `perform_trimming = true`      | -                                    |

- Pre-processing steps that split or rebuild reads can be skipped with `perform_demultiplexing = false` (not for sci-RNA-seq3, which needs them, nor for ScaleBio, where it has no effect), and fastp trimming with `perform_trimming = false`. Results of each step are written to `demultiplex/`, `fastp/` or `marsseq_reformat/`.
- Protocols that pool several samples or sublibraries reuse the same barcodes for different cells. Give each its own `sample` name and combine them after mapping, see [Combining samples](SAMPLESHEET.md#combining-samples).
- Reads can also be given as unaligned CRAMs, see [CRAM input](SAMPLESHEET.md#cram-input).
- The commercial pipelines (Cell Ranger, split-pipe, BD Rhapsody) can run alongside for comparison, see [supported external pipelines](INSTALLATION_EXTERNAL_PIPELINES.md).


## Barcode whitelists

`bc_whitelist` takes one or more paths **or** URLs (not mixed), separated by whitespace. URLs are downloaded and unzipped automatically.

- Each file is a plain list with one barcode per line, except for sci-RNA-seq3 and ScaleBio (see their sections).
- Protocols with a split cell barcode take one file per barcode segment, **in the order listed below**.

| `protocol`                    | `bc_whitelist`                                                                 |
| ----------------------------- | ------------------------------------------------------------------------------ |
| `10xv2` - `10xv4`, `oak_v1`   | the 10x whitelist matching the chemistry                                       |
| `parse_biosciences_WT_mini`   | `bc_data_n26_R1_v3_4` `bc_data_v1` `bc_data_R3_v4`                             |
| `parse_biosciences_WT`        | `bc_data_n107_R1_v3_4` `bc_data_v1` `bc_data_R3_v4`                            |
| `bd_rhapsody_v1`              | `A96_cell_key1` `A96_cell_key2` `A96_cell_key3`                                |
| `bd_rhapsody_enhancedbeads`   | `A96_cell_key1` `A96_cell_key2` `A96_cell_key3` (V1 beads), or `B384_cell_key1` `B384_cell_key2` `B384_cell_key3` (V2/V3 beads) |
| `sciRNAseq3`                  | the sci-rocket barcode definition file, see [sci-RNA-seq3](#sci-rna-seq3)      |
| `marsseq_v1`, `marsseq_v2`    | none, leave empty                                                              |
| `scalebio_quantum_rna`        | bead 1, bead 2, bead 3, RT, PCR lists, see [ScaleBio](#barcode-lists)          |
| `scalebio_rna_v1_1`           | ligation, RT, PCR (i5) lists, see [ScaleBio](#barcode-lists)                   |
| `seqspec`                     | `bc_whitelist`, or else the onlists of the seqspec file's barcode regions      |


## 10x Genomics

**Samplesheet.** Read 2 in `fastq_cDNA`, read 1 (cell barcode + UMI) in `fastq_BC_UMI`.

No pre-processing is needed. Set `perform_cellranger = true` to also run Cell Ranger.


## OAK

OAK splits an overloaded 10x 3' run into aliquots, each with its own i5 index (and optionally its own i7 index). Because the same cell barcode can occur in several aliquots, every aliquot is demultiplexed and processed as a separate sample.

**Samplesheet.** One row per aliquot, all pointing to the same FASTQs, with the aliquot's index sequences in `p5` (required) and `p7` (optional), see [OAK samplesheet](SAMPLESHEET.md#oak).

**Pre-processing.** `OAK_DEMUX` selects the reads of each aliquot on its i7 (if given) and i5 index, allowing `oak_demux_max_mismatches` mismatches (default `1`). The read counts after each step are written to `demultiplex/<sample>/<sample>_oak_demux.log`. Each aliquot is then mapped as a 10x 3' v3.1 sample.


## Parse Biosciences

The two protocols share the same chemistry and differ only in kit size, and therefore in the first-round (RT) whitelist, see [Barcode whitelists](#barcode-whitelists).

Parse pools several samples into one library, which are told apart by the RT wells they were loaded into. The pipeline therefore splits the FASTQs per sample before mapping.

**Samplesheet.** Read 1 in `fastq_cDNA`, read 2 (barcodes + UMI) in `fastq_BC_UMI`, and the sample's wells in `rt`. One row per sublibrary, see [Parse samplesheet](SAMPLESHEET.md#parse-biosciences).

**Demultiplexing.** By default, `PARSEBIO_CUSTOM_DEMUX` matches the RT barcode against the well whitelist given in `bc_whitelist_parse_splitwells` (a CSV with the columns `bci,sequence,uid,well,stype`, as shipped with split-pipe), at the smallest edit distance, at most 2. Set `splitpipe_demultiplex_script` to use split-pipe's own demultiplexing script instead.

**Round-1 barcodes.** Each round-1 well holds an oligo-dT (`T`) and a random-hexamer (`R`) RT primer, and in some kits targeted (`X`) primers, each with its own barcode (the `stype` column of the CSV). split-pipe names cells by well, so these reads form one cell. A mapper given every round-1 barcode would split each cell into a T and an R barcode (in a typical run, about 40% and 60% of its reads). With `parse_collapse_bc1 = true` (the default), `PARSEBIO_COLLAPSE_BC1` therefore rewrites the round-1 barcode in the barcode read to its well's T barcode before mapping:

- the barcode is matched against `bc_whitelist_parse_splitwells` like in the demultiplexing (smallest edit distance, at most 2; no match with an `N`), so it is also corrected;
- reads whose best matches lie in different wells, unmatched reads and reads too short to hold the barcode keep their barcode;
- the mappers get the T barcodes of the CSV as round-1 whitelist instead of the first `bc_whitelist` file; the round-2 and round-3 files are used as given;
- split-pipe, if enabled, still receives the original reads.

The rewrite runs after the demultiplexing, whichever is used, and also with `perform_demultiplexing = false`.


## BD Rhapsody

**Samplesheet.** Read 2 in `fastq_cDNA`, read 1 (cell label + UMI) in `fastq_BC_UMI`.

**Pre-processing.** Enhanced Beads reads start with 0-3 variable bases before the cell label. For `bd_rhapsody_enhancedbeads`, `RM_VARBASES` removes these with cutadapt, so the barcode positions are fixed during mapping. `bd_rhapsody_v1` reads are used as they are.

The BD Rhapsody pipeline, if enabled, receives the untrimmed reads, as it handles the variable bases itself.


## sci-RNA-seq3

**Samplesheet.** Read 2 in `fastq_cDNA`, read 1 in `fastq_BC_UMI`, and the barcode **names** (not sequences) in `p5`, `p7` and `rt`, see [sci-RNA-seq3 samplesheet](SAMPLESHEET.md#sci-rna-seq3).

**Whitelist.** `bc_whitelist` points to the sci-rocket barcode definition file, not to a list of sequences. The per-segment whitelists used for mapping are generated during demultiplexing, one set per sample: the sample's own p5, p7 and RT barcodes and every ligation barcode. Each sample is mapped and checked against its own set.

**Pre-processing.** `SCIROCKET_DEMUX` keeps only reads whose p5, p7 and RT barcodes match a sample in the samplesheet, and writes a new 48 nt barcode read: p7, p5, ligation and RT barcode (10 nt each, 9 nt ligation barcodes are padded with a `G`), followed by the 8 nt UMI. The mapping layout in `conf/seqtech_parameters.config` refers to this new read. `FASTP` then trims adapters and low-quality bases.

> [!NOTE]
> Demultiplexing cannot be skipped: `perform_demultiplexing = false` stops the run, as without it there is neither the 48 nt barcode read nor a whitelist to map against.


## MARS-seq

**Samplesheet.** Read 1 in `fastq_cDNA`, read 2 in `fastq_BC_UMI`.

**Whitelist.** None: cells are called from the raw matrix by the pipeline's cell-calling step (`cellfilter_method`).

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

**Samplesheet.** The cDNA read in `fastq_cDNA`, the other read in `fastq_BC_UMI`, the index reads in `fastq_indices`, the RT wells in `rt` and optionally the PCR wells in `p5`. The samples of a library share its FASTQs, see [ScaleBio samplesheet](SAMPLESHEET.md#scalebio).

#### Read layout

A cell is identified by barcodes spread over several reads. Each one is corrected to its list on its own:

| | `scalebio_quantum_rna` | `scalebio_rna_v1_1` |
| --- | --- | --- |
| cDNA | read 1 | read 2 |
| UMI | read 2, nt 1-8 | read 1, 8 nt after the splint |
| RT barcode (sample) | read 2, nt 9-16 | read 1, 10 nt after the UMI |
| other cell barcodes | index 1: bead barcodes at nt 1-8, 13-20, 25-32 | read 1: ligation barcode (8-9 nt) before the splint `TCAGAGC` |
| PCR barcode | index 2 (i5), 8 nt; a well may hold a pool of sequences | index 2 (i5), 10 nt |
| cell barcode | bead 1 + bead 2 + bead 3 + RT + PCR | ligation + RT + PCR |

#### Barcode lists

`bc_whitelist` takes the lists of the kit, in this order:
- QuantumScale: bead 1, bead 2, bead 3, RT, PCR.
- v1.1: ligation, RT, PCR.

Each list is a headerless text file with one barcode per line:
- an id, then its sequence, separated by tabs or spaces, for example `1A  ACGTACGT`;
- a pooled well lists all its sequences after the id;
- a list of bare sequences is also read, numbered by line.

Get the sequences from your kit's documentation or from ScaleBio. Check with ScaleBio whether your licence lets you reuse the lists that come with their own software.

#### Demultiplexing

First, `SCALEBIO_PREPARE_BARCODES` checks the lists: no sequence may appear twice, and every barcode must fit its segment. Then `SCALEBIO_DEMUX` processes each lane (or Ultima CRAM) of each library. For every read it:
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

> [!NOTE]
> The mappers get no whitelist for QuantumScale. The barcodes are already corrected, so STARsolo takes the 16 nt cell barcode as it is (`--soloCBwhitelist None`), and alevin-fry uses its knee permit list. v1.1 passes the per-segment lists to both mappers.

#### Metrics

Per lane, `demultiplex/scalebio/<library>/<chunk>.metrics.json` gives:
- reads per category: `Pass`, `LinkerError`, `BarcodeError`, `UmiFiltered`, `TooShort`, `Unassigned` (valid barcodes, but no sample);
- per barcode: `Exact`, `Corrected`, `Ambiguous`, `NoMatch`, `Error`;
- reads per sample and per RT well.

`demultiplex/scalebio/<sample>/<sample>.scalebio_demux.json` sums them for each sample. Set `scalebio_write_unassigned = true` to also keep the unassigned reads.

#### Ultima CRAMs

Illumina-style CRAMs are converted to FASTQ first, like for any protocol, see [CRAM input](SAMPLESHEET.md#cram-input). Ultima CRAMs (`cram_type = "ultima"`) are read as they are. Each file is one task, and a library's samples again share its files. The read is the cDNA, and the barcodes come from SAM tags:
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

#### Mapping

These settings follow ScaleBio's own ScaleRna pipeline:
- STARsolo counts `GeneFull_Ex50pAS`;
- it distributes multi-gene reads with `star_soloMultiMappers = "PropUnique"`, so the analysed matrix is `raw/UniqueAndMult-PropUnique.mtx` (see [STARsolo multi-gene reads](CONFIGURATION_PARAMETERS.md#starsolo-multi-gene-reads));
- a read may map to at most 6 loci, where STAR's default is 10 (`--outFilterMultimapNmax 6`, the protocol's `star_extraargs`). Reads at more loci are unmapped, as *mapped to too many loci*. For another limit, give the option with that value in your own `star_extraargs`;
- there is no minimum alignment score (`star_outFilterScoreMin = 0`; most other protocols use 30). A read still has to align over two thirds of its length, as STAR requires by default;
- samples read from Ultima CRAMs are counted on the reverse strand, because Ultima reads are antisense to the gene. Set `star_soloStrand` for another strand.

Two settings are those of the other protocols, and differ from ScaleRna's:

| Setting | This pipeline | ScaleRna | To follow ScaleRna |
|---------|---------------|----------|--------------------|
| UMI filtering and deduplication | `MultiGeneUMI_CR`, `1MM_CR` (Cell Ranger's rules) | `-`, `1MM_All` (STAR's defaults) | `star_soloUMIfiltering = "-"` and `star_soloUMIdedup = "1MM_All"`, together |
| Cell calling | `cellfilter_method`, by default the second derivative | its own cellFinder; STAR's filter is off | not available |

On one QuantumScale library (*C. elegans*, Ultima, 29 million reads) the pipeline, then still with `"EM"`, 10 loci and a minimum score of 30, counted 0.2% more reads on a single gene and 0.8% fewer UMIs than ScaleRna. Neither figure depends on the multi-gene method.


## seqspec

For protocols not listed above, set `protocol = "seqspec"` and point `seqspec_file` to a [seqspec](https://github.com/pachterlab/seqspec) file describing the reads. A worked example on the 10x test data is in [examples/seqspec](../examples/seqspec/conf).

**Samplesheet.** The read with the barcode and UMI regions in `fastq_BC_UMI`, the other read with a cDNA region in `fastq_cDNA`. `seqspec` does not take CRAMs.

#### What is read from the file

`SEQSPEC_LAYOUT` (`bin/seqspec_layout.py`) reads the spec's `rna` modality once per run, before any index is built, and writes the layout to `<outdir>/seqspec/seqspec_layout.env`:

- The **barcode read** is the read that covers the `barcode` region(s) and the `umi` region. The **cDNA read** is the other read that covers a `cdna` region.
- **STARsolo**:
  - one barcode region gives `--soloType CB_UMI_Simple` with its start and length, plus `--soloBarcodeReadLength 0`, because fastp may shorten barcode reads;
  - several barcode regions (a split barcode, e.g. between linkers) give `--soloType CB_UMI_Complex` with one `--soloCBposition` per region and `--soloUMIposition`. The whitelist matching then becomes `1MM`, because the protocol's `1MM_multi_Nbase_pseudocounts` works only for `CB_UMI_Simple`.
- **alevin-fry**: `alevin_bc_geometry` lists the barcode regions in the barcode read (e.g. `1[1-9,14-22,27-35]`), and `alevin_umi_geometry` gives the UMI. `alevin_read_geometry` starts where the cDNA starts in the cDNA read, which is usually `2[1-end]`.
- **Whitelists**: if `bc_whitelist` is unset, each barcode region's onlist is its whitelist, in read order. Remote onlists are downloaded and decompressed. Local ones are read relative to the seqspec file's folder, as seqspec does. If a barcode region has no onlist, set `bc_whitelist`. A single barcode without any onlist runs without a whitelist, which also skips the barcode check.

`star_soloTypestring`, `star_soloCBmatchWLtype` and the `alevin_*_geometry` parameters still win when you set them.

#### Limits

The run stops with a message when:
- the barcodes and the UMI are not all in one read;
- there is no UMI, or more than one UMI region;
- no other read covers the cDNA;
- a region up to the end of the barcodes and UMI has a variable length (`min_len` ≠ `max_len`), which leaves no fixed positions.

`seqspec check` runs too, but its findings are only printed. Only the `rna` modality is read.

seqspec's own `seqspec index -t starsolo` is not used. It prints `--soloType` itself, which the pipeline also added, so STAR got the parameter twice and stopped (*duplicate parameter "soloType"*). It also writes the barcode's end position as its length, which is wrong when the barcode does not start the read, and it has no form for split barcodes.

#### Trimming

Off by default for `seqspec`: the reads are mapped as they are, and STARsolo's own `CellRanger4` clipping removes the TSO and poly-A. Set `perform_trimming = true` to trim with fastp first. fastp keeps the pairs in step, so it drops a whole pair when either read fails its quality filter, the barcode read included. On the demo it dropped 42 of 2,375 pairs and trimmed no adapter from either read. Without trimming, the demo's counts are identical to those of `protocol = "10xv3"` on the same reads. The STARsolo settings are those of 10x v3 (`CellRanger4` adapter clipping, `EmptyDrops_CR`); change them through the `star_*` parameters if your protocol differs.


## Known limitations

**Homopolymer UMIs.** A UMI made of one repeated base (`GGGGGGGG`) is a sequencing or synthesis artefact, not a molecule tag, and STARsolo is meant to drop it. With the split-barcode layouts (`--soloType CB_UMI_Complex`: Parse, BD Rhapsody, sci-RNA-seq3, ScaleBio v1.1), STAR 2.7.11b drops only all-A UMIs ([alexdobin/STAR#2715](https://github.com/alexdobin/STAR/issues/2715), item 1g). All-C, all-G and all-T UMIs are counted, adding at most one UMI per cell and gene each.
- ScaleBio is not affected: `SCALEBIO_DEMUX` already drops these reads.
- In the Parse test data, 0.9 % of reads have an all-G UMI, but 99 % of those also lack a valid cell barcode; the counts change by less than 0.01 %. Runs with many poly-G reads, and BD Rhapsody, whose UMI ends the barcode read, may be affected more.
