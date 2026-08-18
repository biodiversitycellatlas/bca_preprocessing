# Protocol-specific steps & whitelists

---

## Table of Contents

1. [Parse Biosciences](#parse-biosciences)
2. [BD Rhapsody](#bd-rhapsody)
3. [MARS-seq](#mars-seq)
4. [sci-RNA-seq3](#sci-rna-seq3)

---

## Parse Biosciences

Protocol `"parse_biosciences_WT_mini"`:
    Barcode whitelists needed:
        - bc_data_n26_R1_v3_4
        - bc_data_v1
        - bc_data_R3_v4

Protocol `"parse_biosciences_WT"`:
    Barcode whitelists needed:
        - bc_data_n107_R1_v3_4
        - bc_data_v1
        - bc_data_R3_v4


The two Parse protocols describe the same chemistry and differ only in kit size, and therefore
only in which bc1 (RT) whitelist you pass:

| Protocol                    | bc1 whitelist            | bc2 / bc3 whitelists              |
| --------------------------- | ------------------------ | --------------------------------- |
| `parse_biosciences_WT_mini` | `bc_data_n26_R1_v3_4`    | `bc_data_v1` / `bc_data_R3_v4`    |
| `parse_biosciences_WT`      | `bc_data_n107_R1_v3_4`   | `bc_data_v1` / `bc_data_R3_v4`    |

The larger WT kit ships more RT barcodes than the WT Mini kit, so it can carry more samples and
more cells per run. Everything else is identical: both protocols share the same read layout and
the same `star_soloTypestring` and `alevin_*_geometry` values in
`conf/seqtech_parameters.config`, and both take the same route through the subworkflow. Only the
bc1 whitelist you list in `bc_whitelist` changes.

Parse Biosciences pools many samples into one library and separates them afterwards by the
first barcode round (bc1, the RT well barcode). The pipeline therefore splits each FASTQ pair
per sample *before* mapping, driven by the `rt` column of the samplesheet, which lists the
wells belonging to that sample (e.g. `A1-A6`).

**Samplesheet.** Fill the `rt` column with the group-well definition of the sample, and put the
cDNA read in `fastq_cDNA` (R1) and the barcode read in `fastq_BC_UMI` (R2):

> Wells are specified in blocks, ranges, or individually like this:<br/> >&nbsp;&nbsp;&nbsp;&nbsp;'A1:C6' specifies a block as [top-left]:[bottom-right]; A1-A6, B1-B6, C1-C6.<br/> >&nbsp;&nbsp;&nbsp;&nbsp;'A1-B6' specifies a range as [start]-[end]; A1-A12, B1-6.<br/> >&nbsp;&nbsp;&nbsp;&nbsp;'C4' specifies a single well.<br/>
> Multiple selections are joined by commas (no space), e.g. 'A1-A6,B1:D3,C4'

The BC/UMI read (`fastq_BC_UMI` column, R2) has a fixed v3 layout, with the three barcode
rounds separated by constant linkers:

| Positions (1-based) | Content                    |
| ------------------- | -------------------------- |
| 1-10                | UMI                        |
| 11-18               | bc3 (barcode round 3)      |
| 19-30               | linker `ATGAGGGGTCAG`      |
| 31-38               | bc2 (barcode round 2)      |
| 39-50               | linker `TCCAACCACCTC`      |
| 51-58               | bc1 (RT / well barcode)    |

`parse_workflow` (`subworkflows/local/pre-processing/preprocs_parse_biosciences.nf`) picks one
of three routes:

- **`PARSEBIO_CUSTOM_DEMUX`** (default, used when `splitpipe_demultiplex_script` is not set).
  `bin/parsebio_custom_demux.py` streams the read pair, slices bc1 out of positions 51-58 of
  the BC/UMI read, and matches it against `bc_whitelist_parse_splitwells` — a CSV of
  `bci,sequence,uid,well,stype` — allowing a Hamming distance of up to 2. A pair is kept when
  the well of its best-matching barcode falls inside the sample's `rt` well range. The step
  writes `<sample>_group_<sample>_R1/R2.fastq.gz` and reports how many reads matched at edit
  distance 0, 1 and 2, plus the unmatched and too-short counts.
- **`PARSEBIO_PIPELINE_DEMUX`**, when `splitpipe_demultiplex_script` points at Parse's own
  demultiplexing script, which performs the same split with the vendor's implementation.
- **No split at all**, when `perform_demultiplexing` is set to `false`; the reads are passed
  straight through.

Unlike MARS-seq, the reads themselves are never rewritten: demultiplexing only decides which
pairs end up in which sample's FASTQ, so the coordinates above apply unchanged to the original
read. The three whitelists in `bc_whitelist` must be given in the same order as
`--soloCBposition`, i.e. bc1, bc2, bc3.

If `splitpipe_installation` points at a split-pipe installation (and `splitpipe_conda_env` is
set), `PARSEBIO_PIPELINE_MKREF` and `PARSEBIO_PIPELINE` additionally run split-pipe
`-m all --chemistry v3` on the same split FASTQs, so the pipeline's results can be compared
against the commercial one.

## BD Rhapsody

Protocol `"bd_rhapsody_v1"`:
    Barcode whitelists needed:
        - A96_cell_key1
        - A96_cell_key2
        - A96_cell_key3

Protocol `"bd_rhapsody_enhancedbeads"`:
    Barcode whitelists needed:
        - A96_cell_key1
        - A96_cell_key2
        - A96_cell_key3

Protocol `"bd_rhapsody_enhancedbeads"`, working with V2/V3:
    Barcode whitelists needed:
        - B384_cell_key1
        - B384_cell_key2
        - B384_cell_key3


BD Rhapsody builds its cell label from three sections (CLS1, CLS2, CLS3) carried by the BC/UMI
read (`fastq_BC_UMI` column, R1), separated by constant linkers and followed by the UMI. The
two supported bead versions place them differently:

| Protocol                      | CLS1 | CLS2  | CLS3  | UMI   | Bases used |
| ----------------------------- | ---- | ----- | ----- | ----- | ---------- |
| `bd_rhapsody_v1`              | 1-9  | 22-30 | 44-52 | 53-60 | 1-60       |
| `bd_rhapsody_enhancedbeads`   | 1-9  | 14-22 | 27-35 | 36-43 | 1-43       |

(positions are 1-based; the gaps between the sections are fixed linker sequences.)

`bd_rhapsody_workflow` (`subworkflows/local/pre-processing/preprocs_bd_rhapsody.nf`) treats the
two versions differently:

- **`bd_rhapsody_v1`**: the bead layout is fixed, so there is nothing to rewrite and the reads
  are passed through untouched.
- **`bd_rhapsody_enhancedbeads`**: these reads start with 0-3 *variable* bases in front of
  CLS1, which would shift every downstream coordinate. `RM_VARBASES`
  (`modules/local/tools/cutadapt/main.nf`) removes them with three anchored 5' adapters,
  `^A`, `^GT` and `^TCA`; a read that carries no variable bases matches none of them and is
  left as it is. Only the BC/UMI read is trimmed — the cDNA read travels through cutadapt's
  paired output unchanged — and the result is written as
  `noVB_<sample>_R1/R2_001.fastq.gz`.

The positions in the table above, and the `star_soloTypestring` and `alevin_*_geometry` values
in `conf/seqtech_parameters.config`, therefore apply *after* variable-base removal. Both
protocols map with `--soloType CB_UMI_Complex` and `--soloCBmatchWLtype 1MM`, so the three
whitelists in `bc_whitelist` must be given in CLS1, CLS2, CLS3 order.

If `rhapsody_installation` points at a BD Rhapsody installation, `BDRHAP_PIPELINE_MKREF`,
`BDRHAP_PIPELINE_YAML` and `BDRHAP_PIPELINE` also run the commercial pipeline for comparison.
Note that these are deliberately given the **untrimmed** reads: the BD Rhapsody pipeline
detects and handles the variable bases itself.

## MARS-seq

Protocols `"marsseq_v1"` and `"marsseq_v2"`:
    Barcode whitelists: none. Both run with `--soloCBwhitelist None` and
    `--soloCBmatchWLtype Exact`; cells are called from the raw matrix by the
    pipeline's own cell-calling step (`cellfilter_method`).

MARS-seq splits its reads differently from every other supported protocol: read 1 carries
the batch barcode *and* the mapping region, read 2 carries the cell barcode and the UMI,
and both reads contain blocks that must be ignored.

**Samplesheet.** Read 1 goes in the `fastq_cDNA` column even though it is not pure cDNA, and
read 2 in `fastq_BC_UMI`; no other columns are needed. The pre-processing step below rebuilds
both reads before mapping, so the "wrong-looking" assignment is resolved by the pipeline.

| Version   | Read 1 (`fastq_cDNA` column) | Read 2 (`fastq_BC_UMI` column) |
| --------- | ---------------------------- | ------------------------------ |
| `marsseq_v1` | `5I 4B 45M 5I` (59 nt)    | `6W 4R 5I` (15 nt)             |
| `marsseq_v2` | `5I 4B 45M 5I` (59 nt)    | `7W 8R 5I` (20 nt)             |

Where `B` = batch barcode, `W` = cell barcode, `R` = UMI, `M` = mapping region and
`I` = ignore.

`MARSSEQ_BUILD_READS` (`subworkflows/local/pre-processing/preprocs_marsseq.nf`) rewrites
every read pair before mapping:

- the cDNA read becomes the 45 nt mapping region of read 1, with all ignore bases removed;
- the BC/UMI read becomes the batch barcode of read 1 followed by the cell barcode and UMI
  of read 2, with the trailing ignore bases removed.

The resulting BC/UMI read is therefore one contiguous block:

| Version   | BC/UMI read | Cell barcode (CB) | UMI            |
| --------- | ----------- | ----------------- | -------------- |
| `marsseq_v1` | 14 nt    | positions 1-10 (4 nt batch + 6 nt cell) | positions 11-14 |
| `marsseq_v2` | 19 nt    | positions 1-11 (4 nt batch + 7 nt cell) | positions 12-19 |

The batch barcode is part of the cell barcode on purpose: cells from different batches
pooled into one FASTQ then receive distinct barcodes.

The read layouts live in `conf/seqtech_parameters.config` as `marsseq_read1_design` and
`marsseq_read2_design`, next to the `star_soloTypestring` and `alevin_*_geometry` values
that describe the same coordinates, so the two cannot drift apart.

Note on `marsseq_v1`: its UMI is only 4 nt, allowing 256 distinct molecules per cell and
gene. UMI deduplication is therefore set to `Exact` (a 1-mismatch method would collapse
distinct molecules), and counts saturate at the high end. This is a property of the v1
chemistry, not of the pipeline.

## sci-RNA-seq3

Protocol `"sciRNAseq3"`:
    Barcode whitelists: `bc_whitelist` points at the sci-rocket barcode *definition* file, not at
    a plain list of sequences. The per-run p5, p7, ligation and RT whitelists are generated by
    the demultiplexing step and handed to the mapper automatically.

**Samplesheet.** sci-RNA-seq3 is the only protocol that needs the `p5`, `p7` and `rt` columns.
Fill them with the *index* of each barcode, not its sequence; the indexes are defined in the
sci-RNA-seq3 barcode whitelist file, where you can match a sequence to the barcode name used in
the samplesheet. The p5 and p7 sequences must also be present in the FASTQ headers — see
[this bcl2fastq script](../assets/bcl2fastq2.sh) for how to demultiplex the run so that they are.

`sciRNAseq3_nogather_workflow`
(`subworkflows/local/pre-processing/preprocs_sciRNAseq3_nogather.nf`) runs `SCIROCKET_DEMUX`,
which reads p5 and p7 out of the read header and the ligation, RT and UMI blocks out of the
barcode read, then writes a single synthetic 48 nt read per pair:

| Positions (1-based) | Content              |
| ------------------- | -------------------- |
| 1-10                | p7                   |
| 11-20               | p5                   |
| 21-30               | ligation barcode     |
| 31-40               | RT barcode           |
| 41-48               | UMI                  |

The ligation barcode comes in a 9 nt and a 10 nt flavour; the 9 nt version is padded with a
`G` so the block is always 10 nt and the downstream coordinates stay fixed. A read is only kept
when its p5 + p7 + RT combination matches a sample in the samplesheet, which is why those three
columns are required.

`FASTP` then trims adapters and low-quality bases. The `star_soloTypestring` and
`alevin_*_geometry` values in `conf/seqtech_parameters.config` address this synthetic read, not
the original one. Setting `perform_demultiplexing = false` skips the rebuild, which also leaves
the generated whitelists undefined.
