# Protocol-specific steps & whitelists

---

## Table of Contents

1. [BD Rhapsody](#bd-rhapsody)
2. [Parse Biosciences](#parse-biosciences)
3. [MARS-seq](#mars-seq)

---

### Parse Biosciences

Protocol `"parse_biosciences_WT_mini"`:
    Barcode whitelists needed:
        - bc_data_n26_R1_v3_4
        - bc_data_v1
        - bc_data_R3_v4

Protocol `"parse_biosciences_WT_mini"`:
    Barcode whitelists needed:
        - bc_data_n107_R1_v3_4
        - bc_data_v1
        - bc_data_R3_v4


### BD Rhapsody

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


### MARS-seq

Protocols `"marsseq_v1"` and `"marsseq_v2"`:
    Barcode whitelists: none. Both run with `--soloCBwhitelist None` and
    `--soloCBmatchWLtype Exact`; cells are called from the raw matrix by the
    pipeline's own cell-calling step (`cellfilter_method`).

MARS-seq splits its reads differently from every other supported protocol: read 1 carries
the batch barcode *and* the mapping region, read 2 carries the cell barcode and the UMI,
and both reads contain blocks that must be ignored.

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
