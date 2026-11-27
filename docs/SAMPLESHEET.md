# Samplesheet

The samplesheet lists the samples of a run, where their reads are, and how they are combined. It is a comma-separated file with a header row, given to the pipeline as `input`.

### Table of Contents

1. [Overview](#overview)
2. [Columns](#columns)
3. [Columns per protocol](#columns-per-protocol)
4. [CRAM input](#cram-input)
5. [Several rows for one sample](#several-rows-for-one-sample)
6. [Combining samples](#combining-samples)
7. [Examples](#examples)

---

## Overview

```csv
sample,fastq_cDNA,fastq_BC_UMI,expected_cells
spis_10x,data/spis_10x_R2.fastq.gz,data/spis_10x_R1.fastq.gz,5000
```

- One row per sample, or per lane of a sample, see [Several rows for one sample](#several-rows-for-one-sample).
- One samplesheet per species and reference. All its samples are mapped to the same `ref_fasta` and `ref_gtf`.
- Columns can be in any order, and columns a protocol does not use can be left out or left empty.
- Relative paths are read from the samplesheet's own folder. Globs and `~/` are allowed. See [How paths are resolved](CONFIGURATION_PARAMETERS.md#how-paths-are-resolved).
- The samplesheet is checked against [`assets/schema_input.json`](../assets/schema_input.json) when the run starts (`validate_params = true`, the default), and a copy is written to `<outdir>/pipeline_info/`.

More samplesheets are in [`examples/`](../examples/README.md) and under [Examples](#examples).


## Columns

| Column | Required | Content |
|--------|----------|---------|
| `sample` | Yes | Sample name, without spaces. Rows with the same name are merged before mapping, see [Several rows for one sample](#several-rows-for-one-sample). |
| `fastq_cDNA` | Yes, unless `cram` | The FASTQ holding the cDNA, ending in `.fastq.gz`, `.fq.gz`, `.fastq` or `.fq`. Which read this is depends on the protocol, see [Columns per protocol](#columns-per-protocol). |
| `fastq_BC_UMI` | Yes, unless `cram` | The FASTQ holding the cell barcode and UMI. |
| `fastq_indices` | OAK, ScaleBio | The index read FASTQs. A glob such as `Undetermined_S0_I*` gives both. |
| `cram` | No | Unaligned CRAM or BAM files, in place of the three FASTQ columns. See [CRAM input](#cram-input). |
| `expected_cells` | Yes | A positive whole number: the cells you expect. STARsolo's cell filter and CellSweep use it, and the second-derivative cell calling uses it only to break ties. Both reports show it next to the cells called. See [Cell calling](CONFIGURATION_PARAMETERS.md#cell-calling). |
| `manual_cutoff` | With `cellfilter_method = "manual_cutoff"` | The UMI threshold of the sample's cells. With that method every sample needs one, or the run stops. With any other method the column is ignored, with a warning. |
| `p5` | OAK, sci-RNA-seq3, ScaleBio | Depends on the protocol, see [Columns per protocol](#columns-per-protocol). |
| `p7` | OAK, sci-RNA-seq3 | Depends on the protocol, see [Columns per protocol](#columns-per-protocol). |
| `rt` | Parse, sci-RNA-seq3, ScaleBio; optional for SPLiT-seq | The sample's reverse-transcription wells or barcodes, in the protocol's notation. |
| `qc_report` | No | A group name, without spaces. Samples with the same group are shown together in `dashboard.html`. See [Combining samples](#combining-samples). |
| `filtering_report` | No | A group name, without spaces. Samples with the same group are pooled into one Metacell2 run and one entry of `filtering_report.html`. See [Combining samples](#combining-samples). |


## Columns per protocol

| `protocol` | `fastq_cDNA` | `fastq_BC_UMI` | `fastq_indices` | `p5` | `p7` | `rt` |
|------------|--------------|----------------|-----------------|------|------|------|
| `10xv2` - `10xv4` | R2 | R1 | - | - | - | - |
| `oak_v1` | R2 | R1 | I1 and I2 | i5 sequence | i7 sequence (optional) | - |
| `parse_biosciences_*` | R1 | R2 | - | - | - | RT wells |
| `splitseq_v1` | R1 | R2 | - | - | - | RT wells (optional) |
| `bd_rhapsody_*` | R2 | R1 | - | - | - | - |
| `sciRNAseq3` | R2 | R1 | - | p5 name | p7 name | RT names |
| `marsseq_v1`, `marsseq_v2` | R1 | R2 | - | - | - | - |
| `scalebio_quantum_rna` | R1 | R2 | I1 and I2 | PCR wells (optional) | - | RT wells |
| `scalebio_rna_v1_1` | R2 | R1 | I2 | PCR wells (optional) | - | RT wells |
| `seqspec` | the read with the `cdna` region | the read with the barcode and `umi` regions | - | - | - | - |

What the pipeline does with these columns is described per protocol in [Protocols & barcode whitelists](PROTOCOLS_AND_WHITELISTS.md).


### 10x Genomics, BD Rhapsody, MARS-seq

Only the two FASTQ columns, as in the table.


### OAK

One row per aliquot, each with its own `sample` name. All rows of a channel point to the same FASTQs; the aliquots are told apart by their index sequences:

- `fastq_indices`: both index reads, e.g. `/path/Undetermined_S0_I*`;
- `p5`: the aliquot's i5 index sequence (required);
- `p7`: the aliquot's i7 index sequence (optional).


### Parse Biosciences

`rt` lists the RT wells the sample was loaded into:

- `A1:C6` is a block, from top-left to bottom-right: A1-A6, B1-B6, C1-C6;
- `A1-B6` is a range, from start to end: A1-A12, B1-B6;
- `C4` is a single well;
- selections are joined by commas without spaces. The commas would split the CSV field, so put the value in double quotes: `"A1-A6,B1:D3,C4"`.

Every sublibrary reuses the same barcode combinations for other cells. Give a sample one row per sublibrary, each with its own `sample` name and the same `rt` wells, and combine them with [`qc_report` and `filtering_report`](#combining-samples).


### SPLiT-seq

Read 2 must be at least 94 nt, see [SPLiT-seq](PROTOCOLS_AND_WHITELISTS.md#split-seq). Give each sublibrary (i7 index) its own row and its own `sample` name; the same barcodes in two sublibraries are two different cells.

`rt` is optional. Left empty, the row takes every read of its sublibrary. To split a sublibrary into the samples loaded into different RT wells, give one row per sample, all pointing to the sublibrary's FASTQs, with the sample's wells in `rt`, in the [Parse notation](#parse-biosciences) on the 48 RT wells `A1`-`D12`:

```csv
sample,fastq_cDNA,fastq_BC_UMI,rt
sub1_larva,sub1_R1.fastq.gz,sub1_R2.fastq.gz,A1-B12
sub1_adult,sub1_R1.fastq.gz,sub1_R2.fastq.gz,C1-D12
sub2_larva,sub2_R1.fastq.gz,sub2_R2.fastq.gz,A1-B12
sub2_adult,sub2_R1.fastq.gz,sub2_R2.fastq.gz,C1-D12
```

Combine the sublibraries of a sample with [`qc_report` and `filtering_report`](#combining-samples).


### sci-RNA-seq3

`p5`, `p7` and `rt` hold barcode **names** from the sci-rocket barcode definition file (`bc_whitelist`), not sequences, e.g. `p5 = A02`, `p7 = D09`, `rt = P02-A01:P02-H07`. The p5 and p7 sequences must be in the FASTQ headers, see [this bcl2fastq script](../assets/bcl2fastq2.sh).


### ScaleBio

A library pools many samples. Give every sample its own row with the **same** FASTQs; the pipeline reads a library once per lane, not once per sample.

- `fastq_indices`: QuantumScale needs I1 and I2, v1.1 needs I2. Index reads are recognised by `_I1_`/`_I2_` in their names, as written by bcl-convert with `CreateFastqForIndexReads,1`. The lanes of a read are paired by name, so name them the same apart from the read token, e.g. `_L001_R1_` and `_L001_I1_`. v1.1 takes a single index read, given as `_I1_` or from a single-index CRAM, as the i5.
- `rt`: the sample's RT wells, e.g. `1A-2H`, `3A;3B` or `1A;5A-5H`. Ranges run column by column, so `1A-2C` is all of column 1 plus `2A-2C`. Parts are separated by `;`.
- `p5` (optional): the sample's PCR wells, as named in the PCR list, e.g. `QSR-1;QSR-2` or a range `QSR-1-QSR-6`. Left empty, the sample takes every PCR well.
- No RT × PCR well combination may belong to two samples.
- Extended Throughput: the i7 pool (plate) is not part of the cell barcode, so give each plate's FASTQs their own sample names, e.g. `sampleX_AP1` and `sampleX_AP2`, and combine them with [`qc_report` and `filtering_report`](#combining-samples).


### seqspec

The two FASTQ columns follow the seqspec file: the read that covers the barcode and UMI regions goes in `fastq_BC_UMI`, the other read that covers the cDNA in `fastq_cDNA`. See [seqspec](PROTOCOLS_AND_WHITELISTS.md#seqspec).


## CRAM input

Instead of FASTQs, a row can give its reads as unaligned CRAM (or BAM) files in a `cram` column, in place of `fastq_cDNA`, `fastq_BC_UMI` and `fastq_indices`. A glob is allowed, for example one CRAM per lane. A sample's rows must all use FASTQs or all use CRAMs.

```csv
sample,cram,expected_cells,p5,p7,rt
sampleA,run1/lib_S1_L00*.cram,5000,,,1A-2C
```

`cram_type` says what the files hold:

- **`"illumina"`** (the default, any protocol): both reads of every pair, with the index reads in the `BC` (sequence) and `QT` (quality) tags. Files written by `samtools import --i1 --i2`, or by Picard `FastqToSam`/`IlluminaBasecallsToSam`, are like this.
  - `CRAM_TO_FASTQ` converts each file with `samtools fastq` into bcl-convert-named FASTQs (`<file>_R1_001.fastq.gz`, `_R2_`, `_I1_`, `_I2_`). The protocol's `cdna_read` (`conf/seqtech_parameters.config`) says which read holds the cDNA.
  - The number of index reads is counted from the separator in `BC`, so `ACGT-TTGG` is two. Two indices joined without a separator need `cram_index_format`, for example `"i8i8"` (see `samtools fastq --index-format`).
  - A file with single-end reads stops the run, so that no reads are lost.
- **`"ultima"`** (ScaleBio only): Ultima Genomics reads, whose barcodes the Ultima trimmer moved into SAM tags. See [ScaleBio: Ultima CRAMs](PROTOCOLS_AND_WHITELISTS.md#ultima-crams).

Unaligned CRAMs need no reference. For reference-compressed CRAMs, set `cram_reference` to the FASTA they were written against. `seqspec` does not take CRAMs.


## Several rows for one sample

Rows with the same `sample` name are merged before mapping: their reads are mapped as one sample. Use this for rows that hold **the same cells**, such as the lanes of a library or a top-up sequencing of it.

- `expected_cells`, `manual_cutoff`, `p5`, `p7` and `rt` are taken from the sample's first row.
- `qc_report` and `filtering_report` must be the same in every row of a sample.
- A sample's rows must all use FASTQs or all use CRAMs.

Rows that hold **different cells** must not share a name, see [Combining samples](#combining-samples).


## Combining samples

Samplesheet rows can be combined at three points. Which one is right depends on whether the rows hold the same cells.

| Mechanism | Combined | Use it when the rows hold |
|-----------|----------|---------------------------|
| The same `sample` name on several rows | The FASTQs, before mapping | **The same cells**: several lanes, or a top-up sequencing of one library. |
| `qc_report` column | Only the numbers in `dashboard.html` | **Different cells** whose quality you want to judge together. |
| `filtering_report` column | The cell-called objects, just before Metacell2 | **Different cells** you want to filter and analyse as one set of metacells. |

Rows that hold different cells must not share a `sample` name. Their barcodes overlap, so merging their reads before mapping turns two cells with the same barcode into one, and mixes their empty droplets. This happens with:
- one sample loaded on several 10x or BD channels;
- the sublibraries of a Parse, SPLiT-seq, ScaleBio or sci-RNA-seq3 experiment, which reuse the same barcode combinations;
- the plates of a ScaleBio Extended Throughput run;
- the aliquots of an overloaded 10x channel (OAK).

Give such rows their own `sample` names and combine them with the two columns instead. Both are optional. An empty value leaves a sample on its own, and the two columns are independent of each other:

```csv
sample,fastq_cDNA,fastq_BC_UMI,expected_cells,rt,qc_report,filtering_report
nvec_well1_larva_rep1,...,...,2500,A1-A6,nvec_well1,nvec_larva
nvec_well1_larva_rep2,...,...,2500,A1-A6,nvec_well1,nvec_larva
nvec_well2_larva,...,...,2500,B1-B6,nvec_well2,nvec_larva
nvec_well3_adult,...,...,2500,C1-C6,nvec_well3,nvec_adult
nvec_well4_adult,...,...,2500,D1-D6,nvec_well4,nvec_adult
```

Here the two replicates of well 1 are shown together in the dashboard. The metacells and the filtering are done per stage: wells 1 and 2 as `nvec_larva`, wells 3 and 4 as `nvec_adult`.

The `filtering_report` groups need `perform_metacells = true`; without it they are ignored, with a warning. A group's outputs are named like a sample's, `<group><suffix>` (e.g. `nvec_larva_starsolo`), from `metacells/` onwards, see [output](output.md).


### Which mechanism to use

| Situation | `qc_report` | `filtering_report` |
|-----------|-------------|--------------------|
| One sample on several 10x or BD channels | the sample | the sample, or its condition |
| Sublibraries of one Parse, SPLiT-seq, ScaleBio or sci-RNA-seq3 experiment, one row per sublibrary | the sample | the sample |
| ScaleBio Extended Throughput plates (`sampleX_AP1`, `sampleX_AP2`) | `sampleX` | `sampleX` |
| Aliquots of one overloaded 10x channel (OAK) | the channel | the sample or condition |
| The wells of a sample split into rows, to see the quality of each well | the well or replicate | the condition |
| Biological replicates of one condition | the condition, to compare the replicates | the condition, for more cells per Metacell2 run and one filtering applied to all |
| Stages or tissues meant to share metacells (an atlas) | - | one group for all of them |

Do not group:
- **Two sequencing runs of the same library under `filtering_report`.** The same cells would enter Metacell2 twice. Give the runs one `sample` name instead. The pooling step warns when two samples share more than half of their cell barcodes. Unrelated captures share far fewer, a few percent at most.
- **Samples that need very different filtering.** The report's thresholds act on the pooled metacells, so one setting applies to the whole group.
- **Samples of different species or references.** Each needs its own samplesheet anyway, and the pooling refuses samples counted against different genes.

Pooled metacells can be made of one sample only, from a batch effect or from a cell type only that sample has. The report's *Main sample %* column shows this per metacell.


### Rules

Checked when the run starts:
- every row of a sample gives the same groups;
- a group cannot have the name of a sample, or end in an analytical-run suffix such as `_starsolo`;
- the name of a sample in a `filtering_report` group cannot contain `:`.

A group of one sample is allowed, and runs under the group's name.


### Changing the groups of a finished run

The groups are read from the samplesheet when the run starts, and are not part of any step before Metacell2. So the groups can be changed after a run, for example after looking at the per-sample dashboard. Edit the two columns and rerun with `-resume`. Only these steps run again:
- the pooling;
- Metacell2 for the groups that changed, and for samples taken out of a group;
- the two reports.

A [`post_mapping`](usage.md#re-running-after-mapping-post_mapping) run reads the groups from the current samplesheet as well.

A metacell selection made before the groups changed no longer applies to the changed groups. Export a new one from the new report, see [Filtering UMI matrices](usage.md#filtering-umi-matrices-two-passes).


## Examples

**10x, one library sequenced on two lanes.** The rows hold the same cells, so they share a name and are mapped as one sample.

```csv
sample,fastq_cDNA,fastq_BC_UMI,expected_cells
spis_adult,fastq/spis_S1_L001_R2_001.fastq.gz,fastq/spis_S1_L001_R1_001.fastq.gz,8000
spis_adult,fastq/spis_S1_L002_R2_001.fastq.gz,fastq/spis_S1_L002_R1_001.fastq.gz,8000
```

**Parse, two samples over two sublibraries.** Each sample has one row per sublibrary, with the same `rt` wells, and the groups put the sublibraries back together after mapping.

```csv
sample,fastq_cDNA,fastq_BC_UMI,expected_cells,rt,qc_report,filtering_report
larva_sl1,sl1/R1.fastq.gz,sl1/R2.fastq.gz,5000,A1-A6,larva,larva
larva_sl2,sl2/R1.fastq.gz,sl2/R2.fastq.gz,5000,A1-A6,larva,larva
adult_sl1,sl1/R1.fastq.gz,sl1/R2.fastq.gz,5000,"B1:C6,D1",adult,adult
adult_sl2,sl2/R1.fastq.gz,sl2/R2.fastq.gz,5000,"B1:C6,D1",adult,adult
```

**OAK, one channel split into three aliquots.** Every row points to the same FASTQs; the i5 index tells the aliquots apart.

```csv
sample,fastq_cDNA,fastq_BC_UMI,fastq_indices,expected_cells,p5,qc_report,filtering_report
ch1_aliquot1,run/Undetermined_S0_R2_001.fastq.gz,run/Undetermined_S0_R1_001.fastq.gz,run/Undetermined_S0_I*,10000,ACGTACGTAC,ch1,embryo
ch1_aliquot2,run/Undetermined_S0_R2_001.fastq.gz,run/Undetermined_S0_R1_001.fastq.gz,run/Undetermined_S0_I*,10000,TGCATGCATG,ch1,embryo
ch1_aliquot3,run/Undetermined_S0_R2_001.fastq.gz,run/Undetermined_S0_R1_001.fastq.gz,run/Undetermined_S0_I*,10000,GATCGATCGA,ch1,embryo
```

**sci-RNA-seq3, two samples on different RT plates.**

```csv
sample,fastq_cDNA,fastq_BC_UMI,expected_cells,p5,p7,rt
sample1,fastq/run_R2.fastq.gz,fastq/run_R1.fastq.gz,2500,A02,D09,P01-A01:P01-H12
sample2,fastq/run_R2.fastq.gz,fastq/run_R1.fastq.gz,2500,A02,D09,P02-A01:P02-H12
```

**ScaleBio QuantumScale, two samples in one library.** Both rows give the same FASTQs; the RT wells split them, and the second sample is limited to two PCR wells.

```csv
sample,fastq_cDNA,fastq_BC_UMI,fastq_indices,expected_cells,rt,p5
liver,lib/QSR_S1_L001_R1_001.fastq.gz,lib/QSR_S1_L001_R2_001.fastq.gz,lib/QSR_S1_L001_I*,20000,1A-6H,
kidney,lib/QSR_S1_L001_R1_001.fastq.gz,lib/QSR_S1_L001_R2_001.fastq.gz,lib/QSR_S1_L001_I*,20000,7A-12H,QSR-1;QSR-2
```

**Manual cutoffs for a `post_mapping` rerun.** With `cellfilter_method = "manual_cutoff"`, every sample needs a `manual_cutoff`. The FASTQs no longer have to exist for a `post_mapping` run, see [Re-running after mapping](usage.md#re-running-after-mapping-post_mapping).

```csv
sample,fastq_cDNA,fastq_BC_UMI,expected_cells,manual_cutoff
spis_adult,fastq/spis_R2.fastq.gz,fastq/spis_R1.fastq.gz,8000,450
spis_larva,fastq/spis_larva_R2.fastq.gz,fastq/spis_larva_R1.fastq.gz,8000,300
```
