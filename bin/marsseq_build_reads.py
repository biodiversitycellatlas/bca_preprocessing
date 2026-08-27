#!/usr/bin/env python3
"""
Rebuild MARS-seq read pairs into the cDNA / BC+UMI split the pipeline expects.

Both reads are described by a fixed-layout design string, e.g. "5I4B45M5I", where each
segment is a length followed by a one-letter code:

    B = batch barcode
    W = cell barcode
    R = UMI
    I = ignore
    M = mapping region

Two FASTQs are written per sample:

    <sample>_marsseq_cDNA.fastq.gz     the M segments of read 1
    <sample>_marsseq_BC_UMI.fastq.gz   the B segments of read 1, followed by every
                                       non-ignore segment of read 2 (W then R)

The batch barcode therefore ends up directly in front of the cell barcode, so batch,
cell barcode and UMI form one contiguous block that STARsolo can address with
CB_UMI_Simple and alevin with a single bc/umi geometry.

Input reading auto-detects gzip vs. plain text, and also
prefers an external decompressor (pigz/gzip) over the stdlib gzip module.
"""

import argparse
import gzip
import io
import os
import re
import shutil
import subprocess
import sys
from typing import List, Tuple

# Segment codes that carry information (everything except "ignore")
SEGMENT_CODES = "BIMRW"

# How many FASTQ records to buffer in memory before flushing to the compressor.
# Tune via --batch-size if needed; higher = fewer write() syscalls, more RAM.
DEFAULT_BATCH_SIZE = 50_000

# Default output compression level. gzip/pigz default is 6; for intermediate
# files that are read once immediately afterwards, a lower level trades a
# larger file for much faster writing.
DEFAULT_COMPRESS_LEVEL = 4


# --------------------------- Utilities ---------------------------


def parse_design(design: str) -> List[Tuple[int, str]]:
    """Parse a design string such as '5I4B45M5I' into [(5,'I'), (4,'B'), (45,'M'), (5,'I')]."""
    design = design.strip().upper().replace(" ", "")
    if not design:
        raise ValueError("Empty read design string.")

    segments = [(int(length), code) for length, code in re.findall(r"(\d+)([A-Z])", design)]

    # Re-assembling the segments must reproduce the input, otherwise something was dropped
    if "".join(f"{length}{code}" for length, code in segments) != design:
        raise ValueError(f"Cannot parse read design '{design}'. Expected e.g. '5I4B45M5I'.")

    for length, code in segments:
        if code not in SEGMENT_CODES:
            raise ValueError(f"Unknown segment code '{code}' in design '{design}'. Use one of: {SEGMENT_CODES}.")
        if length < 1:
            raise ValueError(f"Segment '{length}{code}' in design '{design}' has no length.")

    return segments


def slice_ranges(segments: List[Tuple[int, str]], codes: str) -> List[Tuple[int, int]]:
    """Return the [start, end) offsets of every segment whose code is in 'codes', in read order."""
    ranges = []
    offset = 0
    for length, code in segments:
        if code in codes:
            ranges.append((offset, offset + length))
        offset += length
    return ranges


def design_length(segments: List[Tuple[int, str]]) -> int:
    """Total number of bases the design describes."""
    return sum(length for length, _ in segments)


def take(sequence: str, ranges: List[Tuple[int, int]]) -> str:
    """Concatenate the given ranges of a sequence or quality string."""
    return "".join(sequence[start:end] for start, end in ranges)


def open_fastq(path: str):
    """Open a FASTQ for reading as text, gzipped or not (detected from the magic bytes)."""
    with open(path, "rb") as probe:
        is_gzipped = probe.read(2) == b"\x1f\x8b"
    if is_gzipped:
        # Use an external decompressor when available
        decompressor = shutil.which("pigz") or shutil.which("gzip")
        if decompressor:
            proc = subprocess.Popen(
                [decompressor, "-dc", path],
                stdout=subprocess.PIPE,
                bufsize=1024 * 1024,
            )
            return io.TextIOWrapper(proc.stdout, encoding="utf-8")
        return gzip.open(path, "rt")
    return open(path, "rt", buffering=1024 * 1024)


class _SubprocessWriter:
    """
    Minimal writable text-mode wrapper around a subprocess's stdin pipe, so
    calling code can just do `.write(str)` and `.close()` regardless of
    whether compression is done externally or via the stdlib fallback.
    """

    def __init__(self, proc: subprocess.Popen, out_fh):
        self._proc = proc
        self._out_fh = out_fh
        self._stdin = io.TextIOWrapper(proc.stdin, encoding="utf-8")

    def write(self, data: str) -> None:
        self._stdin.write(data)

    def close(self) -> None:
        self._stdin.close()
        self._proc.wait()
        self._out_fh.close()
        if self._proc.returncode not in (0, None):
            raise RuntimeError(f"Compressor subprocess exited with code {self._proc.returncode}.")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()


def open_gzip_writer(path: str, compress_level: int, threads: int = 1):
    """
    Open a gzip-compressed output stream backed by an external `pigz`/`gzip`
    process instead of Python's zlib wrapper. Falls back to the stdlib gzip
    module if neither binary is on PATH.

    Returns a writable text-mode file-like object. Caller must close() it,
    which also waits for the subprocess to finish and flush its output file.
    """
    pigz_path = shutil.which("pigz")
    gzip_path = shutil.which("gzip")

    out_fh = open(path, "wb")

    if pigz_path:
        cmd = [pigz_path, "-p", str(max(1, threads)), f"-{compress_level}", "-c"]
        proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=out_fh, bufsize=1024 * 1024)
        return _SubprocessWriter(proc, out_fh)
    if gzip_path:
        cmd = [gzip_path, f"-{compress_level}", "-c"]
        proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=out_fh, bufsize=1024 * 1024)
        return _SubprocessWriter(proc, out_fh)

    # No external compressor available: fall back to the stdlib, still buffered.
    out_fh.close()
    return gzip.open(path, "wt", compresslevel=compress_level)


def read_fastq(handle):
    """Yield (header, sequence, plus, quality) tuples from an open FASTQ handle."""
    while True:
        header = handle.readline()
        if not header:
            return
        sequence = handle.readline().rstrip("\n")
        plus = handle.readline().rstrip("\n")
        quality = handle.readline().rstrip("\n")
        if not quality:
            raise ValueError("Truncated FASTQ record: the file does not end on a complete record.")
        yield header.rstrip("\n"), sequence, plus, quality


def read_id(header: str) -> str:
    """Read name without the '@' and without the trailing comment."""
    return header[1:].split(" ", 1)[0].split("\t", 1)[0]


# --------------------------- Main ---------------------------


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample_id", required=True, help="Sample ID, used to name the output files.")
    parser.add_argument("--fq1", required=True, help="Read 1 FASTQ (batch barcode + mapping region).")
    parser.add_argument("--fq2", required=True, help="Read 2 FASTQ (cell barcode + UMI).")
    parser.add_argument("--read1-design", required=True, help="Read 1 layout, e.g. '5I4B45M5I'.")
    parser.add_argument("--read2-design", required=True, help="Read 2 layout, e.g. '6W4R5I'.")
    parser.add_argument("--output", default=".", help="Output directory (default: current directory).")
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE,
        help=f"Number of FASTQ records to buffer per write batch (default: {DEFAULT_BATCH_SIZE}).")
    parser.add_argument("--compress-level", type=int, default=DEFAULT_COMPRESS_LEVEL, choices=range(1, 10), metavar="[1-9]",
        help=f"gzip/pigz compression level for output files (default: {DEFAULT_COMPRESS_LEVEL}; lower is faster, higher is smaller).")
    parser.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", os.environ.get("NXF_TASK_CPUS", 1)) or 1),
        help="Threads to hand to pigz for output compression (default: from Slurm/Nextflow env or 1).")
    args = parser.parse_args()

    read1_segments = parse_design(args.read1_design)
    read2_segments = parse_design(args.read2_design)

    # cDNA read: the mapping region of read 1
    cdna_ranges = slice_ranges(read1_segments, "M")

    # BC/UMI read: batch barcode from read 1, then everything informative from read 2
    batch_ranges = slice_ranges(read1_segments, "B")
    bc_umi_ranges = slice_ranges(read2_segments, "WR")

    if not cdna_ranges:
        raise ValueError(f"Read 1 design '{args.read1_design}' has no mapping region (M).")
    if not bc_umi_ranges:
        raise ValueError(f"Read 2 design '{args.read2_design}' has no cell barcode (W) or UMI (R).")

    read1_length = design_length(read1_segments)
    read2_length = design_length(read2_segments)
    cdna_length = sum(end - start for start, end in cdna_ranges)
    bc_umi_length = sum(end - start for start, end in batch_ranges + bc_umi_ranges)

    os.makedirs(args.output, exist_ok=True)
    cdna_path = os.path.join(args.output, f"{args.sample_id}_marsseq_cDNA.fastq.gz")
    bc_umi_path = os.path.join(args.output, f"{args.sample_id}_marsseq_BC_UMI.fastq.gz")
    log_path = os.path.join(args.output, f"{args.sample_id}_marsseq_reformat.log")

    total = kept = too_short = 0

    in1 = open_fastq(args.fq1)
    in2 = open_fastq(args.fq2)
    out_cdna = open_gzip_writer(cdna_path, args.compress_level, args.threads)
    out_bc_umi = open_gzip_writer(bc_umi_path, args.compress_level, args.threads)

    try:
        reads1 = read_fastq(in1)
        reads2 = read_fastq(in2)

        cdna_buffer: List[str] = []
        bc_umi_buffer: List[str] = []
        batch_size = max(1, args.batch_size)

        for record1 in reads1:
            record2 = next(reads2, None)
            if record2 is None:
                raise ValueError(f"{args.fq2} ended before {args.fq1}: the FASTQs are not paired.")

            header1, sequence1, _, quality1 = record1
            header2, sequence2, _, quality2 = record2
            total += 1

            if read_id(header1) != read_id(header2):
                raise ValueError(
                    f"Read pair {total} is out of sync: '{read_id(header1)}' vs '{read_id(header2)}'."
                )

            # Reads shorter than their design cannot be sliced reliably, so drop the pair
            if len(sequence1) < read1_length or len(sequence2) < read2_length:
                too_short += 1
                continue

            cdna_buffer.append(
                f"{header1}\n{take(sequence1, cdna_ranges)}\n+\n{take(quality1, cdna_ranges)}\n"
            )
            bc_umi_buffer.append(
                "{}\n{}{}\n+\n{}{}\n".format(
                    header2,
                    take(sequence1, batch_ranges), take(sequence2, bc_umi_ranges),
                    take(quality1, batch_ranges), take(quality2, bc_umi_ranges),
                )
            )
            kept += 1

            if len(cdna_buffer) >= batch_size:
                out_cdna.write("".join(cdna_buffer))
                out_bc_umi.write("".join(bc_umi_buffer))
                cdna_buffer.clear()
                bc_umi_buffer.clear()

        if cdna_buffer:
            out_cdna.write("".join(cdna_buffer))
            out_bc_umi.write("".join(bc_umi_buffer))

        if next(reads2, None) is not None:
            raise ValueError(f"{args.fq1} ended before {args.fq2}: the FASTQs are not paired.")
    finally:
        out_cdna.close()
        out_bc_umi.close()
        in1.close()
        in2.close()

    summary = [
        f"sample_id\t{args.sample_id}",
        f"read1_design\t{args.read1_design}\t({read1_length} nt expected)",
        f"read2_design\t{args.read2_design}\t({read2_length} nt expected)",
        f"cDNA_read_length\t{cdna_length}",
        f"BC_UMI_read_length\t{bc_umi_length}",
        f"read_pairs_total\t{total}",
        f"read_pairs_written\t{kept}",
        f"read_pairs_skipped_too_short\t{too_short}",
        f"compress_level\t{args.compress_level}",
        f"batch_size\t{args.batch_size}",
        f"threads\t{args.threads}",
    ]
    with open(log_path, "w") as log:
        log.write("\n".join(summary) + "\n")
    print("\n".join(summary))

    if kept == 0:
        print(f"ERROR: no read pairs survived reformatting for '{args.sample_id}'.", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
