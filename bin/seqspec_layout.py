#!/usr/bin/env python3
"""Barcode/UMI layout of a seqspec file for STARsolo and alevin-fry.

Reads the spec with the seqspec library, finds the read that holds the cell barcode(s)
and the UMI and the read that holds the cDNA, and writes:

  seqspec_layout.env  shell-sourceable settings:
      SEQSPEC_STAR_TYPE     --soloType ... (CB_UMI_Simple for one barcode region,
                            CB_UMI_Complex with --soloCBposition for several)
      SEQSPEC_STAR_CBMATCH  --soloCBmatchWLtype the layout needs ('1MM' for
                            CB_UMI_Complex, empty when the protocol's setting works)
      SEQSPEC_ALEVIN_BC / _UMI / _READ   salmon alevin geometries, read 1 being the
                            barcode read and read 2 the cDNA read, as ALEVIN_FRY passes them
      SEQSPEC_BARCODE_READ / SEQSPEC_CDNA_READ / SEQSPEC_N_BARCODES   for the logs
  onlists.tsv         one row per barcode region in read order: index, region_id,
                      urltype, url (empty when the region has no onlist)
  onlist_<n>.txt      with --fetch-onlists, the remote onlists, downloaded and
                      decompressed (STAR reads only plain-text whitelists). Local
                      onlists are left to the caller, as they are read relative to
                      the spec file's folder.

seqspec's own 'index -t starsolo' is not used: it writes --soloCBlen as the barcode's
end rather than its length, keeps only the first barcode region and has no
CB_UMI_Complex form.
"""

import argparse
import shlex
import subprocess
import sys
from pathlib import Path


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def region_type(region):
    rtype = region.region_type
    return str(getattr(rtype, "value", rtype)).lower()


def read_layout(spec, modality, read):
    """The read's atomic regions as (region, start, stop), 0-based half-open, cut to the read."""
    from seqspec.Region import itx_read, project_regions_to_coordinates
    from seqspec.utils import map_read_id_to_regions

    _, regions = map_read_id_to_regions(spec, modality, read.read_id)
    coords = itx_read(project_regions_to_coordinates(regions), 0, read.max_len)
    by_id = {r.region_id: r for r in regions}
    return [(by_id.get(c.region_id, c), c.start, c.stop) for c in coords]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--spec", required=True, type=Path, help="seqspec YAML file")
    parser.add_argument("--modality", default="rna", help="modality to read (default: rna)")
    parser.add_argument("--fetch-onlists", action="store_true", help="download the remote onlists")
    parser.add_argument("--out-env", default="seqspec_layout.env", type=Path)
    parser.add_argument("--out-onlists", default="onlists.tsv", type=Path)
    parser.add_argument("--onlist-dir", default=".", type=Path)
    parser.add_argument("--no-check", action="store_true", help="skip 'seqspec check'")
    args = parser.parse_args()

    try:
        from seqspec.utils import load_spec, read_remote_list
    except ImportError as err:
        fail(f"the seqspec Python package is not available ({err})")

    if not args.spec.is_file():
        fail(f"seqspec file '{args.spec}' does not exist")

    # seqspec's own validation, reported but not enforced: it also flags things the
    # pipeline does not use (e.g. FASTQ files listed in the spec that are not present)
    if not args.no_check:
        check = subprocess.run(["seqspec", "check", str(args.spec)], capture_output=True, text=True)
        findings = (check.stdout + check.stderr).strip()
        if findings:
            print("seqspec check (warnings only):")
            print("\n".join(f"  {line}" for line in findings.splitlines()))

    try:
        spec = load_spec(args.spec)
    except Exception as err:  # pydantic/yaml errors carry the useful message
        fail(f"could not read seqspec file '{args.spec}': {err}")

    modalities = list(spec.modalities)
    if args.modality not in modalities:
        fail(f"the seqspec file has no '{args.modality}' modality (it has: {', '.join(modalities) or 'none'})")

    reads = spec.get_seqspec(args.modality)
    if not reads:
        fail(f"the seqspec file defines no reads for modality '{args.modality}'")

    layouts = {}
    for read in reads:
        try:
            layouts[read.read_id] = read_layout(spec, args.modality, read)
        except Exception as err:
            fail(f"could not place the regions of read '{read.read_id}': {err}")

    def regions_of(read_id, wanted):
        return [(r, s, e) for r, s, e in layouts[read_id] if region_type(r) == wanted]

    # ---------------------------------------------------------------- barcode read
    barcode_reads = [r for r in reads if regions_of(r.read_id, "barcode")]
    umi_reads = [r for r in reads if regions_of(r.read_id, "umi")]
    if not barcode_reads:
        fail("no read of the seqspec file covers a 'barcode' region")
    if len(barcode_reads) > 1:
        fail("cell barcodes are spread over reads "
             f"{', '.join(r.read_id for r in barcode_reads)}; STARsolo and alevin-fry need "
             "every barcode segment and the UMI in one read")
    bc_read = barcode_reads[0]
    if not umi_reads:
        fail("no read of the seqspec file covers a 'umi' region; the pipeline counts UMIs")
    if bc_read not in umi_reads or len(umi_reads) > 1:
        fail(f"the barcode(s) are in read '{bc_read.read_id}' but the UMI in "
             f"{', '.join(r.read_id for r in umi_reads)}; both must be in the same read")

    barcodes = regions_of(bc_read.read_id, "barcode")
    umis = regions_of(bc_read.read_id, "umi")
    if len(umis) > 1:
        fail(f"read '{bc_read.read_id}' has {len(umis)} UMI regions; one is supported")

    # Fixed positions need every region up to the last barcode/UMI base to have one length
    last_stop = max(e for _, _, e in barcodes + umis)
    for region, start, stop in layouts[bc_read.read_id]:
        if start >= last_stop:
            break
        if region.min_len != region.max_len:
            fail(f"region '{region.region_id}' ({region_type(region)}, {region.min_len}-{region.max_len} nt) "
                 f"in read '{bc_read.read_id}' has a variable length before the end of the barcode/UMI, "
                 "so the barcode positions are not fixed")
    for region, start, stop in barcodes + umis:
        if stop - start < region.max_len:
            fail(f"read '{bc_read.read_id}' ({bc_read.max_len} nt) ends inside region "
                 f"'{region.region_id}' ({region_type(region)}, {region.max_len} nt)")

    # ---------------------------------------------------------------- cDNA read
    cdna_reads = [r for r in reads if r is not bc_read and regions_of(r.read_id, "cdna")]
    if not cdna_reads:
        fail("no read other than the barcode read covers a 'cdna' region")
    if len(cdna_reads) > 1:
        print(f"WARNING: several reads cover cDNA ({', '.join(r.read_id for r in cdna_reads)}); "
              f"the first, '{cdna_reads[0].read_id}', is taken as the cDNA read")
    cdna_read = cdna_reads[0]
    cdna_start = min(s for _, s, _ in regions_of(cdna_read.read_id, "cdna"))

    # ---------------------------------------------------------------- STARsolo
    umi_region, umi_start, umi_stop = umis[0]
    if len(barcodes) == 1:
        _, cb_start, cb_stop = barcodes[0]
        star_type = (f"--soloType CB_UMI_Simple --soloCBstart {cb_start + 1} --soloCBlen {cb_stop - cb_start} "
                     f"--soloUMIstart {umi_start + 1} --soloUMIlen {umi_stop - umi_start} "
                     # 0: no read-length check, as fastp may shorten barcode reads
                     "--soloBarcodeReadLength 0")
        star_cbmatch = ""
    else:
        # startAnchor_startPosition_endAnchor_endPosition, anchor 0 = read start, 0-based inclusive
        cb_pos = " ".join(f"0_{s}_0_{e - 1}" for _, s, e in barcodes)
        star_type = (f"--soloType CB_UMI_Complex --soloCBposition {cb_pos} "
                     f"--soloUMIposition 0_{umi_start}_0_{umi_stop - 1}")
        star_cbmatch = "1MM"  # the 1MM_multi* types exist only for CB_UMI_Simple

    # ---------------------------------------------------------------- alevin-fry
    alevin_bc = "1[" + ",".join(f"{s + 1}-{e}" for _, s, e in barcodes) + "]"
    alevin_umi = f"1[{umi_start + 1}-{umi_stop}]"
    alevin_read = f"2[{cdna_start + 1}-end]"

    # ---------------------------------------------------------------- onlists
    if bc_read.strand == "neg":
        print(f"WARNING: barcode read '{bc_read.read_id}' is on the negative strand; "
              "its onlists are used as written")
    base = args.spec.parent
    onlist_rows = []
    for i, (region, _, _) in enumerate(barcodes, start=1):
        onlist = region.onlist
        urltype = (onlist.urltype or "local") if onlist else ""
        url = (onlist.url or onlist.filename) if onlist else ""
        onlist_rows.append((i, region.region_id, urltype, url))
        if onlist and urltype != "local" and args.fetch_onlists:
            out = args.onlist_dir / f"onlist_{i}.txt"
            try:
                barcodes_list = read_remote_list(onlist, str(base))
            except Exception as err:
                fail(f"could not download the onlist of region '{region.region_id}' from {url}: {err}")
            if not barcodes_list:
                fail(f"the onlist of region '{region.region_id}' ({url}) is empty")
            args.onlist_dir.mkdir(parents=True, exist_ok=True)
            out.write_text("\n".join(barcodes_list) + "\n")
            print(f"Onlist {i} ({region.region_id}): {len(barcodes_list)} barcodes from {url}")

    with open(args.out_onlists, "w") as fh:
        for row in onlist_rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    settings = {
        "SEQSPEC_STAR_TYPE": star_type,
        "SEQSPEC_STAR_CBMATCH": star_cbmatch,
        "SEQSPEC_ALEVIN_BC": alevin_bc,
        "SEQSPEC_ALEVIN_UMI": alevin_umi,
        "SEQSPEC_ALEVIN_READ": alevin_read,
        "SEQSPEC_BARCODE_READ": bc_read.read_id,
        "SEQSPEC_CDNA_READ": cdna_read.read_id,
        "SEQSPEC_N_BARCODES": str(len(barcodes)),
    }
    with open(args.out_env, "w") as fh:
        for key, value in settings.items():
            fh.write(f"{key}={shlex.quote(value)}\n")

    # Summary for the task log
    print(f"Barcode/UMI read: {bc_read.read_id} ({bc_read.max_len} nt, strand {bc_read.strand})")
    for region, s, e in layouts[bc_read.read_id]:
        print(f"  {s + 1:>4}-{e:<4} {region_type(region):<10} {region.region_id}")
    print(f"cDNA read: {cdna_read.read_id} (cDNA from base {cdna_start + 1})")
    for key, value in settings.items():
        print(f"{key}={value}")
    for i, rid, urltype, url in onlist_rows:
        print(f"Onlist {i} ({rid}): {url or 'none'}{f' [{urltype}]' if urltype else ''}")


if __name__ == "__main__":
    main()
