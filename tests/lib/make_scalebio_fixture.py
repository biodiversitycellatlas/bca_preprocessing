#!/usr/bin/env python3
"""
Synthetic ScaleBio-layout reads for tests/checks/scalebio_demux.sh.

Every barcode is random (no ScaleBio barcode lists are used). Each segment's list keeps
its sequences at least 4 substitutions apart, so one substitution is always corrected
to its own barcode and two never reach another one. The reads cycle through the cases
below, and truth.tsv records, for every read, the category the demultiplexer must give
it and, for the reads it keeps, the sample, the synthetic barcode read and the trimmed
cDNA length.

  clean        all barcodes exact
  mm1          one substitution in one barcode            -> Pass (corrected)
  mm2          two substitutions in one barcode           -> BarcodeError
  nbase        an N in one barcode                        -> Pass (corrected)
  umi_n        an N in the UMI                            -> UmiFiltered
  umi_homo     a homopolymer UMI                          -> UmiFiltered
  polya        poly-A tail inside the cDNA                -> Pass, trimmed
  polyg/polyt  poly-G (QuantumScale) / poly-T (v1.1) tail -> Pass, trimmed
  short        poly-A close to the start of the cDNA      -> TooShort
  rt_other     an RT well that no sample uses             -> Unassigned
  pcr_other    a PCR well outside the sample's PCR wells  -> Unassigned
  pool         a non-first sequence of a pooled PCR well  -> Pass (canonical in CB)
  nolinker     (v1.1) splint replaced                     -> LinkerError
  linker_mm    (v1.1) one substitution in the splint      -> Pass

Usage: make_scalebio_fixture.py --layout quantum|v1.1 --outdir DIR [--reads N]
       [--seed S] [--i5-revcomp]
"""

import argparse
import gzip
import os
import random

B = "ACGT"
SPLINT = "TCAGAGC"
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def ham(a, b):
    return sum(x != y for x, y in zip(a, b)) + abs(len(a) - len(b))


def code(rng, n, length, taken=(), min_d=4):
    """n random sequences pairwise >= min_d apart (and from 'taken')."""
    out = []
    pool = list(taken)
    while len(out) < n:
        s = "".join(rng.choice(B) for _ in range(length))
        if all(ham(s, t) >= min_d for t in pool if len(t) == len(s)):
            out.append(s)
            pool.append(s)
    return out


def no_runs(rng, n):
    """Random bases with no run of 3 or more."""
    s = []
    while len(s) < n:
        b = rng.choice(B)
        if len(s) >= 2 and s[-1] == s[-2] == b:
            continue
        s.append(b)
    return "".join(s)


def mutate(rng, s, k, alphabet=B):
    pos = rng.sample(range(len(s)), k)
    s = list(s)
    for p in pos:
        s[p] = rng.choice([b for b in alphabet if b != s[p]])
    return "".join(s)


def wells():
    return [f"{c}{r}" for c in range(1, 13) for r in "ABCDEFGH"]


def write_list(path, entries):
    with open(path, "w") as fh:
        for bc_id, seqs in entries:
            fh.write(bc_id + "\t" + "\t".join(seqs) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--layout", required=True, choices=["quantum", "v1.1"])
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--reads", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--i5-revcomp", action="store_true")
    a = ap.parse_args()
    rng = random.Random(a.seed)
    os.makedirs(a.outdir, exist_ok=True)
    quantum = a.layout == "quantum"

    # ---- barcode lists ------------------------------------------------------
    lists = {}
    rt_len, pcr_len = (8, 8) if quantum else (10, 10)
    lists["rt"] = list(zip(wells(), [[s] for s in code(rng, 96, rt_len)]))
    pcr_seqs = code(rng, 12 * (4 if quantum else 1), pcr_len)
    if quantum:   # pooled wells: 4 i5 sequences each
        lists["pcr"] = [(f"P-{i + 1}", pcr_seqs[4 * i:4 * i + 4]) for i in range(12)]
        for n in (1, 2, 3):
            lists[f"bead{n}"] = [(f"{i + 1:02d}", [s]) for i, s in enumerate(code(rng, 96, 8))]
        segments = ["bead1", "bead2", "bead3", "rt", "pcr"]
        width = {"bead1": 8, "bead2": 8, "bead3": 8, "rt": 8, "pcr": 8}
    else:
        lists["pcr"] = [(f"P-{i + 1}", [s]) for i, s in enumerate(pcr_seqs)]
        lig8 = code(rng, 24, 8)
        lig9 = code(rng, 24, 9)
        lists["lig"] = [(f"L{i + 1}", [s]) for i, s in enumerate(lig8 + lig9)]
        segments = ["lig", "rt", "pcr"]
        width = {"lig": 9, "rt": 10, "pcr": 10}

    # Canonical = first sequence, padded with the first base that keeps them unique
    canon = {}
    for seg in segments:
        for pad in B:
            c = [seqs[0] + pad * (width[seg] - len(seqs[0])) for _, seqs in lists[seg]]
            if len(set(c)) == len(c):
                break
        canon[seg] = c
    wl_paths = []
    for n, seg in enumerate(segments, 1):
        p = os.path.join(a.outdir, f"list_{n}_{seg}.tsv")
        write_list(p, lists[seg])
        wl_paths.append(p)
    with open(os.path.join(a.outdir, "whitelists.txt"), "w") as fh:
        fh.write(" ".join(wl_paths) + "\n")

    # ---- samples ------------------------------------------------------------
    # A: RT 1A-2C (column-major: column 1 then 2A-2C), every PCR well
    # B: RT 3A;3B;5A-5H, PCR P-1..P-6     C: RT 3A;3B, PCR P-7..P-12
    w = wells()
    samples = {
        "sampleA": ("1A-2C", "", w[:11], [x for x, _ in lists["pcr"]]),
        "sampleB": ("3A;3B;5A-5H", "P-1-P-6", ["3A", "3B"] + w[32:40], [f"P-{i}" for i in range(1, 7)]),
        "sampleC": ("3A;3B", "P-7-P-12", ["3A", "3B"], [f"P-{i}" for i in range(7, 13)]),
    }
    with open(os.path.join(a.outdir, "samples.tsv"), "w") as fh:
        for s, (rt_spec, pcr_spec, _, _) in samples.items():
            fh.write(f"{s}\t{rt_spec}\t{pcr_spec}\n")
    assign = {}
    for s, (_, _, rts, pcrs) in samples.items():
        for r in rts:
            for p in pcrs:
                assign[(r, p)] = s
    used_rt = {r for r, _ in assign}
    free_rt = [x for x in w if x not in used_rt]
    rt_index = {x: i for i, (x, _) in enumerate(lists["rt"])}
    pcr_index = {x: i for i, (x, _) in enumerate(lists["pcr"])}

    cases = ["clean", "mm1", "mm2", "nbase", "umi_n", "umi_homo", "polya",
             "polyg" if quantum else "polyt", "short", "rt_other", "pcr_other", "pool"]
    if not quantum:
        cases += ["nolinker", "linker_mm"]
    if not quantum:
        cases.remove("pool")

    # ---- reads --------------------------------------------------------------
    names = ["R1", "R2", "I1", "I2"] if quantum else ["R1", "R2", "I2"]
    fhs = {n: gzip.open(os.path.join(a.outdir, f"lib_S1_L001_{n}_001.fastq.gz"), "wt") for n in names}
    truth = open(os.path.join(a.outdir, "truth.tsv"), "w")
    truth.write("read\tcase\tcategory\tsample\tbc_read\tcdna_len\n")

    def put(n, rid, seq):
        fhs[n].write(f"@{rid}\n{seq}\n+\n{'I' * len(seq)}\n")

    for i in range(a.reads):
        case = cases[i % len(cases)]
        rid = f"r{i}"
        # choose an assigned (rt, pcr) unless the case asks otherwise
        sample = rng.choice(list(samples))
        _, _, rts, pcrs = samples[sample]
        rt, pcr = rng.choice(rts), rng.choice(pcrs)
        if case == "rt_other":
            rt = rng.choice(free_rt)
        if case == "pcr_other":
            sample = "sampleB"
            rt, pcr = "5C", "P-9"          # 5C is sampleB's only, P-9 not in sampleB's wells
        pcr_choice = 0
        if case == "pool":
            pcr_choice = rng.randrange(1, 4)
        idx = {"rt": rt_index[rt], "pcr": pcr_index[pcr]}
        for seg in segments:
            if seg not in idx:
                idx[seg] = rng.randrange(len(lists[seg]))
        obs = {seg: lists[seg][idx[seg]][1][pcr_choice if seg == "pcr" else 0] for seg in segments}
        if case in ("mm1", "mm2", "nbase"):
            seg = rng.choice(segments)
            k = 2 if case == "mm2" else 1
            obs[seg] = mutate(rng, obs[seg], k, "N" if case == "nbase" else B)
        umi = no_runs(rng, 8)
        if case == "umi_n":
            umi = umi[:3] + "N" + umi[4:]
        if case == "umi_homo":
            umi = "G" * 8

        # cDNA: 60 nt without runs; a tail starting at 'cut' where the case asks for one
        cdna = no_runs(rng, 60)
        cut = 60
        tail = {"polya": "A", "polyg": "G", "polyt": "T", "short": "A"}.get(case)
        if tail:
            n = 12 if tail == "G" else 8
            cut = 10 if case == "short" else rng.randrange(20, 45)
            while cdna[cut - 1] == tail:
                cdna = no_runs(rng, 60)
            cdna = cdna[:cut] + tail * n + cdna[cut + n:]
            cdna = cdna[:60]

        i2 = obs["pcr"]
        if a.i5_revcomp:
            i2 = rc(i2)
        if quantum:
            r1 = cdna
            r2 = umi + obs["rt"] + no_runs(rng, 30)
            i1 = obs["bead1"] + "ACGT" + obs["bead2"] + "TGCA" + obs["bead3"] + "GT"
            put("I1", rid, i1)
        else:
            splint = SPLINT
            if case == "nolinker":
                splint = "GGTTCCA"
            if case == "linker_mm":
                splint = SPLINT[:3] + ("T" if SPLINT[3] != "T" else "A") + SPLINT[4:]
            lig = obs["lig"]
            # The splint one base off overlaps itself or the ligation barcode and differs
            # from TCAGAGC at 2+ positions, so it is only ever found where it was put
            r1 = lig + splint + umi + obs["rt"] + no_runs(rng, 20)
            r2 = cdna
        put("R1", rid, r1)
        put("R2", rid, r2)
        put("I2", rid, i2)

        # expectation
        if case == "nolinker":
            cat = "LinkerError"
        elif case == "mm2":
            cat = "BarcodeError"
        elif case in ("umi_n", "umi_homo"):
            cat = "UmiFiltered"
        elif case == "short":
            cat = "TooShort"
        elif (rt, pcr) not in assign:
            cat = "Unassigned"
        else:
            cat = "Pass"
        if cat == "Pass":
            if quantum:
                # the cell's list entries as one number (mixed radix, list sizes), 16 base-4 digits
                n = 0
                for seg in segments:
                    n = n * len(lists[seg]) + idx[seg]
                cb = "".join("ACGT"[(n >> (2 * k)) & 3] for k in range(15, -1, -1))
            else:
                cb = "".join(canon[seg][idx[seg]] for seg in segments)
            truth.write(f"{rid}\t{case}\t{cat}\t{assign[(rt, pcr)]}\t{cb}{umi}\t{cut}\n")
        else:
            truth.write(f"{rid}\t{case}\t{cat}\t\t\t\n")

    for fh in fhs.values():
        fh.close()
    truth.close()


if __name__ == "__main__":
    main()
