#!/usr/bin/env python3
"""Pick one bin's matched pool from binned candidates.

Step 6 of the matched-SNP design (Choo, Nov 18: "set up one bin first, then
automate"). Input is the candidate table after liftover, distance, MAF and
assign_bins.py, binned with the SAME edges file as the leads. This script:
  1. prints how many candidates fall in every bin (pool size before the draw)
  2. chooses the bin: --bin D4_M4, or (default) the bin holding the most leads
     in --leads-binned
  3. keeps candidates whose MAF came from 1000 Genomes (as for the leads)
  4. shuffles them (fixed seed) and takes them in that order, skipping any
     within --min-spacing bp of one already taken (independent SNPs, no shared LD)
     until --n are chosen
  5. compares the chosen pool with the leads in that bin (MAF and distance),
     as the match-quality check before any model run

Output: the pool in the lead-table format (rsid chr pos a1 a2 ... beta) plus
pool_bin, pool_seed, pool_order, so it can go straight to ism_with_cache.py.

Usage (from the repo root)
  python3 phase1/select_bin_pool.py phase1/data/pool_candidates_scz_20k_hg38.bins_tss.tsv \\
      --leads-binned phase1/data/scz_pgc3_lead_snps_hg38.bins_tss.tsv \\
      --n 20 --seed 2026 --out phase1/data/pool_scz_tss_<bin>_n20.tsv
  ("<bin>" in --out is replaced by the chosen bin.)

Requires: Python 3 standard library.
"""

import argparse
import bisect
import csv
import random
import statistics as st
import sys
from datetime import datetime
from pathlib import Path


def read(path):
    with open(path, newline="", encoding="utf-8") as f:
        rd = csv.reader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
        head = next(rd)
        return head, [r for r in rd if r]


def bin_key(b):
    d, m = b.split("_")
    return int(d[1:]), int(m[1:])


def main():
    ap = argparse.ArgumentParser(description="Draw one bin's matched pool.")
    ap.add_argument("candidates")
    ap.add_argument("--leads-binned", required=True)
    ap.add_argument("--bin")
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--min-spacing", type=int, default=1_000_000)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    ch, crows = read(args.candidates)
    lh, lrows = read(args.leads_binned)
    for h, name in ((ch, args.candidates), (lh, args.leads_binned)):
        miss = [c for c in ("bin", "maf_used", "maf_source", "dist_used_bp", "chr", "pos") if c not in h]
        if miss:
            sys.exit(f"{name} lacks {miss}: run assign_bins.py on it first")
    C = {c: i for i, c in enumerate(ch)}
    L = {c: i for i, c in enumerate(lh)}

    lines = []

    def log(msg=""):
        print(msg)
        lines.append(msg)

    log(f"select_bin_pool.py  {datetime.now().isoformat(timespec='seconds')}")
    log(f"Candidates: {args.candidates} ({len(crows):,})")
    log(f"Leads: {args.leads_binned} ({len(lrows)})")

    lead_n, cand_n, cand_1kg = {}, {}, {}
    for r in lrows:
        if r[L["bin"]]:
            lead_n[r[L["bin"]]] = lead_n.get(r[L["bin"]], 0) + 1
    for r in crows:
        b = r[C["bin"]]
        if b:
            cand_n[b] = cand_n.get(b, 0) + 1
            if r[C["maf_source"]] == "1000G_EUR":
                cand_1kg[b] = cand_1kg.get(b, 0) + 1
    bins = sorted(set(lead_n) | set(cand_n), key=bin_key)
    log("")
    log("Per bin: leads / candidates (with 1000G MAF):")
    nd = max(bin_key(b)[0] for b in bins)
    nm = max(bin_key(b)[1] for b in bins)
    log("        " + "".join(f"{'M' + str(m):>16}" for m in range(1, nm + 1)))
    for d in range(1, nd + 1):
        cells = []
        for m in range(1, nm + 1):
            b = f"D{d}_M{m}"
            cells.append(f"{lead_n.get(b, 0):>4} / {cand_1kg.get(b, 0):>6,}")
        log(f"  D{d}   " + "".join(f"{c:>16}" for c in cells))
    tot = sum(cand_1kg.values())
    log(f"  candidates binned with 1000G MAF: {tot:,}; smallest bin {min(cand_1kg.get(b, 0) for b in bins):,}")

    target = args.bin or max(sorted(lead_n, key=bin_key), key=lambda b: lead_n[b])
    if target not in bins:
        sys.exit(f"Bin {target} not found")
    log("")
    log(f"Target bin: {target} ({'given' if args.bin else 'most leads'}; {lead_n.get(target, 0)} leads)")

    pool = [r for r in crows if r[C["bin"]] == target and r[C["maf_source"]] == "1000G_EUR"]
    other = sum(1 for r in crows if r[C["bin"]] == target) - len(pool)
    log(f"  candidates in bin: {len(pool) + other:,}; with 1000G MAF: {len(pool):,}")
    rnd = random.Random(args.seed)
    order = list(range(len(pool)))
    rnd.shuffle(order)
    chosen, taken, skipped = [], {}, 0
    for k in order:
        r = pool[k]
        c, p = r[C["chr"]], int(r[C["pos"]])
        v = taken.setdefault(c, [])
        j = bisect.bisect_left(v, p)
        if (j < len(v) and v[j] - p < args.min_spacing) or (j > 0 and p - v[j - 1] < args.min_spacing):
            skipped += 1
            continue
        v.insert(j, p)
        chosen.append(r)
        if len(chosen) == args.n:
            break
    log(f"  drawn {len(chosen)} (seed {args.seed}); skipped {skipped} within {args.min_spacing:,} bp of one already drawn")
    if len(chosen) < args.n:
        log(f"  WARNING: only {len(chosen)} of {args.n} could be drawn")

    def summ(vals):
        return f"median {st.median(vals):.3f} [{min(vals):.3f}-{max(vals):.3f}]" if vals else "n/a"

    lm = [float(r[L["maf_used"]]) for r in lrows if r[L["bin"]] == target]
    ld = [int(r[L["dist_used_bp"]]) / 1e3 for r in lrows if r[L["bin"]] == target]
    pm = [float(r[C["maf_used"]]) for r in chosen]
    pd = [int(r[C["dist_used_bp"]]) / 1e3 for r in chosen]
    log("  match check (leads in bin vs pool):")
    log(f"    MAF       leads {summ(lm)}   pool {summ(pm)}")
    log(f"    dist (kb) leads {summ(ld)}   pool {summ(pd)}")

    out = Path(args.out.replace("<bin>", target))
    with open(out, "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE, escapechar="\\")
        wr.writerow(ch + ["pool_bin", "pool_seed", "pool_order"])
        for i, r in enumerate(chosen, 1):
            wr.writerow(r + [target, str(args.seed), str(i)])
    log(f"Written: {out}")
    lp = out.with_name((out.name[:-4] if out.name.endswith(".tsv") else out.name) + ".log.txt")
    lp.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Log: {lp}")


if __name__ == "__main__":
    main()
