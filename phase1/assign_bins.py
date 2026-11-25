#!/usr/bin/env python3
"""Assign SNPs to matching bins (European MAF x distance to gene).

Step 5 of the matched-SNP design (Choo, Nov 18): each risk SNP is later
normalised by non-associated SNPs from its own bin, so this script decides the
bins. It reads a features table made by add_matching_features.py (distance)
then add_eur_maf.py (MAF), and writes the table back with bin columns.

Two distance definitions, chosen with --distance (open question for Choo,
slides_session_F slide 4):
  tss   |pc_tss_dist_bp|: distance to the nearest protein-coding TSS, cut into
        --n-dist quantile bins. No tie at 0, so the bins are balanced.
  body  pc_body_dist_bp: Choo's sketch. About two thirds of leads lie inside a
        gene (distance exactly 0), so quantiles collapse; instead bin 1 is
        "inside a gene" and the intergenic SNPs are cut into --n-dist - 1
        quantile bins.
MAF is cut into --n-maf quantile bins.

Edges: by default they are quantiles of this table (--edges-from-table). They
are written to a JSON file (--save-edges) so the same edges can be applied to
another table with --edges (e.g. SCZ edges applied to T1D leads, or later
edges from genome-wide common SNPs applied to every trait). A value equal to an
edge goes to the upper bin. Bins are numbered from 1 (nearest to a gene,
rarest).

MAF source: eur_maf (1000 Genomes EUR). Where it is blank (SNP not in the
panel), the GWAS's own a1 frequency is folded to MAF instead
(--fallback-freq-col, default a1_freq_controls or a1_freq), and maf_source
says so. Agreement between the two was r = 0.997 for SCZ and T1D leads.

Columns added:
  maf_used  maf_source  dist_used_bp  maf_bin  dist_bin  bin   (bin = "D<d>_M<m>")

Usage (from the repo root)
  # SCZ: edges from the SCZ leads, TSS distance; save the edges
  python3 phase1/assign_bins.py phase1/data/scz_pgc3_lead_snps_hg38.features.tsv \\
      --distance tss --save-edges phase1/data/bins_scz_leads_tss.json
  # T1D with the SCZ edges, to see how another trait fills the same bins
  python3 phase1/assign_bins.py phase1/data/t1d_chiou2021_lead_snps_hg38.features.tsv \\
      --distance tss --edges phase1/data/bins_scz_leads_tss.json \\
      --out phase1/data/t1d_chiou2021_lead_snps_hg38.bins_tss_sczedges.tsv
  Default output: <input stem without .features>.bins_<distance>.tsv (+ .log.txt)

Requires: numpy.
"""

import argparse
import csv
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

ADDED = ["maf_used", "maf_source", "dist_used_bp", "maf_bin", "dist_bin", "bin"]


def quantile_edges(x, n):
    """n-1 inner edges at the k/n quantiles (linear interpolation)."""
    return [float(v) for v in np.quantile(np.asarray(x, dtype=float), [k / n for k in range(1, n)])]


def to_bin(x, edges):
    return int(np.searchsorted(np.asarray(edges, dtype=float), x, side="right")) + 1


def main():
    ap = argparse.ArgumentParser(description="Assign MAF x distance matching bins.")
    ap.add_argument("features_file")
    ap.add_argument("--distance", choices=["tss", "body"], required=True)
    ap.add_argument("--n-maf", type=int, default=5)
    ap.add_argument("--n-dist", type=int, default=5)
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--edges", help="JSON edges file from an earlier --save-edges run")
    g.add_argument("--edges-from-table", action="store_true", help="quantiles of this table (default)")
    ap.add_argument("--save-edges")
    ap.add_argument("--fallback-freq-col")
    ap.add_argument("--out")
    args = ap.parse_args()

    src = Path(args.features_file)
    if not src.exists():
        sys.exit(f"Not found: {src}")
    if args.out:
        out = Path(args.out)
    else:
        stem = src.name[:-4] if src.name.endswith(".tsv") else src.name
        stem = stem[:-9] if stem.endswith(".features") else stem
        out = src.with_name(f"{stem}.bins_{args.distance}.tsv")
    if out.resolve() == src.resolve():
        sys.exit("--out would overwrite the input")
    log_path = out.with_name(out.name[:-4] + ".log.txt" if out.name.endswith(".tsv") else out.name + ".log.txt")

    lines = []

    def log(msg=""):
        print(msg)
        lines.append(msg)

    with open(src, newline="", encoding="utf-8") as f:
        rd = csv.reader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
        head = next(rd)
        rows = [r for r in rd if r]
    need = ["eur_maf", "pc_tss_dist_bp" if args.distance == "tss" else "pc_body_dist_bp"]
    miss = [c for c in need if c not in head]
    if miss:
        sys.exit(f"{src} lacks {', '.join(miss)}: run add_matching_features.py and add_eur_maf.py first")
    clash = [c for c in ADDED if c in head]
    if clash:
        sys.exit(f"{src} already has {', '.join(clash)}")
    fb = args.fallback_freq_col or next((c for c in ("a1_freq_controls", "a1_freq") if c in head), None)
    H = {c: i for i, c in enumerate(head)}
    rid = (lambda r: r[H["rsid"]]) if "rsid" in H else (lambda r: f"{r[H['chr']]}:{r[H['pos']]}")

    log(f"assign_bins.py  {datetime.now().isoformat(timespec='seconds')}")
    log(f"Input: {src}  ({len(rows)} SNPs); distance = {args.distance}")

    # features per SNP
    maf, msrc, dist, dropped = [], [], [], []
    for r in rows:
        m, s = None, ""
        if r[H["eur_maf"]] != "":
            m, s = float(r[H["eur_maf"]]), "1000G_EUR"
        elif fb and r[H[fb]] not in ("", "NA", "nan"):
            f1 = float(r[H[fb]])
            m, s = min(f1, 1 - f1), f"gwas_fallback:{fb}"
        dcol = "pc_tss_dist_bp" if args.distance == "tss" else "pc_body_dist_bp"
        d = abs(int(r[H[dcol]])) if r[H[dcol]] != "" else None
        maf.append(m)
        msrc.append(s)
        dist.append(d)
        if m is None or d is None:
            dropped.append(f"  {rid(r)}: " + ("no MAF" if m is None else "no distance"))
    for i, s in enumerate(msrc):
        if s.startswith("gwas_fallback"):
            log(f"MAF fallback: {rid(rows[i])} uses {s.split(':', 1)[1]} (MAF {maf[i]:.3f})")
    if dropped:
        log(f"No bin (missing feature) for {len(dropped)} SNP(s):")
        for d in dropped:
            log(d)
    ok = [i for i in range(len(rows)) if maf[i] is not None and dist[i] is not None]

    # edges
    if args.edges:
        E = json.loads(Path(args.edges).read_text())
        if E["distance"] != args.distance:
            sys.exit(f"{args.edges} has distance '{E['distance']}', not '{args.distance}'")
        maf_e, dist_e = E["maf_edges"], E["dist_edges_bp"]
        log(f"Edges from {args.edges} (made from {E.get('source')}, n = {E.get('n_snps')})")
    else:
        maf_e = quantile_edges([maf[i] for i in ok], args.n_maf)
        if args.distance == "tss":
            dist_e = quantile_edges([dist[i] for i in ok], args.n_dist)
        else:
            inter = [dist[i] for i in ok if dist[i] > 0]
            if len(inter) < args.n_dist - 1:
                sys.exit("Too few intergenic SNPs to cut into bins")
            # bin 1 = inside a gene (d = 0). Edge 0.5 puts every d >= 1 above it.
            dist_e = [0.5] + (quantile_edges(inter, args.n_dist - 1) if args.n_dist > 2 else [])
        E = {"distance": args.distance, "maf_edges": maf_e, "dist_edges_bp": dist_e,
             "source": str(src), "n_snps": len(ok), "created": datetime.now().isoformat(timespec="seconds"),
             "rule": "inner edges; a value equal to an edge goes to the upper bin; bins numbered from 1"}
        log(f"Edges: quantiles of this table ({len(ok)} SNPs)")
    n_m, n_d = len(maf_e) + 1, len(dist_e) + 1
    log("  MAF edges: " + " / ".join(f"{e:.3f}" for e in maf_e))
    if args.distance == "body":
        log("  distance edges (bp): inside a gene | " + " / ".join(f"{e:,.0f}" for e in dist_e[1:]))
    else:
        log("  |TSS| edges (bp): " + " / ".join(f"{e:,.0f}" for e in dist_e))
    if args.save_edges:
        Path(args.save_edges).write_text(json.dumps(E, indent=2) + "\n")
        log(f"  saved to {args.save_edges}")

    # assign
    grid = np.zeros((n_d, n_m), dtype=int)
    out_rows = []
    okset = set(ok)
    for i, r in enumerate(rows):
        if i in okset:
            mb, db = to_bin(maf[i], maf_e), to_bin(dist[i], dist_e)
            grid[db - 1, mb - 1] += 1
            out_rows.append(r + [f"{maf[i]:.6f}", msrc[i], str(dist[i]), str(mb), str(db), f"D{db}_M{mb}"])
        else:
            out_rows.append(r + ["", msrc[i], "" if dist[i] is None else str(dist[i]), "", "", ""])
    with open(out, "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE, escapechar="\\")
        wr.writerow(head + ADDED)
        wr.writerows(out_rows)

    # grid
    def dlabel(d):
        if args.distance == "body" and d == 1:
            return "D1 inside gene"
        lo = 0 if d == 1 else dist_e[d - 2]
        hi = dist_e[d - 1] if d <= len(dist_e) else None
        if hi is None:
            return f"D{d} > {lo / 1e3:,.0f} kb"
        if lo in (0, 0.5):
            return f"D{d} < {hi / 1e3:,.0f} kb" + (" (outside)" if lo == 0.5 else "")
        return f"D{d} {lo / 1e3:,.0f}-{hi / 1e3:,.0f} kb"
    log("")
    log("Leads per bin (rows: distance, 1 = nearest; columns: MAF, 1 = rarest):")
    log(f"  {'':<22}" + "".join(f"{'M' + str(m):>6}" for m in range(1, n_m + 1)) + f"{'total':>8}")
    for d in range(1, n_d + 1):
        log(f"  {dlabel(d):<22}" + "".join(f"{v:>6}" for v in grid[d - 1]) + f"{grid[d - 1].sum():>8}")
    log(f"  {'total':<22}" + "".join(f"{v:>6}" for v in grid.sum(0)) + f"{grid.sum():>8}")
    nz = grid[grid > 0]
    log(f"Cells: {grid.size}; empty {int((grid == 0).sum())}; with 1 lead {int((grid == 1).sum())}; "
        f"leads per occupied cell median {int(np.median(nz)) if nz.size else 0}, max {int(grid.max())}")
    log(f"Written: {out}")
    log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Log: {log_path}")


if __name__ == "__main__":
    main()
