#!/usr/bin/env python3
"""
compare_score_runs.py
---------------------
Compare two ISM score files for the same SNPs, track by track: e.g. the same SNPs
scored with two different model weights (fold 0 vs all folds), input lengths, or
software versions. Answers "does this change alter the scores, and does it alter
anything we conclude from them?"

Per SNP x modality (real tracks only; padding dropped by phase2/metadata/track_classes.tsv):
    n_tracks          tracks compared
    spearman          rank agreement of variant_effect_score across tracks
    pearson           linear agreement
    slope             least-squares slope of B on A (1 = same scale)
    sign_agree        share of tracks with the same sign, among tracks where
                      |score| >= --min-effect in A or B (near-zero scores are
                      skipped because their sign is noise)
    median_abs_diff, max_abs_diff
    top_overlap       overlap of the 20 most-affected tracks (|score|) in A and B

Needs only numpy (no scipy).

Usage (repo root)
    python phase2/compare_score_runs.py \\
        --a fold0=results/ism/pilot_ism_scores.tsv \\
        --b allfolds=results/ism/pilot_allfolds_ism_scores.tsv \\
        --out results/summaries/weights_comparison_pilot.tsv
"""
import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np


def _ranks(x):
    """Average ranks (ties share the mean rank), as in Spearman's correlation."""
    order = np.argsort(x, kind="mergesort")
    r = np.empty(len(x))
    r[order] = np.arange(1, len(x) + 1)
    vals, inv, counts = np.unique(x, return_inverse=True, return_counts=True)
    sums = np.bincount(inv, weights=r)
    return (sums / counts)[inv]


def pearson(a, b):
    if np.ptp(a) == 0 or np.ptp(b) == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def spearman(a, b):
    return pearson(_ranks(a), _ranks(b))


def load(spec):
    name, _, path = spec.partition("=")
    if not path:
        sys.exit(f"Use NAME=PATH, got {spec!r}")
    d = defaultdict(dict)
    with open(path, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            d[(r["rsid"], r["modality"].upper())][int(r["track_index"])] = float(r["variant_effect_score"])
    return name, path, d


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--a", required=True, metavar="NAME=PATH")
    ap.add_argument("--b", required=True, metavar="NAME=PATH")
    ap.add_argument("--labels", default="phase2/metadata/track_classes.tsv")
    ap.add_argument("--min-effect", type=float, default=0.05)
    ap.add_argument("--top", type=int, default=20)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    padding = {(r["output_type"], int(r["track_index"])) for r in
               csv.DictReader(open(args.labels, newline=""), delimiter="\t") if r["is_padding"] == "True"}
    na, pa, A = load(args.a)
    nb, pb, B = load(args.b)
    common = sorted(set(A) & set(B))
    only = sorted(set(A) ^ set(B))
    if not common:
        sys.exit("No SNP x modality in common between the two files.")
    if only:
        print(f"NOTE: {len(only)} SNP x modality present in only one file (skipped): {only[:6]}")

    rows = []
    for key in common:
        rsid, mod = key
        tracks = sorted(t for t in set(A[key]) & set(B[key]) if (mod, t) not in padding)
        a = np.array([A[key][t] for t in tracks])
        b = np.array([B[key][t] for t in tracks])
        big = (np.abs(a) >= args.min_effect) | (np.abs(b) >= args.min_effect)
        k = min(args.top, len(tracks))
        top_a = set(np.argsort(-np.abs(a))[:k])
        top_b = set(np.argsort(-np.abs(b))[:k])
        slope = float(np.polyfit(a, b, 1)[0]) if np.ptp(a) > 0 else float("nan")
        rows.append({
            "rsid": rsid, "modality": mod, "n_tracks": len(tracks),
            "spearman": spearman(a, b), "pearson": pearson(a, b), "slope": slope,
            "n_tracks_effect": int(big.sum()),
            "sign_agree": float((np.sign(a[big]) == np.sign(b[big])).mean()) if big.any() else float("nan"),
            "median_abs_diff": float(np.median(np.abs(a - b))), "max_abs_diff": float(np.max(np.abs(a - b))),
            "max_abs_a": float(np.max(np.abs(a))), "max_abs_b": float(np.max(np.abs(b))),
            f"top{k}_overlap": len(top_a & top_b) / k,
        })

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()})

    print(f"A = {na} ({pa})\nB = {nb} ({pb})")
    print(f"{len(common)} SNP x modality compared; min effect for sign agreement {args.min_effect}")
    for mod in sorted({r['modality'] for r in rows}):
        rr = [r for r in rows if r["modality"] == mod]
        def q(k):
            v = [r[k] for r in rr if not np.isnan(r[k])]
            return np.percentile(v, [0, 50, 100]) if v else (float("nan"),) * 3
        topk = [k for k in rr[0] if k.startswith("top")][0]
        print(f"\n{mod} ({len(rr)} SNPs): median [min - max] across SNPs")
        for k in ("spearman", "pearson", "slope", "sign_agree", topk, "median_abs_diff", "max_abs_diff"):
            lo, med, hi = q(k)
            if np.isnan(med):
                print(f"  {k:<16}      n/a  (no tracks with |score| >= {args.min_effect})")
            else:
                print(f"  {k:<16} {med:8.4f}  [{lo:.4f} - {hi:.4f}]")
        worst = sorted(rr, key=lambda r: r["spearman"])[:3]
        print("  lowest Spearman: " + ", ".join(f"{r['rsid']} {r['spearman']:.3f}" for r in worst))
    print(f"\nWritten: {out}")


if __name__ == "__main__":
    main()
