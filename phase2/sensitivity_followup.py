#!/usr/bin/env python3
"""
sensitivity_followup.py
-----------------------
Two follow-up questions on the D3 output (no model needed; reads results/sensitivity/):

1. When the window changes (101 or 2001 bp instead of 501), which tracks flip sign?
   Are they tracks with tiny scores, which don't matter, or tracks with big effects?
2. For each SNP, are the brain tracks' signs stable across windows (101 / 501 / 2001)?

Also prints the 1 bp vs 128 bp units line from the log.

Usage (from the repo root)
    python phase2/sensitivity_followup.py
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

SCORES = Path("results/sensitivity/sensitivity_scores.tsv")
LOG = Path("results/sensitivity/sensitivity.log.txt")
CLASSES = Path("phase2/metadata/track_classes.tsv")
for p in (SCORES, LOG, CLASSES):
    if not p.exists():
        sys.exit(f"Not found: {p} (run from the repo root, after scoring_sensitivity.py)")

lab = {int(r["track_index"]): r for r in csv.DictReader(open(CLASSES, encoding="utf-8"), delimiter="\t")
       if r["output_type"] == "DNASE"}
sc = defaultdict(dict)                       # (rsid, window) -> {track: score}
for r in csv.DictReader(open(SCORES, encoding="utf-8"), delimiter="\t"):
    if r["family"] == "window" and r["setting"] in ("101", "501", "2001"):
        sc[(r["rsid"], r["setting"])][int(r["track_index"])] = float(r["score"])

print("Units (from the log):")
for line in open(LOG, encoding="utf-8"):
    if line.startswith("units"):
        print("  " + line.rstrip())

sign = lambda v: (v > 0) - (v < 0)
for rs in dict.fromkeys(k[0] for k in sc):
    base = sc[(rs, "501")]
    ranked = sorted(base, key=lambda t: -abs(base[t]))
    print(f"\n=== {rs} ===")
    print("Sign agreement with 501 bp, by size of the 501 bp score:")
    for w in ("101", "2001"):
        other = sc[(rs, w)]
        cells = []
        for name, tracks in (("top 20", ranked[:20]), ("top 21-100", ranked[20:100]),
                             ("rest, |s|>0.001", [t for t in ranked[100:] if abs(base[t]) > 1e-3])):
            agree = sum(sign(base[t]) == sign(other[t]) for t in tracks)
            cells.append(f"{name}: {agree}/{len(tracks)}")
        print(f"  {w:>4} bp   " + "   ".join(cells))
    print("Brain tracks (score at 101 / 501 / 2001 bp):")
    for t in sorted((t for t in base if lab[t]["analysis_set"] == "brain"), key=lambda t: -abs(base[t])):
        vals = [sc[(rs, w)][t] for w in ("101", "501", "2001")]
        flag = "" if len({sign(v) for v in vals if abs(v) > 1e-3}) <= 1 else "   <- sign changes"
        name = lab[t]["biosample_name"].replace("\\'", "'")    # metadata escapes apostrophes
        print(f"  {vals[0]:+.3f} {vals[1]:+.3f} {vals[2]:+.3f}  {name}{flag}")