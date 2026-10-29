#!/usr/bin/env python3
"""
score_tracks.py
---------------
Join ISM variant-effect scores to the track labels, drop padding, and write one
tidy table (one row per SNP x track) plus per-SNP means and variances for each
tissue set. Everything downstream (plots, specificity statistics, controls)
reads the tidy table, so this is the single place where scores meet labels.

Inputs
    --scores NAME=PATH [NAME=PATH ...]
        ISM score files written by run_ism_pilot.py (current format, with
        a1/a2/hg38_ref_is/risk_oriented_score columns). NAME labels the set,
        e.g. pilot=results/ism/pilot_ism_scores.tsv control=results/ism/control_ism_scores.tsv
    --classes  phase2/metadata/track_classes.tsv  (from build_track_classes.py)

Outputs (default folder results/summaries/)
    track_scores_long.tsv   one row per set x SNP x real track
    track_scores_by_set.tsv per set x SNP x modality x analysis_set:
                            n tracks, mean / variance / median / max of |score|,
                            mean risk-oriented score
    track_scores.log.txt    inputs (with md5), row counts, padding dropped

Checks (the script stops if any fails)
    - every score row matches a track in the label file (channel numbering intact)
    - padding channels sit only at the end of each head
    - after dropping padding, every SNP has exactly the head's number of real
      tracks (e.g. 305 DNase), so a truncated or duplicated run cannot pass
    - each SNP has one set of alleles across all its rows

|score| is used for magnitude summaries. The direction of effect is carried in
risk_oriented_score (positive = the schizophrenia risk allele opens chromatin);
see docs/allele_conventions.md.

Usage (from the repo root)
    python phase2/score_tracks.py \\
        --scores pilot=results/ism/pilot_ism_scores.tsv control=results/ism/control_ism_scores.tsv
"""

import argparse
import csv
import datetime
import hashlib
import math
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

REQUIRED = ["rsid", "chr", "pos", "a1", "a2", "beta", "hg38_ref", "hg38_alt", "hg38_ref_is",
            "palindromic", "modality", "track_index", "variant_effect_score", "risk_oriented_score"]
SNP_FIELDS = ["chr", "pos", "a1", "a2", "beta", "hg38_ref", "hg38_alt", "hg38_ref_is", "palindromic"]
LABEL_FIELDS = ["biosample_name", "biosample_type", "ontology_curie", "stage", "life_stage_group",
                "group", "analysis_set", "origin", "cell_class", "brain_region", "histone_mark"]
SETS = ["brain", "neural_model", "excluded", "background"]


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_classes(path):
    rows = list(csv.DictReader(open(path, newline="", encoding="utf-8"), delimiter="\t"))
    by_key = {(r["output_type"], int(r["track_index"])): r for r in rows}
    real = Counter(r["output_type"] for r in rows if r["is_padding"] != "True")
    # padding must be a contiguous block at the end of each head
    for head in {r["output_type"] for r in rows}:
        flags = [r["is_padding"] == "True" for r in sorted((r for r in rows if r["output_type"] == head),
                                                        key=lambda r: int(r["track_index"]))]
        if True in flags and not all(flags[flags.index(True):]):
            sys.exit(f"{path}: padding in {head} is not a single block at the end")
    return by_key, real


def summarise(values):
    a = [abs(v) for v in values]
    return dict(n=len(a), mean_abs=statistics.fmean(a) if a else math.nan,
                var_abs=statistics.pvariance(a) if len(a) > 1 else math.nan,
                median_abs=statistics.median(a) if a else math.nan,
                max_abs=max(a) if a else math.nan)


def fmt(x):
    return "" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.6f}"


def main():
    ap = argparse.ArgumentParser(description="Join ISM scores to track labels; tidy table + per-set summaries.")
    ap.add_argument("--scores", nargs="+", required=True, metavar="NAME=PATH")
    ap.add_argument("--classes", default="phase2/metadata/track_classes.tsv")
    ap.add_argument("--outdir", default="results/summaries")
    args = ap.parse_args()

    inputs = []
    for item in args.scores:
        if "=" not in item:
            sys.exit(f"--scores expects NAME=PATH, got {item!r}")
        name, path = item.split("=", 1)
        if not Path(path).exists():
            sys.exit(f"Score file not found: {path}")
        inputs.append((name, path))
    if len({n for n, _ in inputs}) != len(inputs):
        sys.exit("Each --scores NAME must be different")
    if not Path(args.classes).exists():
        sys.exit(f"Label file not found: {args.classes} (run phase2/build_track_classes.py)")

    classes, n_real = load_classes(args.classes)
    long_rows, problems, log_lines = [], [], []

    for set_name, path in inputs:
        reader = csv.DictReader(open(path, newline="", encoding="utf-8"), delimiter="\t")
        missing = [c for c in REQUIRED if c not in (reader.fieldnames or [])]
        if missing:
            hint = (" This looks like the old format (true_ref/true_alt); re-run run_ism_pilot.py."
                    if "true_ref" in (reader.fieldnames or []) else "")
            sys.exit(f"{path} is missing columns {missing}.{hint}")
        per_snp_head = Counter()
        seen = set()
        snp_info = {}
        n_in = n_pad = 0
        for r in reader:
            n_in += 1
            head = r["modality"].upper()
            key = (head, int(r["track_index"]))
            lab = classes.get(key)
            if lab is None:
                problems.append(f"{path}: {r['rsid']} {head} track {key[1]} has no label")
                continue
            if (r["rsid"], key) in seen:
                problems.append(f"{path}: duplicate row {r['rsid']} {head} track {key[1]}")
            seen.add((r["rsid"], key))
            info = tuple(r[f] for f in SNP_FIELDS)
            if snp_info.setdefault(r["rsid"], info) != info:
                problems.append(f"{path}: {r['rsid']} has inconsistent alleles/position across rows")
            if lab["is_padding"] == "True":
                n_pad += 1
                continue
            per_snp_head[(r["rsid"], head)] += 1
            s = float(r["variant_effect_score"])
            ro = r["risk_oriented_score"]
            long_rows.append({"set": set_name, "rsid": r["rsid"], **{f: r[f] for f in SNP_FIELDS},
                              "modality": head, "track_index": key[1],
                              **{f: lab[f] for f in LABEL_FIELDS},
                              "variant_effect_score": r["variant_effect_score"], "abs_score": fmt(abs(s)),
                              "risk_oriented_score": ro})
        for (rsid, head), n in sorted(per_snp_head.items()):
            if n != n_real[head]:
                problems.append(f"{path}: {rsid} has {n} real {head} tracks, expected {n_real[head]}")
        log_lines.append(f"{set_name:<10} {path}  md5 {md5(path)}  rows {n_in}  padding dropped {n_pad}  "
                         f"real kept {n_in - n_pad}  SNPs {len(snp_info)}")
        print(f"{set_name:<10} {n_in:>6} rows  {n_pad:>5} padding dropped  {n_in - n_pad:>6} kept  {len(snp_info)} SNPs")

    if problems:
        sys.exit("STOPPED:\n  " + "\n  ".join(problems[:20]) +
                 (f"\n  ... and {len(problems) - 20} more" if len(problems) > 20 else ""))

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    long_path = out / "track_scores_long.tsv"
    with open(long_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(long_rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(long_rows)

    # ---- per SNP x modality x analysis_set summaries
    groups = defaultdict(list)
    meta = {}
    for r in long_rows:
        groups[(r["set"], r["rsid"], r["modality"], r["analysis_set"])].append(r)
        meta[(r["set"], r["rsid"])] = r
    summary = []
    for (sname, rsid, head, aset), rows in groups.items():
        st = summarise([float(r["variant_effect_score"]) for r in rows])
        risk = [float(r["risk_oriented_score"]) for r in rows if r["risk_oriented_score"] != ""]
        m = meta[(sname, rsid)]
        summary.append({"set": sname, "rsid": rsid, "modality": head, "analysis_set": aset,
                        "hg38_ref_is": m["hg38_ref_is"], "palindromic": m["palindromic"],
                        "n": st["n"], "mean_abs": fmt(st["mean_abs"]), "var_abs": fmt(st["var_abs"]),
                        "median_abs": fmt(st["median_abs"]), "max_abs": fmt(st["max_abs"]),
                        "mean_risk_oriented": fmt(statistics.fmean(risk)) if risk else ""})
    order = {s: i for i, s in enumerate(SETS)}
    summary.sort(key=lambda r: (r["set"], r["modality"], r["rsid"], order.get(r["analysis_set"], 9)))
    sum_path = out / "track_scores_by_set.tsv"
    with open(sum_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(summary)

    with open(out / "track_scores.log.txt", "w", encoding="utf-8") as f:
        f.write(f"date    : {datetime.datetime.now().isoformat(timespec='seconds')}\n")
        f.write(f"labels  : {args.classes}  md5 {md5(args.classes)}\n")
        for line in log_lines:
            f.write(f"scores  : {line}\n")
        f.write(f"long table rows: {len(long_rows)}\n")

    # ---- console view: DNase, mean |score| by set
    idx = {(r["set"], r["rsid"], r["modality"], r["analysis_set"]): r for r in summary}
    for head in ("DNASE",):
        print(f"\n{head}: mean |score| (variance) by track set")
        print(f"{'set':<9}{'SNP':<13}" + "".join(f"{s:>22}" for s in SETS))
        for sname, rsid in sorted({(r['set'], r['rsid']) for r in summary if r['modality'] == head}):
            cells = []
            for s in SETS:
                r = idx.get((sname, rsid, head, s))
                cells.append(f"{float(r['mean_abs']):.4f} ({float(r['var_abs']):.4f}) n={r['n']:<3}" if r else "-")
            print(f"{sname:<9}{rsid:<13}" + "".join(f"{c:>22}" for c in cells))
    print(f"\nLONG   : {long_path}  ({len(long_rows)} rows)\nBY SET : {sum_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())