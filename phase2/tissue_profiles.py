#!/usr/bin/env python3
"""
tissue_profiles.py
------------------
Per-SNP tissue profiles from per-track NORMALISED scores (Choo, Nov 5).

The problem it solves
    Raw variant scores are not comparable between tracks: each track has its own
    dynamic range (AlphaGenome docs: scores from "different tracks and modalities
    ... are on different scales"). So each score is first converted to a quantile
    WITHIN ITS OWN TRACK, against a reference set of variants:

        q(SNP, track) = share of reference variants whose |score| in that track is
                        smaller (ties count half)

    q = 0.95 in putamen means "a bigger effect in putamen than 95% of the reference
    variants have in putamen". Quantiles are on the same 0-1 scale in every track,
    so tissues can be compared.

Reference set
    --reference NAME   which score set (the `set` column of track_scores_long.tsv)
                       defines each track's distribution.
      * Final analysis: a large set of common, non-associated null SNPs (AlphaGenome
        itself uses ~300K common gnomAD variants). Then q answers "how unusual is
        this effect for a common variant in this tissue?"
      * Interim (before the null set exists): the SCZ lead SNPs themselves. Then q
        answers "how unusual is this effect compared with the other SCZ leads in
        this tissue?". This removes track-scale differences and supports
        per-SNP tissue profiles, but CANNOT test whether SCZ variants are more
        brain-active than common variants in general (every track's quantiles
        average 0.5 across the reference set by construction).

Per SNP x tissue category
    mean quantile over the category's main-analysis tracks (tissues, primary cells,
    stem-cell-derived; cell lines excluded unless --include-cell-lines), each
    category weighted equally. Categories from phase2/metadata/track_categories.tsv.

Per SNP summary
    brain_q          mean quantile over CNS tracks
    brain_prenatal_q / brain_postnatal_q   same, split by life stage
    body_mean_q      mean over the non-neural body categories (equal weight)
    top_body_category, top_body_q          the body category with the highest mean quantile
                     (only categories with >= --min-top tracks, default 3)
    specificity      brain_q - body_mean_q
    brain_rank       rank of CNS among all categories (1 = highest)
    profile_class    descriptive, threshold --high (default 0.75):
                       brain_specific  brain_q >= high, top_body_q <  high
                       shared          brain_q >= high, top_body_q >= high
                       other_tissue    brain_q <  high, top_body_q >= high
                       none            neither
    mean_risk_oriented_brain   mean raw risk-oriented score over CNS tracks (direction:
                               + = risk allele opens chromatin)

Figures (results/profiles/)
    <mod>_heatmap.png      SNPs x categories, mean quantile (sorted by specificity)
    <mod>_slope.png        each SNP: body tissues (mean) -> brain
    <mod>_scatter.png      brain q vs top body-category q, four profile classes
    <mod>_closeups.png     track-level quantiles by category for the most brain-high
                           and most brain-low SNPs

Usage (from the repo root)
    python phase2/tissue_profiles.py --sets lead control --reference lead \\
        [--annotation phase1/data/scz_pgc3_lead_snps_hg38.annotation.tsv]
"""

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

NEURAL_CATS = ["CNS", "neural_stem_cell_derived", "eye_retina", "PNS_neural_crest"]
NOT_BODY = set(NEURAL_CATS) | {"CNS_non_neural", "eye_anterior", "pluripotent", "other",
                                "stem_cell_derived_other"}
INK, INK2, MUTED, GRID, AXIS, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
BLUE, ORANGE, GREEN, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#898781"
CLASS_COLOR = {"brain_specific": BLUE, "shared": GREEN, "other_tissue": ORANGE, "none": GREY}


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def quantile_against(values, ref_sorted):
    """Share of reference values below each value; ties count half."""
    lo = np.searchsorted(ref_sorted, values, side="left")
    hi = np.searchsorted(ref_sorted, values, side="right")
    return (lo + hi) / 2 / len(ref_sorted)


def fmt(v):
    return "" if v is None or (isinstance(v, float) and math.isnan(v)) else (f"{v:.4f}" if isinstance(v, float) else str(v))


def main():
    ap = argparse.ArgumentParser(description="Per-track quantile normalisation and per-SNP tissue profiles.")
    ap.add_argument("--long", default="results/summaries/track_scores_long.tsv")
    ap.add_argument("--categories", default="phase2/metadata/track_categories.tsv")
    ap.add_argument("--sets", nargs="+", default=["lead"], help="score sets to profile")
    ap.add_argument("--reference", default="lead", help="score set that defines each track's distribution")
    ap.add_argument("--modality", default="DNASE")
    ap.add_argument("--annotation", help="optional SNP annotation TSV (from phase1/annotate_snps.py)")
    ap.add_argument("--include-cell-lines", action="store_true")
    ap.add_argument("--min-tracks", type=int, default=2, help="drop categories with fewer tracks")
    ap.add_argument("--high", type=float, default=0.75)
    ap.add_argument("--min-top", type=int, default=3,
                    help="categories need at least this many tracks to be the 'top body category' "
                         "(a mean over 1-2 tracks is noisy and wins the max too often)")
    ap.add_argument("--n-closeups", type=int, default=4)
    ap.add_argument("--outdir", default="results/profiles")
    args = ap.parse_args()

    for p in (args.long, args.categories):
        if not Path(p).exists():
            sys.exit(f"Not found: {p}")
    cats = {(r["output_type"], int(r["track_index"])): r for r in read_tsv(args.categories)}

    # ---- load scores for this modality
    scores = defaultdict(dict)         # (set, rsid) -> {track: (abs, risk_oriented)}
    snp_info = {}
    for r in read_tsv(args.long):
        if r["modality"] != args.modality:
            continue
        key = (r["set"], r["rsid"])
        ro = r["risk_oriented_score"]
        scores[key][int(r["track_index"])] = (abs(float(r["variant_effect_score"])),
                                              float(ro) if ro != "" else float("nan"))
        snp_info[key] = r
    have_sets = sorted({s for s, _ in scores})
    for s in set(args.sets) | {args.reference}:
        if s not in have_sets:
            sys.exit(f"Set {s!r} not in {args.long} (sets present: {have_sets})")

    # ---- tracks used: main-analysis tracks with a category, present for every SNP
    all_tracks = sorted(set.intersection(*(set(v) for v in scores.values())))
    use = []
    for t in all_tracks:
        c = cats.get((args.modality, t))
        if c is None:
            sys.exit(f"{args.modality} track {t} has no category (rerun assign_tissue_categories.py)")
        if c["main_analysis"] == "True" or args.include_cell_lines:
            use.append(t)
    track_cat = {t: cats[(args.modality, t)]["tissue_category"] for t in use}
    stage = {t: cats[(args.modality, t)]["life_stage_group"] for t in use}
    n_per_cat = defaultdict(int)
    for t in use:
        n_per_cat[track_cat[t]] += 1
    categories = [c for c in sorted(n_per_cat, key=lambda c: (c not in NEURAL_CATS,
                                                             NEURAL_CATS.index(c) if c in NEURAL_CATS else 0, c))
                  if n_per_cat[c] >= args.min_tracks and c not in ("other",)]
    if "CNS" not in categories:
        sys.exit("No CNS tracks in this modality - a brain profile is not possible (e.g. ATAC has none).")
    body = [c for c in categories if c not in NOT_BODY]
    body_top = [c for c in body if n_per_cat[c] >= args.min_top]

    # ---- per-track reference distributions and quantiles
    ref_keys = [k for k in scores if k[0] == args.reference]
    ref_sorted = {t: np.sort([scores[k][t][0] for k in ref_keys]) for t in use}
    profile_keys = [k for k in scores if k[0] in args.sets]
    q = {k: {t: float(quantile_against(np.array([scores[k][t][0]]), ref_sorted[t])[0]) for t in use}
         for k in profile_keys}

    annot = {}
    if args.annotation:
        annot = {r["rsid"]: r for r in read_tsv(args.annotation)}

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    mod = args.modality.lower()

    # ---- profiles and summary
    prof_rows, summ = [], []
    for k in profile_keys:
        sname, rs = k
        cm = {c: float(np.mean([q[k][t] for t in use if track_cat[t] == c])) for c in categories}
        brain_tr = [t for t in use if track_cat[t] == "CNS"]
        pre = [q[k][t] for t in brain_tr if stage[t] == "prenatal"]
        post = [q[k][t] for t in brain_tr if stage[t] == "postnatal"]
        top_body = max(body_top, key=lambda c: cm[c])
        body_mean = float(np.mean([cm[c] for c in body]))
        ranked = sorted(categories, key=lambda c: -cm[c])
        b, tb = cm["CNS"], cm[top_body]
        cls = ("brain_specific" if b >= args.high and tb < args.high else "shared" if b >= args.high
               else "other_tissue" if tb >= args.high else "none")
        ro = [scores[k][t][1] for t in brain_tr if not math.isnan(scores[k][t][1])]
        a = annot.get(rs, {})
        summ.append({"set": sname, "rsid": rs, "chr": snp_info[k]["chr"], "pos": snp_info[k]["pos"],
                     "nearest_gene": a.get("nearest_gene", ""), "consequence": a.get("most_severe_consequence", ""),
                     "brain_q": b, "brain_prenatal_q": float(np.mean(pre)) if pre else float("nan"),
                     "brain_postnatal_q": float(np.mean(post)) if post else float("nan"),
                     "body_mean_q": body_mean, "top_body_category": top_body, "top_body_q": tb,
                     "specificity": b - body_mean, "brain_rank": ranked.index("CNS") + 1,
                     "n_categories": len(categories), "profile_class": cls,
                     "max_abs_score_brain": max(scores[k][t][0] for t in brain_tr),
                     "max_abs_score_any": max(scores[k][t][0] for t in use),
                     "mean_risk_oriented_brain": float(np.mean(ro)) if ro else float("nan"),
                     "palindromic": snp_info[k]["palindromic"]})
        prof_rows.append({"set": sname, "rsid": rs, **{c: cm[c] for c in categories}})
    summ.sort(key=lambda r: (r["set"] != args.sets[0], -r["specificity"]))
    order = {(r["set"], r["rsid"]): i for i, r in enumerate(summ)}
    prof_rows.sort(key=lambda r: order[(r["set"], r["rsid"])])

    def write(path, rows):
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows({k: fmt(v) for k, v in r.items()} for r in rows)
    write(out / f"{mod}_snp_summary.tsv", summ)
    write(out / f"{mod}_snp_category_profiles.tsv", prof_rows)
    with open(out / f"{mod}_track_quantiles.tsv", "w", encoding="utf-8") as f:
        f.write("set\trsid\ttrack_index\ttissue_category\tbiosample_name\tquantile\tabs_score\n")
        for k in profile_keys:
            for t in use:
                f.write(f"{k[0]}\t{k[1]}\t{t}\t{track_cat[t]}\t{cats[(args.modality, t)]['biosample_name']}"
                        f"\t{q[k][t]:.4f}\t{scores[k][t][0]:.6f}\n")
    with open(out / f"{mod}_categories_used.tsv", "w", encoding="utf-8") as f:
        f.write("tissue_category\tn_tracks\tin_body_mean\n")
        for c in categories:
            f.write(f"{c}\t{n_per_cat[c]}\t{c in body}\n")

    # ---- console
    main_set = [r for r in summ if r["set"] == args.sets[0]]
    print(f"{args.modality}: {len(use)} tracks in {len(categories)} categories "
          f"({', '.join(f'{c} {n_per_cat[c]}' for c in categories)})")
    print(f"Reference for quantiles: set '{args.reference}' ({len(ref_keys)} SNPs)"
          + ("  [INTERIM: profiles are relative to the other SCZ leads]" if args.reference in args.sets else ""))
    cnt = defaultdict(int)
    for r in main_set:
        cnt[r["profile_class"]] += 1
    print(f"Profile classes ({args.sets[0]}, threshold {args.high}): " +
          ", ".join(f"{c} {cnt[c]}" for c in ("brain_specific", "shared", "other_tissue", "none")))
    lab = lambda r: r["rsid"] + (f" ({r['nearest_gene']})" if r["nearest_gene"] else "")
    print("\nMost brain-specific:")
    for r in main_set[:10]:
        print(f"  {lab(r):<28} brain {r['brain_q']:.2f}  body mean {r['body_mean_q']:.2f}  "
              f"top body {r['top_body_category']} {r['top_body_q']:.2f}  [{r['profile_class']}]")
    print("Least brain (most other-tissue):")
    for r in main_set[-5:]:
        print(f"  {lab(r):<28} brain {r['brain_q']:.2f}  body mean {r['body_mean_q']:.2f}  "
              f"top body {r['top_body_category']} {r['top_body_q']:.2f}  [{r['profile_class']}]")

    plots(summ, prof_rows, categories, q, use, track_cat, args, out, mod)
    print(f"\nWritten to {out}/")
    return 0


# ============================================================================ figures
def plots(summ, prof_rows, categories, q, use, track_cat, args, out, mod):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    plt.rcParams.update({"font.size": 8, "axes.edgecolor": AXIS, "xtick.color": INK2, "ytick.color": INK2,
                         "axes.labelcolor": INK2, "figure.facecolor": SURF, "axes.facecolor": SURF,
                         "savefig.facecolor": SURF, "axes.spines.top": False, "axes.spines.right": False})
    main_set = [r for r in summ if r["set"] == args.sets[0]]
    others = [r for r in summ if r["set"] != args.sets[0]]
    label = lambda r: r["rsid"] + (f" {r['nearest_gene']}" if r["nearest_gene"] else "")
    interim = args.reference in args.sets
    ref_note = (f"Quantiles within each track relative to the {args.reference} SNPs (interim: not yet a common-variant null)"
                if interim else f"Quantiles within each track relative to the '{args.reference}' null SNPs")

    # 1. heatmap
    rows = [r for r in prof_rows if r["set"] == args.sets[0]]
    M = np.array([[r[c] for c in categories] for r in rows])
    cmap = LinearSegmentedColormap.from_list("seq", ["#f7f6f2", "#bcd6f3", "#2a78d6", "#0d3a73"])
    h = max(4, 0.12 * len(rows) + 1.8)
    fig, ax = plt.subplots(figsize=(0.42 * len(categories) + 3.2, h))
    im = ax.imshow(M, aspect="auto", cmap=cmap, vmin=0, vmax=1, interpolation="nearest")
    ax.set_xticks(range(len(categories)))
    ax.set_xticklabels([c.replace("_", " ") for c in categories], rotation=60, ha="right")
    ax.set_yticks(range(len(rows)))
    lab = {(r["set"], r["rsid"]): label(r) for r in summ}
    ax.set_yticklabels([lab[(r["set"], r["rsid"])] for r in rows], fontsize=max(3, min(7, 400 / len(rows))))
    for x in (0.5, len([c for c in categories if c in NEURAL_CATS]) - 0.5):
        ax.axvline(x, color="white", linewidth=1.5)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("mean within-track quantile")
    ax.set_title(f"{args.modality}: tissue profile of each SNP (sorted: most brain-specific at top)",
                 loc="left", fontsize=9, color=INK)
    fig.text(0.01, 0.002, ref_note, fontsize=6.5, color=INK2)
    fig.tight_layout(rect=(0, 0.01, 1, 1))
    fig.savefig(out / f"{mod}_heatmap.png", dpi=200)
    plt.close(fig)

    # 2. slope plot: body mean -> brain
    fig, ax = plt.subplots(figsize=(5.2, 6.4))
    d = sorted(main_set, key=lambda r: r["specificity"])
    for r in d:
        up = r["specificity"] > 0.1
        down = r["specificity"] < -0.1
        col = BLUE if up else ORANGE if down else GRID
        ax.plot([0, 1], [r["body_mean_q"], r["brain_q"]], color=col, linewidth=1.6 if (up or down) else 0.8,
                alpha=0.9 if (up or down) else 0.7, zorder=3 if (up or down) else 2)
    lab_pts = sorted(d[-6:] + d[:4], key=lambda r: r["brain_q"])
    ys, gap = [], 0.028
    for r in lab_pts:                      # spread labels so they don't overlap
        y = r["brain_q"] if not ys else max(r["brain_q"], ys[-1] + gap)
        ys.append(y)
    shift = max(0, ys[-1] - 0.99) if ys else 0
    for r, y in zip(lab_pts, ys):
        ax.plot([1.0, 1.03], [r["brain_q"], y - shift], color=AXIS, linewidth=0.5)
        ax.text(1.04, y - shift, label(r), fontsize=6.5, va="center",
                color=BLUE if r["specificity"] > 0 else ORANGE)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["non-neural body tissues\n(mean of categories)", "brain (CNS)"])
    ax.set_xlim(-0.15, 1.55)
    ax.set_ylim(0, 1)
    ax.set_ylabel("mean within-track quantile")
    ax.grid(True, axis="y", color=GRID, linewidth=0.6)
    n_up = sum(r["specificity"] > 0.1 for r in d)
    n_dn = sum(r["specificity"] < -0.1 for r in d)
    ax.set_title(f"{args.sets[0]} SNPs: {n_up} higher in brain, {n_dn} lower (by > 0.1)",
                 loc="left", fontsize=9, color=INK)
    fig.text(0.01, 0.005, ref_note, fontsize=6, color=INK2)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    fig.savefig(out / f"{mod}_slope.png", dpi=200)
    plt.close(fig)

    # 3. scatter: brain vs top body category
    fig, ax = plt.subplots(figsize=(6.2, 5.6))
    ax.axhline(args.high, color=AXIS, linewidth=1)
    ax.axvline(args.high, color=AXIS, linewidth=1)
    for cls, col in CLASS_COLOR.items():
        pts = [r for r in main_set if r["profile_class"] == cls]
        ax.scatter([r["top_body_q"] for r in pts], [r["brain_q"] for r in pts], s=22, color=col,
                   edgecolor="white", linewidth=0.5, label=f"{cls.replace('_', ' ')} ({len(pts)})", zorder=3)
    if others:
        ax.scatter([r["top_body_q"] for r in others], [r["brain_q"] for r in others], s=24, marker="^",
                   facecolor="none", edgecolor=INK2, linewidth=0.8, label=f"{others[0]['set']} SNPs", zorder=4)
    for j, r in enumerate(sorted(main_set, key=lambda r: -r["specificity"])[:5]):
        ax.annotate(label(r), (r["top_body_q"], r["brain_q"]), fontsize=6, color=INK2,
                    xytext=(-70, 10 + 9 * j), textcoords="offset points",
                    arrowprops=dict(arrowstyle="-", color=AXIS, linewidth=0.5))
    for j, r in enumerate(sorted(main_set, key=lambda r: r["specificity"])[:3]):
        ax.annotate(f"{label(r)} [{r['top_body_category']}]", (r["top_body_q"], r["brain_q"]), fontsize=6,
                    color=ORANGE, xytext=(-150, -12 - 9 * j), textcoords="offset points",
                    arrowprops=dict(arrowstyle="-", color=AXIS, linewidth=0.5))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("strongest non-neural tissue category (mean quantile)")
    ax.set_ylabel("brain / CNS (mean quantile)")
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    ax.set_title("Impactful in brain? vs impactful elsewhere?", loc="left", fontsize=9, color=INK)
    fig.text(0.01, 0.005, ref_note, fontsize=6, color=INK2)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    fig.savefig(out / f"{mod}_scatter.png", dpi=200)
    plt.close(fig)

    # 4. close-ups: most brain-high and most brain-low
    n = args.n_closeups
    picks = main_set[:n] + main_set[-n:][::-1]
    fig, axes = plt.subplots(2, n, figsize=(2.9 * n, 6.2), sharey=True, squeeze=False)
    rng = np.random.default_rng(0)
    for ax, r in zip(axes.flat, picks):
        k = (r["set"], r["rsid"])
        for i, c in enumerate(categories):
            ys = [q[k][t] for t in use if track_cat[t] == c]
            xs = i + rng.uniform(-0.22, 0.22, len(ys))
            ax.scatter(xs, ys, s=7, color=BLUE if c == "CNS" else (ORANGE if c in NEURAL_CATS else GREY),
                       alpha=0.85, linewidth=0, zorder=3)
            ax.plot([i - 0.3, i + 0.3], [np.mean(ys)] * 2, color=INK, linewidth=1.4, zorder=4)
        ax.set_xticks(range(len(categories)))
        ax.set_xticklabels([c.replace("_", " ") for c in categories], rotation=90, fontsize=5.5)
        ax.set_ylim(-0.02, 1.02)
        ax.grid(True, axis="y", color=GRID, linewidth=0.5)
        ax.set_title(f"{label(r)}\nbrain {r['brain_q']:.2f} · body {r['body_mean_q']:.2f}",
                     loc="left", fontsize=7.5, color=INK)
    axes[0, 0].set_ylabel("within-track quantile\n(most brain-specific)")
    axes[1, 0].set_ylabel("within-track quantile\n(least brain-specific)")
    fig.text(0.01, 0.004, ref_note + ".  Blue = CNS, orange = other neural, grey = body; bar = category mean.",
             fontsize=6.5, color=INK2)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    fig.savefig(out / f"{mod}_closeups.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())