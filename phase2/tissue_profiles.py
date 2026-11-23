#!/usr/bin/env python3
"""
tissue_profiles.py
------------------
Per-SNP tissue profiles from per-track NORMALISED scores (Choo, Nov 5 and Nov 9).

The problem it solves
    Raw variant scores are not comparable between tracks: each track has its own
    dynamic range (AlphaGenome docs: scores from "different tracks and modalities
    ... are on different scales"). So each score is first normalised WITHIN ITS OWN
    TRACK against a reference set of variants scored by the same model and code.
    Two methods (--method):

    ratio  (control-GWAS denominator; Choo, Nov 9; use this for the preliminary results)
        D_t     = mean |score| of the reference (control-trait) SNPs in track t
        n(i, t) = risk_oriented_score(i, t) / D_t
        |n| = 1 means "as large as a typical control-trait lead in this track";
        the sign is the direction (+ = the risk allele opens chromatin).
        Reference = the control GWAS leads (e.g. --reference t1d). Later upgrade:
        per-SNP matched controls (1-5 per SNP, matched on EUR allele frequency and
        distance to the nearest gene) instead of one denominator per track.

    quantile  (Nov 5 interim method; unchanged)
        q(i, t) = share of reference variants whose |score| in track t is smaller
                  (ties count half). 0-1 in every track.
        With a small reference set the quantiles are coarse (20 SNPs -> steps of
        0.05, and anything above all of them is 1.0), so ratio is preferred there.

Per SNP x tissue category (main-analysis tracks: tissues, primary cells, stem-cell-
derived; cell lines excluded unless --include-cell-lines), each category weighted
equally. Categories from phase2/metadata/track_categories.tsv.
    size         mean over the category's tracks of |n|   (ratio)  or q (quantile)
    direction    mean signed n over the category's tracks  (ratio only)
    concordance  share of the category's tracks whose n has the majority sign (ratio only)

Per SNP summary
    brain_size (brain_q in quantile mode), prenatal / postnatal split
    body_mean    mean size over the non-neural body categories (equal weight)
    top_body_category, top_body_size   body category with the largest size
                 (only categories with >= --min-top tracks, default 3)
    specificity  brain size - body mean
    brain_rank   rank of CNS among all categories (1 = largest)
    profile_class, threshold --high (default 1.0 for ratio, 0.75 for quantile):
                   brain_specific  brain >= high, top body <  high
                   shared          brain >= high, top body >= high
                   other_tissue    brain <  high, top body >= high
                   none            neither
    ratio only:  brain_direction, brain_concordance, top_body_direction,
                 bootstrap 95% intervals for brain_size and specificity
                 (reference SNPs resampled with replacement, --bootstrap times)

Guards (ratio)
    Tracks whose denominator D_t is in the bottom --low-pct percent (default 5) are
    flagged: a tiny denominator makes ratios explode. --drop-low-denominator removes
    them; run with and without and compare. Tracks with D_t = 0 are always dropped.
    Palindromic SNPs keep their size; their direction assumes the + strand (flagged).

Outputs (results/profiles/ by default; file names carry the method)
    <mod>_<method>_snp_summary.tsv            one row per SNP
    <mod>_<method>_snp_category_profiles.tsv  size per SNP x category
    <mod>_ratio_snp_category_direction.tsv    direction and concordance per SNP x category
    <mod>_<method>_track_values.tsv           per SNP x track value
    <mod>_<method>_categories_used.tsv        categories, track counts, body or not
    <mod>_ratio_denominators.tsv              D_t per track, category, low flag
    <mod>_ratio_denominator_by_category.tsv   D_t summary per category (control-trait bias check)
    figures: <mod>_<method>_heatmap.png, _slope.png, _scatter.png, _closeups.png,
             <mod>_ratio_denominators.png

Usage (from the repo root)
    # control-GWAS ratio (preliminary results): pilots, with random SCZ SNPs as a check
    python phase2/tissue_profiles.py --method ratio --sets pilot control --reference t1d \\
        --annotation phase1/data/scz_pgc3_lead_snps_hg38.annotation.tsv
    # Nov 5 interim quantiles
    python phase2/tissue_profiles.py --method quantile --sets lead control --reference lead
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


def nanmean(xs):
    xs = [x for x in xs if not math.isnan(x)]
    return float(np.mean(xs)) if xs else float("nan")


def concordance(xs):
    """Share of values with the majority sign (zeros and NaN ignored)."""
    xs = [x for x in xs if not math.isnan(x) and x != 0]
    if not xs:
        return float("nan")
    pos = sum(x > 0 for x in xs)
    return max(pos, len(xs) - pos) / len(xs)


def main():
    ap = argparse.ArgumentParser(description="Per-track normalisation and per-SNP tissue profiles.")
    ap.add_argument("--method", choices=["ratio", "quantile"], default="quantile",
                    help="ratio = control-GWAS denominator (Nov 9); quantile = Nov 5 interim (default, unchanged)")
    ap.add_argument("--long", default="results/summaries/track_scores_long.tsv")
    ap.add_argument("--categories", default="phase2/metadata/track_categories.tsv")
    ap.add_argument("--sets", nargs="+", default=["lead"], help="score sets to profile")
    ap.add_argument("--reference", default="lead", help="score set that defines each track's scale")
    ap.add_argument("--modality", default="DNASE")
    ap.add_argument("--annotation", help="optional SNP annotation TSV (from phase1/annotate_snps.py)")
    ap.add_argument("--include-cell-lines", action="store_true")
    ap.add_argument("--min-tracks", type=int, default=2, help="drop categories with fewer tracks")
    ap.add_argument("--high", type=float, default=None,
                    help="class threshold (default 1.0 for ratio, 0.75 for quantile)")
    ap.add_argument("--min-top", type=int, default=3,
                    help="categories need at least this many tracks to be the 'top body category' "
                         "(a mean over 1-2 tracks is noisy and wins the max too often)")
    ap.add_argument("--low-pct", type=float, default=5.0,
                    help="ratio: flag tracks whose denominator is in this bottom percentile")
    ap.add_argument("--drop-low-denominator", action="store_true",
                    help="ratio: leave the flagged low-denominator tracks out")
    ap.add_argument("--bootstrap", type=int, default=1000, help="ratio: bootstrap resamples (0 = off)")
    ap.add_argument("--seed", type=int, default=2026, help="ratio: bootstrap seed")
    ap.add_argument("--n-closeups", type=int, default=4)
    ap.add_argument("--outdir", default="results/profiles")
    args = ap.parse_args()
    ratio = args.method == "ratio"
    high = args.high if args.high is not None else (1.0 if ratio else 0.75)

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
    if ratio and args.reference in args.sets:
        sys.exit("--method ratio: the reference (control) set must not also be profiled; "
                 "every track would average 1 for it by construction.")

    # ---- tracks used: main-analysis tracks with a category, present for every SNP
    all_tracks = sorted(set.intersection(*(set(v) for v in scores.values())))
    use = []
    for t in all_tracks:
        c = cats.get((args.modality, t))
        if c is None:
            sys.exit(f"{args.modality} track {t} has no category (rerun assign_tissue_categories.py)")
        if c["main_analysis"] == "True" or args.include_cell_lines:
            use.append(t)

    ref_keys = sorted(k for k in scores if k[0] == args.reference)
    profile_keys = [k for k in scores if k[0] in args.sets]

    # ---- per-track normalisation
    denom, low = {}, set()
    if ratio:
        denom = {t: float(np.mean([scores[k][t][0] for k in ref_keys])) for t in use}
        zero = [t for t in use if denom[t] == 0]
        positive = [denom[t] for t in use if denom[t] > 0]
        cut = float(np.percentile(positive, args.low_pct)) if positive else 0.0
        low = {t for t in use if 0 < denom[t] <= cut}
        drop = set(zero) | (low if args.drop_low_denominator else set())
        all_use = list(use)
        use = [t for t in use if t not in drop]
        size = {k: {t: scores[k][t][0] / denom[t] for t in use} for k in profile_keys}
        signed = {k: {t: scores[k][t][1] / denom[t] for t in use} for k in profile_keys}
    else:
        ref_sorted = {t: np.sort([scores[k][t][0] for k in ref_keys]) for t in use}
        size = {k: {t: float(quantile_against(np.array([scores[k][t][0]]), ref_sorted[t])[0]) for t in use}
                for k in profile_keys}
        signed = None

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

    annot = {}
    if args.annotation:
        annot = {r["rsid"]: r for r in read_tsv(args.annotation)}

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    mod = args.modality.lower()
    stem = f"{mod}_{args.method}"

    # ---- bootstrap (ratio): resample reference SNPs, recompute denominators
    boot = {}
    if ratio and args.bootstrap > 0:
        rng = np.random.default_rng(args.seed)
        R = np.array([[scores[k][t][0] for t in use] for k in ref_keys])                 # n_ref x T
        S = np.array([[scores[k][t][0] for t in use] for k in profile_keys])             # n_snp x T
        idx = rng.integers(0, len(ref_keys), size=(args.bootstrap, len(ref_keys)))
        D = R[idx].mean(axis=1)                                                          # B x T
        D[D == 0] = np.nan
        masks = {c: np.array([track_cat[t] == c for t in use]) for c in categories}
        for i, k in enumerate(profile_keys):
            v = S[i][None, :] / D                                                        # B x T
            cm = {c: np.nanmean(v[:, masks[c]], axis=1) for c in categories}
            brain = cm["CNS"]
            spec = brain - np.mean([cm[c] for c in body], axis=0)
            boot[k] = (np.nanpercentile(brain, [2.5, 97.5]), np.nanpercentile(spec, [2.5, 97.5]))

    # ---- profiles and summary
    prof_rows, dir_rows, summ = [], [], []
    for k in profile_keys:
        sname, rs = k
        cm = {c: float(np.mean([size[k][t] for t in use if track_cat[t] == c])) for c in categories}
        brain_tr = [t for t in use if track_cat[t] == "CNS"]
        pre = [size[k][t] for t in brain_tr if stage[t] == "prenatal"]
        post = [size[k][t] for t in brain_tr if stage[t] == "postnatal"]
        top_body = max(body_top, key=lambda c: cm[c])
        body_mean = float(np.mean([cm[c] for c in body]))
        ranked = sorted(categories, key=lambda c: -cm[c])
        b, tb = cm["CNS"], cm[top_body]
        cls = ("brain_specific" if b >= high and tb < high else "shared" if b >= high
               else "other_tissue" if tb >= high else "none")
        ro = [scores[k][t][1] for t in brain_tr if not math.isnan(scores[k][t][1])]
        a = annot.get(rs, {})
        row = {"set": sname, "rsid": rs, "chr": snp_info[k]["chr"], "pos": snp_info[k]["pos"],
               "nearest_gene": a.get("nearest_gene", ""), "consequence": a.get("most_severe_consequence", "")}
        if ratio:
            cd = {c: nanmean([signed[k][t] for t in use if track_cat[t] == c]) for c in categories}
            cc = {c: concordance([signed[k][t] for t in use if track_cat[t] == c]) for c in categories}
            ci = boot.get(k)
            row.update({"brain_size": b,
                        "brain_size_ci_lo": ci[0][0] if ci else float("nan"),
                        "brain_size_ci_hi": ci[0][1] if ci else float("nan"),
                        "brain_direction": cd["CNS"], "brain_concordance": cc["CNS"],
                        "brain_prenatal_size": float(np.mean(pre)) if pre else float("nan"),
                        "brain_postnatal_size": float(np.mean(post)) if post else float("nan"),
                        "body_mean_size": body_mean, "top_body_category": top_body, "top_body_size": tb,
                        "top_body_direction": cd[top_body], "top_body_concordance": cc[top_body],
                        "specificity": b - body_mean,
                        "specificity_ci_lo": ci[1][0] if ci else float("nan"),
                        "specificity_ci_hi": ci[1][1] if ci else float("nan")})
            dir_rows.append({"set": sname, "rsid": rs,
                             **{f"{c}_direction": cd[c] for c in categories},
                             **{f"{c}_concordance": cc[c] for c in categories}})
        else:
            row.update({"brain_q": b, "brain_prenatal_q": float(np.mean(pre)) if pre else float("nan"),
                        "brain_postnatal_q": float(np.mean(post)) if post else float("nan"),
                        "body_mean_q": body_mean, "top_body_category": top_body, "top_body_q": tb,
                        "specificity": b - body_mean})
        row.update({"brain_rank": ranked.index("CNS") + 1, "n_categories": len(categories),
                    "profile_class": cls,
                    "max_abs_score_brain": max(scores[k][t][0] for t in brain_tr),
                    "max_abs_score_any": max(scores[k][t][0] for t in use),
                    "mean_risk_oriented_brain": float(np.mean(ro)) if ro else float("nan"),
                    "palindromic": snp_info[k]["palindromic"]})
        if ratio:
            row["direction_note"] = ("strand-assumed (palindromic)" if snp_info[k]["palindromic"] == "True"
                                     else "")
        summ.append(row)
        prof_rows.append({"set": sname, "rsid": rs, **{c: cm[c] for c in categories}})
    summ.sort(key=lambda r: (r["set"] != args.sets[0], -r["specificity"]))
    order = {(r["set"], r["rsid"]): i for i, r in enumerate(summ)}
    prof_rows.sort(key=lambda r: order[(r["set"], r["rsid"])])
    dir_rows.sort(key=lambda r: order[(r["set"], r["rsid"])])

    def write(path, rows):
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows({k: fmt(v) for k, v in r.items()} for r in rows)
    write(out / f"{stem}_snp_summary.tsv", summ)
    write(out / f"{stem}_snp_category_profiles.tsv", prof_rows)
    if ratio:
        write(out / f"{stem}_snp_category_direction.tsv", dir_rows)
    vname = "ratio" if ratio else "quantile"
    with open(out / f"{stem}_track_values.tsv", "w", encoding="utf-8") as f:
        f.write(f"set\trsid\ttrack_index\ttissue_category\tbiosample_name\t{vname}"
                + ("\tsigned_ratio" if ratio else "") + "\tabs_score\n")
        for k in profile_keys:
            for t in use:
                f.write(f"{k[0]}\t{k[1]}\t{t}\t{track_cat[t]}\t{cats[(args.modality, t)]['biosample_name']}"
                        f"\t{size[k][t]:.4f}" + (f"\t{fmt(signed[k][t])}" if ratio else "")
                        + f"\t{scores[k][t][0]:.6f}\n")
    with open(out / f"{stem}_categories_used.tsv", "w", encoding="utf-8") as f:
        f.write("tissue_category\tn_tracks\tin_body_mean\n")
        for c in categories:
            f.write(f"{c}\t{n_per_cat[c]}\t{c in body}\n")
    if ratio:
        with open(out / f"{mod}_ratio_denominators.tsv", "w", encoding="utf-8") as f:
            f.write("track_index\ttissue_category\tbiosample_name\tdenominator\tlow_flag\tused\n")
            for t in all_use:
                c = cats[(args.modality, t)]
                f.write(f"{t}\t{c['tissue_category']}\t{c['biosample_name']}\t{denom[t]:.6f}"
                        f"\t{t in low}\t{t in use}\n")
        by_cat = defaultdict(list)
        for t in all_use:
            by_cat[cats[(args.modality, t)]["tissue_category"]].append(denom[t])
        overall = float(np.median([denom[t] for t in all_use]))
        with open(out / f"{mod}_ratio_denominator_by_category.tsv", "w", encoding="utf-8") as f:
            f.write("tissue_category\tn_tracks\tmedian_denominator\tmean_denominator\t"
                    "median_vs_all_tracks\tn_low_flag\n")
            for c in sorted(by_cat, key=lambda c: -np.median(by_cat[c])):
                v = by_cat[c]
                nl = sum(1 for t in all_use if cats[(args.modality, t)]["tissue_category"] == c and t in low)
                f.write(f"{c}\t{len(v)}\t{np.median(v):.6f}\t{np.mean(v):.6f}\t{np.median(v) / overall:.3f}\t{nl}\n")

    # ---- console
    main_set = [r for r in summ if r["set"] == args.sets[0]]
    print(f"{args.modality}: {len(use)} tracks in {len(categories)} categories "
          f"({', '.join(f'{c} {n_per_cat[c]}' for c in categories)})")
    print(f"Method: {args.method}. Reference set '{args.reference}' ({len(ref_keys)} SNPs)"
          + ("  [INTERIM: profiles are relative to the other SCZ leads]"
             if (not ratio and args.reference in args.sets) else ""))
    if ratio:
        print(f"Denominators: median {np.median([denom[t] for t in all_use]):.4f}; "
              f"{len(low)} tracks flagged low (bottom {args.low_pct:g}%)"
              + ("; dropped" if args.drop_low_denominator else "; kept (use --drop-low-denominator to compare)")
              + (f"; {len(zero)} zero-denominator tracks dropped" if zero else ""))
        print(f"Bootstrap: {args.bootstrap} resamples of the reference SNPs, seed {args.seed}")
    cnt = defaultdict(int)
    for r in main_set:
        cnt[r["profile_class"]] += 1
    print(f"Profile classes ({args.sets[0]}, threshold {high}): " +
          ", ".join(f"{c} {cnt[c]}" for c in ("brain_specific", "shared", "other_tissue", "none")))
    lab = lambda r: r["rsid"] + (f" ({r['nearest_gene']})" if r["nearest_gene"] else "")
    bkey, bodyk, topk = (("brain_size", "body_mean_size", "top_body_size") if ratio
                         else ("brain_q", "body_mean_q", "top_body_q"))
    for title, rows_ in (("Most brain-specific:", main_set[:10]), ("Least brain (most other-tissue):", main_set[-5:])):
        print("\n" + title if title.startswith("Most") else title)
        for r in rows_:
            extra = (f"  dir {r['brain_direction']:+.2f} conc {r['brain_concordance']:.2f}"
                     if ratio and not math.isnan(r["brain_direction"]) else "")
            print(f"  {lab(r):<28} brain {r[bkey]:.2f}  body mean {r[bodyk]:.2f}  "
                  f"top body {r['top_body_category']} {r[topk]:.2f}  [{r['profile_class']}]{extra}")

    plots(summ, prof_rows, dir_rows, categories, size, signed, use, track_cat, args, out, stem, high,
          (denom, low, all_use, cats) if ratio else None)
    print(f"\nWritten to {out}/")
    return 0


# ============================================================================ figures
def plots(summ, prof_rows, dir_rows, categories, size, signed, use, track_cat, args, out, stem, high, den):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
    plt.rcParams.update({"font.size": 8, "axes.edgecolor": AXIS, "xtick.color": INK2, "ytick.color": INK2,
                         "axes.labelcolor": INK2, "figure.facecolor": SURF, "axes.facecolor": SURF,
                         "savefig.facecolor": SURF, "axes.spines.top": False, "axes.spines.right": False})
    ratio = args.method == "ratio"
    main_set = [r for r in summ if r["set"] == args.sets[0]]
    others = [r for r in summ if r["set"] != args.sets[0]]
    label = lambda r: r["rsid"] + (f" {r['nearest_gene']}" if r["nearest_gene"] else "")
    if ratio:
        bk, bodyk, topk = "brain_size", "body_mean_size", "top_body_size"
        unit = "size = mean |score| ÷ control mean |score| (1 = typical control lead)"
        ref_note = (f"Each track's scores divided by the mean |score| of the reference set '{args.reference}' "
                    f"in that track ({'low denominators dropped' if args.drop_low_denominator else 'all tracks kept'})")
    else:
        bk, bodyk, topk = "brain_q", "body_mean_q", "top_body_q"
        unit = "mean within-track quantile"
        interim = args.reference in args.sets
        ref_note = (f"Quantiles within each track relative to the {args.reference} SNPs (interim: not yet a common-variant null)"
                    if interim else f"Quantiles within each track relative to the '{args.reference}' SNPs")
    vmax_size = 1.0 if not ratio else max(high * 1.5, max(max(r[bk], r[topk]) for r in summ) * 1.08)

    # 1. heatmap
    rows = [r for r in prof_rows if r["set"] == args.sets[0]]
    lab = {(r["set"], r["rsid"]): label(r) for r in summ}
    h = max(4, 0.12 * len(rows) + 1.8) if len(rows) > 30 else max(3.6, 0.32 * len(rows) + 1.9)
    fig, ax = plt.subplots(figsize=(0.42 * len(categories) + 3.6, h))
    if ratio:
        drows = {(r["set"], r["rsid"]): r for r in dir_rows}
        M = np.array([[drows[(r["set"], r["rsid"])][f"{c}_direction"] for c in categories] for r in rows])
        S = np.array([[r[c] for c in categories] for r in rows])
        lim = max(1.0, float(np.nanpercentile(np.abs(M), 95)))   # colour capped; numbers show the true size
        cmap = LinearSegmentedColormap.from_list("div", [ORANGE, "#f7d3c2", "#f7f6f2", "#bcd6f3", BLUE])
        im = ax.imshow(M, aspect="auto", cmap=cmap, norm=TwoSlopeNorm(0, -lim, lim), interpolation="nearest")
        if len(rows) <= 40:
            for i in range(S.shape[0]):
                for j in range(S.shape[1]):
                    ax.text(j, i, f"{S[i, j]:.1f}", ha="center", va="center", fontsize=5.2,
                            color=INK if abs(M[i, j]) < 0.6 * lim else "white")
        cbl = (f"direction: mean signed ratio (+ risk allele opens, − closes)\n"
               f"colour capped at ±{lim:.1f} (95th percentile); numbers: size (mean |ratio|)")
        title = f"{args.modality}: direction and size by tissue category (sorted: most brain-specific at top)"
    else:
        M = np.array([[r[c] for c in categories] for r in rows])
        cmap = LinearSegmentedColormap.from_list("seq", ["#f7f6f2", "#bcd6f3", "#2a78d6", "#0d3a73"])
        im = ax.imshow(M, aspect="auto", cmap=cmap, vmin=0, vmax=1, interpolation="nearest")
        cbl = unit
        title = f"{args.modality}: tissue profile of each SNP (sorted: most brain-specific at top)"
    ax.set_xticks(range(len(categories)))
    ax.set_xticklabels([c.replace("_", " ") for c in categories], rotation=60, ha="right")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([lab[(r["set"], r["rsid"])] for r in rows], fontsize=max(3, min(7, 400 / len(rows))))
    for x in (0.5, len([c for c in categories if c in NEURAL_CATS]) - 0.5):
        ax.axvline(x, color="white", linewidth=1.5)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label(cbl, fontsize=6.5)
    ax.set_title(title, loc="left", fontsize=9, color=INK)
    fig.text(0.01, 0.002, ref_note, fontsize=6.0, color=INK2)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    fig.savefig(out / f"{stem}_heatmap.png", dpi=200)
    plt.close(fig)

    # 2. slope plot: body mean -> brain
    fig, ax = plt.subplots(figsize=(5.2, 6.4))
    d = sorted(main_set, key=lambda r: r["specificity"])
    fold = 1.5                                  # ratio: brain vs body mean, as a fold change
    def is_up(r):
        return (r[bk] >= fold * r[bodyk]) if ratio else r["specificity"] > 0.1
    def is_down(r):
        return (r[bk] <= r[bodyk] / fold) if ratio else r["specificity"] < -0.1
    for r in d:
        up, down = is_up(r), is_down(r)
        col = BLUE if up else ORANGE if down else GRID
        ax.plot([0, 1], [r[bodyk], r[bk]], color=col, linewidth=1.6 if (up or down) else 0.8,
                alpha=0.9 if (up or down) else 0.7, zorder=3 if (up or down) else 2)
    lab_pts = sorted(d[-6:] + d[:4] if len(d) > 10 else d, key=lambda r: r[bk])
    if ratio:                                  # ratios are multiplicative: log axis, spacing in log10
        vals_all = [v for r in d for v in (r[bk], r[bodyk]) if v > 0]
        ylo, yhi = min(vals_all) / 1.6, max(vals_all) * 1.6
        tf, inv = math.log10, (lambda v: 10 ** v)
        gap = 0.035 * (math.log10(yhi) - math.log10(ylo))
    else:
        ylo, yhi = 0, 1.0
        tf, inv = (lambda v: v), (lambda v: v)
        gap = 0.028
    ys = []
    for r in lab_pts:
        y = tf(r[bk]) if not ys else max(tf(r[bk]), ys[-1] + gap)
        ys.append(y)
    shift = max(0, ys[-1] - (tf(yhi) - gap)) if ys else 0
    for r, y in zip(lab_pts, ys):
        ax.plot([1.0, 1.03], [r[bk], inv(y - shift)], color=AXIS, linewidth=0.5)
        ax.text(1.04, inv(y - shift), label(r), fontsize=6.5, va="center",
                color=BLUE if r["specificity"] > 0 else ORANGE)
    if ratio:
        ax.set_yscale("log")
        ax.axhline(1.0, color=AXIS, linewidth=0.8, linestyle="--")
        ax.text(-0.13, 1.0, "typical control lead", fontsize=6, color=INK2, va="bottom")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["non-neural body tissues\n(mean of categories)", "brain (CNS)"])
    ax.set_xlim(-0.15, 1.55)
    ax.set_ylim(ylo, yhi)
    ax.set_ylabel(unit + (" (log scale)" if ratio else ""))
    ax.grid(True, axis="y", color=GRID, linewidth=0.6)
    n_up = sum(is_up(r) for r in d)
    n_dn = sum(is_down(r) for r in d)
    crit = f"brain ≥ {fold:g}× body mean, or ≤ 1/{fold:g}×" if ratio else "by > 0.1"
    ax.set_title(f"{args.sets[0]} SNPs: {n_up} higher in brain, {n_dn} lower\n({crit})",
                 loc="left", fontsize=9, color=INK)
    fig.text(0.01, 0.005, ref_note, fontsize=5.5, color=INK2)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    fig.savefig(out / f"{stem}_slope.png", dpi=200)
    plt.close(fig)

    # 3. scatter: brain vs top body category
    fig, ax = plt.subplots(figsize=(6.2, 5.6))
    ax.axhline(high, color=AXIS, linewidth=1)
    ax.axvline(high, color=AXIS, linewidth=1)
    for cls, col in CLASS_COLOR.items():
        pts = [r for r in main_set if r["profile_class"] == cls]
        ax.scatter([r[topk] for r in pts], [r[bk] for r in pts], s=22, color=col,
                   edgecolor="white", linewidth=0.5, label=f"{cls.replace('_', ' ')} ({len(pts)})", zorder=3)
    if others:
        ax.scatter([r[topk] for r in others], [r[bk] for r in others], s=24, marker="^",
                   facecolor="none", edgecolor=INK2, linewidth=0.8, label=f"{others[0]['set']} SNPs", zorder=4)
    if len(main_set) <= 12:                    # few SNPs: label every point next to it
        for j, r in enumerate(sorted(main_set, key=lambda r: r[bk])):   # alternate above/below
            ax.annotate(label(r) + ("" if r["profile_class"] != "other_tissue" else f" [{r['top_body_category']}]"),
                        (r[topk], r[bk]), fontsize=6, color=ORANGE if r["specificity"] < 0 else INK2,
                        xytext=(5, 3 if j % 2 else -9), textcoords="offset points")
    else:
        for j, r in enumerate(sorted(main_set, key=lambda r: -r["specificity"])[:5]):
            ax.annotate(label(r), (r[topk], r[bk]), fontsize=6, color=INK2,
                        xytext=(-70, 10 + 9 * j), textcoords="offset points",
                        arrowprops=dict(arrowstyle="-", color=AXIS, linewidth=0.5))
        for j, r in enumerate(sorted(main_set, key=lambda r: r["specificity"])[:3]):
            ax.annotate(f"{label(r)} [{r['top_body_category']}]", (r[topk], r[bk]), fontsize=6,
                        color=ORANGE, xytext=(-150, -12 - 9 * j), textcoords="offset points",
                        arrowprops=dict(arrowstyle="-", color=AXIS, linewidth=0.5))
    if ratio:
        v = [x for r in summ for x in (r[bk], r[topk]) if x > 0]
        lo_, hi_ = min(v) / 1.6, max(v) * 1.6
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.plot([lo_, hi_], [lo_, hi_], color=GRID, linewidth=1, linestyle=":", zorder=1)
        ax.text(hi_ / 1.15, hi_ / 1.6, "brain = body", fontsize=6, color=INK2, ha="right")
        ax.set_xlim(lo_, hi_)
        ax.set_ylim(lo_, hi_)
    else:
        ax.set_xlim(0, vmax_size)
        ax.set_ylim(0, vmax_size)
    ax.set_xlabel(f"strongest non-neural tissue category ({'size, log scale' if ratio else 'mean quantile'})")
    ax.set_ylabel(f"brain / CNS ({'size, log scale' if ratio else 'mean quantile'})")
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    ax.set_title("Impactful in brain? vs impactful elsewhere?", loc="left", fontsize=9, color=INK)
    fig.text(0.01, 0.005, ref_note, fontsize=5.5, color=INK2)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    fig.savefig(out / f"{stem}_scatter.png", dpi=200)
    plt.close(fig)

    # 4. close-ups: most brain-high and most brain-low (track-level values)
    n = min(args.n_closeups, max(1, len(main_set) // 2))
    picks = main_set[:n] + main_set[-n:][::-1]
    fig, axes = plt.subplots(2, n, figsize=(2.9 * n, 6.2), sharey=not ratio, squeeze=False)   # ratio: own scale
    rng = np.random.default_rng(0)
    vals = signed if ratio else size
    allv = [vals[(r["set"], r["rsid"])][t] for r in picks for t in use]
    lim = max(1.0, float(np.nanpercentile(np.abs(allv), 99))) * 1.05 if ratio else 1.02
    for ax, r in zip(axes.flat, picks):
        k = (r["set"], r["rsid"])
        for i, c in enumerate(categories):
            ys = [vals[k][t] for t in use if track_cat[t] == c]
            xs = i + rng.uniform(-0.22, 0.22, len(ys))
            ax.scatter(xs, ys, s=7, color=BLUE if c == "CNS" else (ORANGE if c in NEURAL_CATS else GREY),
                       alpha=0.85, linewidth=0, zorder=3)
            ax.plot([i - 0.3, i + 0.3], [np.nanmean(ys)] * 2, color=INK, linewidth=1.4, zorder=4)
        if ratio:
            ax.axhline(0, color=AXIS, linewidth=0.8)
        ax.set_xticks(range(len(categories)))
        ax.set_xticklabels([c.replace("_", " ") for c in categories], rotation=90, fontsize=5.5)
        if ratio:
            own = [vals[k][t] for t in use if not math.isnan(vals[k][t])]
            pl = max(1.0, max(abs(v) for v in own)) * 1.08
            ax.set_ylim(-pl, pl)
        else:
            ax.set_ylim(-0.02, lim)
        ax.grid(True, axis="y", color=GRID, linewidth=0.5)
        ax.set_title(f"{label(r)}\nbrain {r[bk]:.2f} · body {r[bodyk]:.2f}", loc="left", fontsize=7.5, color=INK)
    yl = "signed ratio per track" if ratio else "within-track quantile"
    axes[0, 0].set_ylabel(f"{yl}\n(most brain-specific)")
    axes[1, 0].set_ylabel(f"{yl}\n(least brain-specific)")
    fig.text(0.01, 0.004, ref_note + ".  Blue = CNS, orange = other neural, grey = body; bar = category mean.",
             fontsize=6.0, color=INK2)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    fig.savefig(out / f"{stem}_closeups.png", dpi=200)
    plt.close(fig)

    # 5. ratio: denominators by category (control-trait bias check)
    if den:
        denom, low, all_use, cats = den
        mod = args.modality.lower()
        by = defaultdict(list)
        for t in all_use:
            by[cats[(args.modality, t)]["tissue_category"]].append(t)
        order = sorted(by, key=lambda c: -np.median([denom[t] for t in by[c]]))
        fig, ax = plt.subplots(figsize=(6.4, 0.24 * len(order) + 1.4))
        med_all = float(np.median([denom[t] for t in all_use]))
        for i, c in enumerate(order):
            v = [denom[t] for t in by[c]]
            ys = i + rng.uniform(-0.2, 0.2, len(v))
            ax.scatter(v, ys, s=8, linewidth=0, alpha=0.85, zorder=3,
                       color=[ORANGE if t in low else (BLUE if c == "CNS" else GREY) for t in by[c]])
            ax.plot([np.median(v)] * 2, [i - 0.3, i + 0.3], color=INK, linewidth=1.4, zorder=4)
        ax.axvline(med_all, color=AXIS, linewidth=0.8, linestyle="--")
        ax.set_yticks(range(len(order)))
        ax.set_yticklabels([f"{c.replace('_', ' ')} ({len(by[c])})" for c in order])
        ax.invert_yaxis()
        ax.set_xscale("log")
        ax.set_xlabel(f"denominator: mean |score| of the {args.reference} SNPs in the track (log scale)")
        ax.grid(True, axis="x", color=GRID, linewidth=0.5)
        ax.set_title(f"{args.modality}: denominator by tissue category\n"
                     f"bar = median · dashed = median of all tracks · orange = bottom {args.low_pct:g}% (flagged)",
                     loc="left", fontsize=8, color=INK)
        fig.tight_layout()
        fig.savefig(out / f"{mod}_ratio_denominators.png", dpi=200)
        plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
