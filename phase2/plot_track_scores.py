#!/usr/bin/env python3
"""
plot_track_scores.py
--------------------
Figures for looking at ISM variant-effect scores BEFORE choosing a summary
statistic (Choo, 10/12: "plot some of the distributions before we start
summarizing"). Reads the tidy table written by phase2/score_tracks.py.

Figures (DNase by default; --modality to change), written to results/figures/:
  1_per_snp_strips.png     one panel per SNP: every track's |score|, split into
                           brain / neural_model / excluded / background, with the
                           mean of each set marked. Pilots first (ordered by
                           brain-to-background ratio), then controls.
  2_means_variances.png    per SNP: brain vs background mean |score| (left) and
                           variance of |score| (right); pilots vs controls.
  3_score_ranges.png       every SNP's full spread of |score| across all real
                           tracks (range and tails), brain tracks overlaid.
  4_top_tracks_<rsid>.png  the tracks with the largest |score| for chosen SNPs,
                           coloured by set: which tissues carry the effect.

|score| is plotted on a log scale shared across panels, so SNPs with small and
large effects can be compared. Scores below 1e-4 are drawn at 1e-4.

Usage (from the repo root)
    python phase2/plot_track_scores.py
    python phase2/plot_track_scores.py --top rs1198588 rs2710323 rs72730413
"""

import argparse
import csv
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

try:
    import matplotlib
except ImportError:
    sys.exit(f"matplotlib is not installed for this Python ({sys.executable}).\n"
             "Activate the environment that has it first, e.g. on the Mac:\n"
             "    source ~/alphagenome-env/bin/activate\n"
             "or install it:  pip install matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SETS = ["brain", "neural_model", "excluded", "background"]
SET_LABEL = {"brain": "brain", "neural_model": "stem-cell\nneural", "excluded": "excluded",
             "background": "background"}
# Validated categorical palette (light): blue, orange, aqua; background in neutral grey
COLOR = {"brain": "#2a78d6", "neural_model": "#eb6834", "excluded": "#1baf7a", "background": "#898781"}
INK, INK2, MUTED, GRID, AXIS, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
SET_COLOR = {"pilot": "#2a78d6", "control": "#898781"}
FLOOR = 1e-4

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 9, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
    "axes.titlecolor": INK, "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.8,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.spines.top": False, "axes.spines.right": False,
})


def load(path, modality):
    rows = [r for r in csv.DictReader(open(path, newline="", encoding="utf-8"), delimiter="\t")
            if r["modality"] == modality]
    if not rows:
        raise SystemExit(f"No {modality} rows in {path}")
    return rows


def grid(ax, axis="y"):
    ax.grid(True, axis=axis, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def jitter(i, n, width=0.32):
    # deterministic spread so the figure is identical every run
    return 0 if n == 1 else (((i * 0.61803398875) % 1) - 0.5) * 2 * width


def per_snp(rows):
    by = defaultdict(lambda: defaultdict(list))
    info = {}
    for r in rows:
        by[r["rsid"]][r["analysis_set"]].append(max(float(r["abs_score"]), FLOOR))
        info[r["rsid"]] = r["set"]
    return by, info


def mean(v):
    return statistics.fmean(v) if v else float("nan")


def ordered_snps(by, info):
    ratio = {s: mean(by[s]["brain"]) / mean(by[s]["background"]) for s in by}
    pil = sorted([s for s in by if info[s] == "pilot"], key=lambda s: -ratio[s])
    con = sorted([s for s in by if info[s] != "pilot"], key=lambda s: -ratio[s])
    return pil + con, ratio


def fig_strips(by, info, modality, out):
    snps, ratio = ordered_snps(by, info)
    ncol = 5
    nrow = math.ceil(len(snps) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.6 * ncol, 2.5 * nrow), sharey=True)
    axes = axes.ravel()
    allv = [v for s in snps for vs in by[s].values() for v in vs]
    lo, hi = FLOOR * 0.8, max(allv) * 1.5
    for ax, s in zip(axes, snps):
        for k, st in enumerate(SETS):
            vals = by[s][st]
            xs = [k + jitter(i, len(vals)) for i in range(len(vals))]
            ax.scatter(xs, vals, s=9, color=COLOR[st], alpha=0.75 if st != "background" else 0.45,
                       linewidths=0, zorder=2)
            if vals:
                m = mean(vals)
                ax.plot([k - 0.38, k + 0.38], [m, m], color=INK, linewidth=1.6, zorder=3)
        ax.set_yscale("log")
        ax.set_ylim(lo, hi)
        ax.set_xticks(range(len(SETS)))
        ax.set_xticklabels([SET_LABEL[x] for x in SETS], fontsize=7)
        for t, st in zip(ax.get_xticklabels(), SETS):
            t.set_color(INK2)
        grid(ax)
        tag = "pilot" if info[s] == "pilot" else "control"
        ax.set_title(f"{s}  ({tag})\nbrain ÷ background = {ratio[s]:.2f}", fontsize=8, loc="left",
                     color=INK if tag == "pilot" else INK2)
    for ax in axes[len(snps):]:
        ax.axis("off")
    fig.supylabel(f"|variant effect score|  ({modality}, log scale, shared across panels)", color=INK2, fontsize=9)
    fig.suptitle(f"{modality}: every track's |score| for each SNP, by track set  "
                 f"(black bar = mean of the set; n = "
                 + ", ".join(f"{len(by[snps[0]][st])} {st.replace('_', ' ')}" for st in SETS) + ")",
                 x=0.01, ha="left", fontsize=10, color=INK)
    fig.tight_layout(rect=(0.015, 0, 1, 0.97))
    fig.savefig(out, dpi=200)
    plt.close(fig)


def fig_means_variances(by, info, modality, out, label_snps):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6))
    for ax, (stat, name) in zip(axes, [(mean, "mean |score|"),
                                       (lambda v: statistics.pvariance(v) if len(v) > 1 else float("nan"),
                                        "variance of |score|")]):
        xs, ys = [], []
        for s in by:
            x, y = stat(by[s]["background"]), stat(by[s]["brain"])
            xs.append(x)
            ys.append(y)
            tag = "pilot" if info[s] == "pilot" else "control"
            ax.scatter([x], [y], s=46, color=SET_COLOR[tag], edgecolors=SURFACE, linewidths=1.5, zorder=3)
            if s in label_snps:
                left = y > x and s != "rs1198588"  # points just above the diagonal: label on the left
                ax.annotate(s, (x, y), xytext=(-6, 4) if left else (5, 4), textcoords="offset points",
                            fontsize=8, color=INK, ha="right" if left else "left")
        lo = min(min(xs), min(ys)) / 1.6
        hi = max(max(xs), max(ys)) * 1.6
        ax.plot([lo, hi], [lo, hi], color=AXIS, linewidth=1, zorder=1)
        ax.text(hi / 1.15, hi / 1.6, "brain = background", ha="right", fontsize=7, color=MUTED, rotation=0)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect("equal")
        ax.set_xlabel(f"background tracks: {name}")
        ax.set_ylabel(f"brain tracks: {name}")
        grid(ax, "both")
        ax.set_title(f"{modality} {name}: brain vs background", loc="left", fontsize=10)
    h = [plt.Line2D([], [], marker="o", linestyle="", color=SET_COLOR[t], markersize=7, label=t + "s")
         for t in ("pilot", "control")]
    axes[1].legend(handles=h, frameon=False, loc="lower right", labelcolor=INK2)
    fig.text(0.01, 0.01, "Above the grey line: brain tracks respond more than background tracks. "
             "Each point is one SNP.", fontsize=8, color=INK2)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(out, dpi=200)
    plt.close(fig)


def fig_ranges(by, info, modality, out):
    snps, ratio = ordered_snps(by, info)
    fig, ax = plt.subplots(figsize=(9, 0.32 * len(snps) + 1.4))
    for k, s in enumerate(reversed(snps)):
        allv = sorted(v for st in SETS for v in by[s][st])
        q = lambda p: allv[min(len(allv) - 1, int(p * (len(allv) - 1)))]
        ax.plot([allv[0], allv[-1]], [k, k], color=GRID, linewidth=1, zorder=1)
        ax.plot([q(0.25), q(0.75)], [k, k], color=AXIS, linewidth=5, solid_capstyle="round", zorder=2)
        ax.scatter([q(0.5)], [k], s=22, color=INK2, zorder=3)
        ax.scatter([q(0.95)], [k], s=26, marker="|", color=INK2, zorder=3)
        b = by[s]["brain"]
        ax.scatter(b, [k] * len(b), s=11, color=COLOR["brain"], alpha=0.8, linewidths=0, zorder=4)
    ax.set_yticks(range(len(snps)))
    ax.set_yticklabels([f"{s} ({'pilot' if info[s] == 'pilot' else 'control'})" for s in reversed(snps)], fontsize=8)
    for t, s in zip(ax.get_yticklabels(), reversed(snps)):
        t.set_color(INK if info[s] == "pilot" else MUTED)
    ax.set_xscale("log")
    ax.set_xlabel(f"|variant effect score|  ({modality}, all {sum(len(v) for v in by[snps[0]].values())} real tracks, log scale)")
    grid(ax, "x")
    ax.set_title(f"{modality}: spread of |score| across all tracks for each SNP", loc="left", fontsize=10)
    fig.text(0.01, 0.01, "Thin line: full range.  Thick bar: middle 50% of tracks.  Grey dot: median.  "
             "Tick: 95th percentile.  Blue dots: the 16 brain tracks.", fontsize=8, color=INK2)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(out, dpi=200)
    plt.close(fig)


def fig_top_tracks(rows, rsid, modality, out, n=30):
    rs = sorted((r for r in rows if r["rsid"] == rsid), key=lambda r: -float(r["abs_score"]))
    if not rs:
        print(f"  (skipped {rsid}: not in table)")
        return
    top = rs[:n]
    allv = [float(r["abs_score"]) for r in rs]
    p95 = sorted(allv)[int(0.95 * (len(allv) - 1))]
    fig, ax = plt.subplots(figsize=(8.5, 0.24 * n + 1.5))
    for k, r in enumerate(reversed(top)):
        v = float(r["abs_score"])
        ax.barh(k, v, height=0.72, color=COLOR[r["analysis_set"]], zorder=2)
        label = r["biosample_name"].replace("\\'", "'")
        extra = r["life_stage_group"] if r["analysis_set"] == "brain" else r["group"].replace("_", " ")
        ax.text(-0.01 * max(allv), k, f"{label}  [{extra}]", ha="right", va="center", fontsize=7, color=INK2)
    ax.axvline(p95, color=INK2, linewidth=1, linestyle=(0, (3, 3)), zorder=1)
    ax.text(p95, n - 0.2, " 95th percentile of all tracks", fontsize=7, color=INK2, va="bottom")
    ax.set_yticks([])
    ax.set_xlim(0, max(allv) * 1.08)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel(f"|variant effect score|  ({modality})")
    grid(ax, "x")
    first = rs[0]
    ax.set_title(f"{rsid} ({first['set']}): the {n} {modality} tracks with the largest effect", loc="left", fontsize=10)
    h = [plt.Rectangle((0, 0), 1, 1, color=COLOR[st]) for st in SETS]
    ax.legend(h, [SET_LABEL[st].replace("\n", " ") for st in SETS], frameon=False, loc="lower right",
              fontsize=8, labelcolor=INK2)
    fig.tight_layout()
    fig.subplots_adjust(left=0.42)
    fig.savefig(out, dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Plot ISM score distributions by track set.")
    ap.add_argument("--table", default="results/summaries/track_scores_long.tsv")
    ap.add_argument("--modality", default="DNASE")
    ap.add_argument("--outdir", default="results/figures")
    ap.add_argument("--top", nargs="*", default=["rs1198588", "rs2710323", "rs2238057", "rs72730413"],
                    help="SNPs to draw a top-tracks figure for")
    args = ap.parse_args()

    rows = load(args.table, args.modality)
    by, info = per_snp(rows)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    m = args.modality.lower()
    fig_strips(by, info, args.modality, out / f"{m}_1_per_snp_strips.png")
    fig_means_variances(by, info, args.modality, out / f"{m}_2_means_variances.png", set(args.top))
    fig_ranges(by, info, args.modality, out / f"{m}_3_score_ranges.png")
    for s in args.top:
        fig_top_tracks(rows, s, args.modality, out / f"{m}_4_top_tracks_{s}.png")
    print("Figures written to", out)
    for p in sorted(out.glob(f"{m}_*.png")):
        print("  ", p)


if __name__ == "__main__":
    main()