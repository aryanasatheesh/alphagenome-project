#!/usr/bin/env python3
"""
scoring_sensitivity.py
----------------------
Do our conclusions depend on HOW the variant effect score is computed?

Our score is AlphaGenome's recommended DNase scorer (CenterMaskScorer, width 501,
DIFF_LOG2_SUM, 1 bp resolution):

    score = log2(1 + sum ALT) - log2(1 + sum REF)    over the 501 bp centred on the variant

For a few SNPs this script runs the model once (ref + alt) and re-scores every
real DNase track under other reasonable choices:

    window       1, 11, 101, 501, 2001, 10001, 100001 bp, whole 131,072 bp input
    pseudocount  0.01, 0.1, 1, 10, 100          (the "+1" before the log; 501 bp)
    aggregation  log2sum   log2 of the two sums (ours)
                 sumlog2   sum over positions of log2(x+1) differences
                 diff_mean difference of mean signal
                 diff_sum  difference of summed signal
                 l2_log1p  Euclidean distance of log(1+x) profiles (unsigned)
                 log2max   log2 ratio of the peak values
    resolution   501 and 2001 bp windows scored on the 128 bp output

For each setting it reports, per SNP:
    spearman_abs     rank correlation of |score| across the 305 tracks with our baseline
    top20_overlap    how many of the baseline's 20 largest-|score| tracks stay in the top 20
    sign_agree       share of tracks (|baseline| > 0.001) whose sign matches (signed scores only)
    mean_abs_brain, mean_abs_background
    brain_bg_ratio   mean |score| brain / mean |score| background
    brain_pct        brain tracks' mean percentile within the background tracks (0.5 = no enrichment)

Checks (the script stops if one fails), for every real DNase track:
    (1) baseline = run_ism_pilot.compute_variant_effect on the same predictions
        (|diff| < 1e-4; this script adds up in float64, the pipeline in float32)
    (2) baseline = the alphagenome_pytorch library scorer, compute_aggregation(DIFF_LOG2_SUM)
        with create_center_mask(width=501)  (|diff| < 1e-4; skipped with a warning if the
        installed version does not have it)
    (3) run_ism_pilot's score today = the committed score in results/ism/*_ism_scores.tsv
        Stops above 1e-3 (that would mean a different sequence, allele or model). Notes
        anything above 1e-6: the files keep 6 dp, and a newer torch can move the last digits.

Outputs (results/sensitivity/)
    sensitivity_summary.tsv   one row per SNP x setting with the metrics above
    sensitivity_scores.tsv    one row per SNP x setting x track: the score
    sensitivity_windows.png   metrics across window widths, one column per SNP
    profile_<rsid>.png        predicted signal for ref and alt around the variant in the most
                              affected brain and background track, with the windows shaded
    sensitivity.log.txt       inputs with md5, checks, 1 bp vs 128 bp unit comparison
    predictions_<rsid>.npz    ref/alt DNase predictions for the real tracks (+/-10 kb at
                              1 bp; whole input at 128 bp) so the analysis can be redone
                              without the model. Large; keep out of git.

Usage (Mac, from the repo root, with ~/alphagenome-env activated; about as long as
running the ISM pilot on the same number of SNPs)
    PYTORCH_ENABLE_MPS_FALLBACK=1 python phase2/scoring_sensitivity.py
    PYTORCH_ENABLE_MPS_FALLBACK=1 python phase2/scoring_sensitivity.py --rsids rs1198588
"""

import argparse
import csv
import datetime
import hashlib
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
SEQ_LEN = 131_072
CENTER = SEQ_LEN // 2          # 0-based index of the variant in the input
BASELINE = ("window", "501")

# (family, setting label, aggregation, window width in bp (None = whole input), pseudocount, resolution)
SETTINGS = (
    [("window", str(w) if w else "131072", "log2sum", w, 1.0, 1)
     for w in (1, 11, 101, 501, 2001, 10001, 100001, None)]
    + [("pseudocount", f"{pc:g}", "log2sum", 501, pc, 1) for pc in (0.01, 0.1, 1.0, 10.0, 100.0)]
    + [("aggregation", k, k, 501, 1.0, 1)
       for k in ("log2sum", "sumlog2", "diff_mean", "diff_sum", "l2_log1p", "log2max")]
    + [("resolution", f"{w}@128bp", "log2sum", w, 1.0, 128) for w in (501, 2001)]
)
KIND = {(s[0], s[1]): s[2] for s in SETTINGS}
SIGNED = {"log2sum", "sumlog2", "diff_mean", "diff_sum", "log2max"}
METRICS = ["spearman_abs", "top20_overlap", "sign_agree", "mean_abs_brain", "mean_abs_background",
           "brain_bg_ratio", "brain_pct"]


# ============================================================================ pure numpy core
def window_slice(width, n_bins, resolution=1, center_bp=CENTER):
    """Bins covered by a `width` bp window centred on the variant.

    Same arithmetic as alphagenome_pytorch create_center_mask: bp positions
    centre - width//2 up to (not including) centre + width//2 + width%2, and at
    coarser resolution every bin that overlaps that range.
    """
    if width is None:
        return slice(0, n_bins)
    start_bp = center_bp - width // 2
    end_bp = center_bp + width // 2 + width % 2
    return slice(max(0, start_bp // resolution), min(n_bins, (end_bp - 1) // resolution + 1))


def aggregate(ref, alt, kind, pc=1.0):
    """ref, alt: [positions, tracks] inside the window. Returns one score per track."""
    ref = ref.astype(np.float64)
    alt = alt.astype(np.float64)
    if kind == "log2sum":
        return np.log2(alt.sum(0) + pc) - np.log2(ref.sum(0) + pc)
    if kind == "sumlog2":
        return np.log2(alt + 1).sum(0) - np.log2(ref + 1).sum(0)
    if kind == "diff_mean":
        return alt.mean(0) - ref.mean(0)
    if kind == "diff_sum":
        return alt.sum(0) - ref.sum(0)
    if kind == "l2_log1p":
        return np.sqrt(((np.log1p(alt) - np.log1p(ref)) ** 2).sum(0))
    if kind == "log2max":
        return np.log2(alt.max(0) + pc) - np.log2(ref.max(0) + pc)
    raise ValueError(f"unknown aggregation {kind!r}")


def all_settings(p1, p128, real):
    """p1 / p128: {'ref','alt'} -> [bins, channels] at 1 bp / 128 bp. Returns {(family, label): scores}."""
    out = {}
    for fam, label, kind, width, pc, res in SETTINGS:
        p = p1 if res == 1 else p128
        sl = window_slice(width, p["ref"].shape[0], res)
        out[(fam, label)] = aggregate(p["ref"][sl][:, real], p["alt"][sl][:, real], kind, pc)
    return out


def ranks(x):
    """0-based ranks; tied values share their average rank."""
    x = np.asarray(x, dtype=np.float64)
    order = np.argsort(x, kind="mergesort")
    r = np.empty(len(x))
    i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and x[order[j + 1]] == x[order[i]]:
            j += 1
        r[order[i:j + 1]] = (i + j) / 2
        i = j + 1
    return r


def spearman(a, b):
    ra, rb = ranks(a), ranks(b)
    ra -= ra.mean()
    rb -= rb.mean()
    den = math.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / den) if den else float("nan")


def top_k(v, k=20):
    return set(np.argsort(-v, kind="mergesort")[:k].tolist())


def metrics(score, base, is_brain, is_bg, signed):
    a, b = np.abs(score), np.abs(base)
    bg_sorted = np.sort(a[is_bg])
    keep = b > 1e-3
    m_brain, m_bg = float(a[is_brain].mean()), float(a[is_bg].mean())
    return {
        "spearman_abs": spearman(a, b),
        "top20_overlap": len(top_k(a) & top_k(b)),
        "sign_agree": (float(np.mean(np.sign(score[keep]) == np.sign(base[keep])))
                       if signed and keep.any() else float("nan")),
        "mean_abs_brain": m_brain,
        "mean_abs_background": m_bg,
        "brain_bg_ratio": m_brain / m_bg if m_bg > 0 else float("nan"),
        # share of background tracks below each brain track (ties count half), averaged over brain tracks
        "brain_pct": float(np.mean([(np.searchsorted(bg_sorted, v, "left") + np.searchsorted(bg_sorted, v, "right"))
                                    / 2 / len(bg_sorted) for v in a[is_brain]])),
    }


def fmt(v):
    if isinstance(v, int):
        return str(v)
    return "" if v is None or math.isnan(v) else f"{v:.6g}"


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ============================================================================ figures
INK, INK2, MUTED, GRID, AXIS, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
BLUE = "#2a78d6"


def _style(plt):
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": AXIS, "xtick.color": MUTED, "ytick.color": MUTED,
                         "axes.labelcolor": INK2, "figure.facecolor": SURF, "axes.facecolor": SURF,
                         "savefig.facecolor": SURF, "axes.spines.top": False, "axes.spines.right": False})


def plot_windows(summary, rsids, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _style(plt)
    rows = [("brain_bg_ratio", "brain ÷ background\nmean |score|", 1.0),
            ("brain_pct", "brain tracks' mean\npercentile in background", 0.5),
            ("spearman_abs", "track ranking vs 501 bp\n(Spearman of |score|)", None),
            ("top20_overlap", "top-20 tracks shared\nwith 501 bp", None)]
    fig, axes = plt.subplots(len(rows), len(rsids), figsize=(2.6 * len(rsids) + 0.8, 1.9 * len(rows)),
                             sharex=True, squeeze=False)
    for j, rs in enumerate(rsids):
        pts = sorted(((int(s["setting"]), s) for s in summary if s["rsid"] == rs and s["family"] == "window"),
                     key=lambda t: t[0])
        xs = [p[0] for p in pts]
        for i, (key, label, ref) in enumerate(rows):
            ax = axes[i][j]
            ys = [float(p[1][key]) for p in pts]
            ax.plot(xs, ys, color=BLUE, linewidth=2, marker="o", markersize=4.5, zorder=3)
            ax.axvline(501, color=INK2, linewidth=1, linestyle=(0, (3, 3)), zorder=1)
            if ref is not None:
                ax.axhline(ref, color=AXIS, linewidth=1.2, zorder=1)
            if key == "top20_overlap":
                ax.set_ylim(-1, 21)
            elif key == "spearman_abs":
                ax.set_ylim(min(0, min(ys)) - 0.05, 1.05)
            ax.set_xscale("log")
            ax.grid(True, color=GRID, linewidth=0.6)
            ax.set_axisbelow(True)
            if j == 0:
                ax.set_ylabel(label)
            if i == 0:
                ax.set_title(rs, loc="left", fontsize=10, color=INK)
            if i == len(rows) - 1:
                ax.set_xlabel("window width (bp)")
    fig.text(0.01, 0.006, "Dashed: the recommended 501 bp window we use.  Grey line: no brain enrichment "
             "(ratio 1, percentile 0.5).  Rightmost point: the whole 131,072 bp input.",
             fontsize=8, color=INK2)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out, dpi=200)
    plt.close(fig)


def plot_profile(p1, real, names, sets, base, rs, out, half=1500):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _style(plt)
    a = np.abs(base)
    picks = []
    for s in ("brain", "background"):
        idx = [i for i in range(len(real)) if sets[i] == s]
        if idx:
            picks.append(max(idx, key=lambda i: a[i]))
    fig, axes = plt.subplots(len(picks), 1, figsize=(9, 2.5 * len(picks)), squeeze=False)
    x = np.arange(-half, half + 1)
    sl = slice(CENTER - half, CENTER + half + 1)
    for ax, i in zip(axes[:, 0], picks):
        ch = real[i]
        ax.axvspan(-1000, 1000, color="#f0efec", zorder=0, linewidth=0)
        ax.axvspan(-250, 250, color=GRID, zorder=0, linewidth=0)
        ax.plot(x, p1["ref"][sl, ch], color=MUTED, linewidth=1.3, label="reference allele", zorder=2)
        ax.plot(x, p1["alt"][sl, ch], color=BLUE, linewidth=1.3, label="alternate allele", zorder=3)
        ax.axvline(0, color=INK2, linewidth=0.8, linestyle=(0, (2, 2)), zorder=4)
        ax.set_title(f"{names[i]}  ({sets[i]})    score at 501 bp = {base[i]:+.3f}",
                     loc="left", fontsize=9, color=INK)
        ax.set_ylabel("predicted DNase")
        ax.grid(True, axis="y", color=GRID, linewidth=0.6)
        ax.legend(frameon=False, fontsize=8, loc="upper right", labelcolor=INK2)
    axes[-1, 0].set_xlabel("position relative to the variant (bp)")
    fig.suptitle(f"{rs}: the signal the score adds up.  Dark band = 501 bp window, light band = 2001 bp.",
                 x=0.01, ha="left", fontsize=10, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out, dpi=200)
    plt.close(fig)


# ============================================================================ main (needs the model)
def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description="Sensitivity of DNase variant scores to scoring choices.")
    ap.add_argument("--rsids", nargs="+", default=["rs1198588", "rs2710323", "rs2238057", "rs72730413"])
    ap.add_argument("--snp-files", nargs="+", default=["phase1/data/scz_pgc3_pilot_10snps_hg38.tsv",
                                                       "phase1/data/scz_pgc3_control_10snps_hg38.tsv"])
    ap.add_argument("--score-files", nargs="+", default=["results/ism/pilot_ism_scores.tsv",
                                                         "results/ism/control_ism_scores.tsv"])
    ap.add_argument("--classes", default="phase2/metadata/track_classes.tsv")
    ap.add_argument("--outdir", default="results/sensitivity")
    ap.add_argument("--no-save-predictions", action="store_true")
    args = ap.parse_args()

    for p in [args.classes, *args.snp_files, *args.score_files]:
        if not Path(p).exists():
            sys.exit(f"Not found: {p}  (run from the repo root)")

    # ---- DNase track labels (real tracks only, in channel order)
    cls = sorted((r for r in read_tsv(args.classes) if r["output_type"] == "DNASE"),
                 key=lambda r: int(r["track_index"]))
    cls = [r for r in cls if r["is_padding"] != "True"]
    real = np.array([int(r["track_index"]) for r in cls])
    names = [r["biosample_name"] for r in cls]
    sets = [r["analysis_set"] for r in cls]
    is_brain = np.array([s == "brain" for s in sets])
    is_bg = np.array([s == "background" for s in sets])
    print(f"DNase tracks: {len(real)} real  ({is_brain.sum()} brain, {is_bg.sum()} background)")

    snps = {r["rsid"]: r for f in args.snp_files for r in read_tsv(f)}
    missing = [r for r in args.rsids if r not in snps]
    if missing:
        sys.exit(f"Not in {args.snp_files}: {missing}")
    committed = {(r["rsid"], int(r["track_index"])): float(r["variant_effect_score"])
                 for f in args.score_files for r in read_tsv(f)
                 if r["modality"].lower() == "dnase" and r["rsid"] in args.rsids}

    # ---- the code we already score with, and the model
    sys.argv = sys.argv[:1]                       # run_ism_pilot reads sys.argv when imported
    sys.path.insert(0, str(REPO))
    import torch
    import run_ism_pilot as ism
    from alphagenome_pytorch import AlphaGenome
    try:
        from alphagenome_pytorch.variant_scoring.aggregations import compute_aggregation, create_center_mask
        from alphagenome_pytorch.variant_scoring.types import AggregationType
        have_lib = True
    except ImportError as e:
        print(f"WARNING: library scorer not importable ({e}); check (3) skipped")
        have_lib = False

    device = "mps" if torch.backends.mps.is_available() else ("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}\nLoading model {ism.WEIGHTS_PATH} ...")
    model = AlphaGenome.from_pretrained(ism.WEIGHTS_PATH, device=device)
    model.eval()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    log = [f"date    : {datetime.datetime.now().isoformat(timespec='seconds')}",
           f"labels  : {args.classes}  md5 {md5(args.classes)}"]
    log += [f"input   : {p}  md5 {md5(p)}" for p in [*args.snp_files, *args.score_files]]
    log += [f"weights : {ism.WEIGHTS_PATH}", f"device  : {device}", ""]

    score_rows, summary = [], []
    for rs in args.rsids:
        snp = snps[rs]
        chrom = snp["chr"][3:] if snp["chr"].lower().startswith("chr") else snp["chr"]
        pos = int(snp["pos"])
        print(f"\n--- {rs}  chr{chrom}:{pos} ---")
        seq = ism.fetch_hg38_sequence(chrom, pos - CENTER, pos + CENTER - 1)   # same window as the ISM run
        if len(seq) != SEQ_LEN:
            sys.exit(f"{rs}: expected {SEQ_LEN} bp from Ensembl, got {len(seq)}")
        a = ism.assign_alleles(seq[CENTER], snp["a1"], snp["a2"], ism.parse_beta(snp["beta"]))
        if a["skip_reason"]:
            sys.exit(f"{rs}: cannot score ({a['skip_reason']})")
        alt_seq = seq[:CENTER] + a["hg38_alt"] + seq[CENTER + 1:]
        print(f"  hg38 ref {a['hg38_ref']} ({a['hg38_ref_is']})  alt {a['hg38_alt']}   running model (ref, alt) ...")
        with torch.inference_mode():
            ref_out = model.predict(ism.one_hot_encode(seq).to(device), organism_index=0)
            alt_out = model.predict(ism.one_hot_encode(alt_seq).to(device), organism_index=0)
        p1 = {k: o["dnase"][1][0].float().cpu().numpy() for k, o in (("ref", ref_out), ("alt", alt_out))}
        p128 = {k: o["dnase"][128][0].float().cpu().numpy() for k, o in (("ref", ref_out), ("alt", alt_out))}

        res = all_settings(p1, p128, real)
        base = res[BASELINE]

        # ---- checks on the baseline
        problems = []
        pipeline = np.asarray(ism.compute_variant_effect(ref_out["dnase"][1], alt_out["dnase"][1], CENTER))[real]
        d_ism = float(np.max(np.abs(base - pipeline)))
        if d_ism > 1e-4:
            problems.append(f"differs from run_ism_pilot.compute_variant_effect by {d_ism:.2e}")
        com = np.array([committed.get((rs, int(t)), np.nan) for t in real])
        n_com = int(np.sum(~np.isnan(com)))
        d_com = float(np.nanmax(np.abs(pipeline - com))) if n_com else float("nan")
        if n_com and n_com != len(real):
            problems.append(f"committed score file has {n_com} of {len(real)} tracks")
        if n_com and d_com > 1e-3:
            problems.append(f"differs from the committed scores by {d_com:.2e} (different sequence, "
                            f"alleles or weights?)")
        elif n_com and d_com > 1e-6:
            print(f"  NOTE: differs from the committed scores by {d_com:.1e}: more than rounding to 6 dp, "
                  f"small enough to be a torch/MPS version difference")
        d_lib = float("nan")
        if have_lib:
            mask = create_center_mask(variant_position=pos, interval_start=pos - CENTER - 1, width=501,
                                      seq_length=SEQ_LEN, resolution=1).unsqueeze(0)
            lib = compute_aggregation(ref_out["dnase"][1].float().cpu(), alt_out["dnase"][1].float().cpu(),
                                      AggregationType.DIFF_LOG2_SUM, mask=mask)
            d_lib = float(np.max(np.abs(base - np.asarray(lib[0].numpy())[real])))
            if d_lib > 1e-4:
                problems.append(f"differs from the library DIFF_LOG2_SUM scorer by {d_lib:.2e}")
        msg = (f"checks  : {rs}  max |diff| baseline vs run_ism_pilot {d_ism:.1e} | baseline vs library "
               f"{'skipped' if math.isnan(d_lib) else f'{d_lib:.1e}'} | run_ism_pilot vs committed file "
               f"{'not found' if not n_com else f'{d_com:.1e} ({n_com} tracks)'}")
        print("  " + msg[10:])
        log.append(msg)
        if problems:
            log.append(f"STOPPED : {rs}: " + "; ".join(problems))
            (out / "sensitivity.log.txt").write_text("\n".join(log) + "\n")
            sys.exit(f"STOPPED: {rs} baseline " + "; ".join(problems))

        # units of the 128 bp output: compare with the 1 bp output over the same 4 bins (512 bp)
        sl128 = window_slice(501, p128["ref"].shape[0], 128)
        s128 = p128["ref"][sl128][:, real].sum(0)
        s1 = p1["ref"][sl128.start * 128:sl128.stop * 128][:, real].sum(0)
        ok = s1 > 1
        if ok.any():
            log.append(f"units   : {rs}  sum of 128 bp bins / sum of 1 bp values over the same "
                       f"{(sl128.stop - sl128.start) * 128} bp: median {np.median(s128[ok] / s1[ok]):.3g} "
                       f"across {int(ok.sum())} tracks (~1 = bins hold sums, ~0.0078 = bins hold means)")

        for (fam, label), sc in res.items():
            m = metrics(sc, base, is_brain, is_bg, KIND[(fam, label)] in SIGNED)
            summary.append({"rsid": rs, "family": fam, "setting": label, "aggregation": KIND[(fam, label)],
                            **{k: fmt(m[k]) for k in METRICS}})
            score_rows.extend((rs, fam, label, int(t), f"{v + 0.0:.6f}") for t, v in zip(real, sc))
        plot_profile(p1, real, names, sets, base, rs, out / f"profile_{rs}.png")
        if not args.no_save_predictions:
            w = slice(CENTER - 10_000, CENTER + 10_001)
            np.savez_compressed(out / f"predictions_{rs}.npz",
                                ref_1bp=p1["ref"][w][:, real], alt_1bp=p1["alt"][w][:, real],
                                ref_128bp=p128["ref"][:, real], alt_128bp=p128["alt"][:, real],
                                track_index=real, first_1bp_offset=-10_000)
        del ref_out, alt_out
        if device == "mps":
            torch.mps.empty_cache()

    with open(out / "sensitivity_scores.tsv", "w", encoding="utf-8") as f:
        f.write("rsid\tfamily\tsetting\ttrack_index\tscore\n")
        f.writelines("\t".join(map(str, r)) + "\n" for r in score_rows)
    with open(out / "sensitivity_summary.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(summary)
    plot_windows(summary, args.rsids, out / "sensitivity_windows.png")
    log.append(f"\nwritten : {len(summary)} summary rows, {len(score_rows)} score rows")
    (out / "sensitivity.log.txt").write_text("\n".join(log) + "\n")
    print_summary(summary, args.rsids)
    print(f"\nWritten to {out}/")
    return 0


def print_summary(summary, rsids):
    head = f"{'family':<12}{'setting':<10}{'spearman':>9}{'top20':>6}{'sign':>6}{'brain/bg':>9}{'brain pct':>10}"
    for rs in rsids:
        print(f"\n{rs}\n{head}")
        for s in summary:
            if s["rsid"] != rs:
                continue
            sign = f"{float(s['sign_agree']):.2f}" if s["sign_agree"] else "-"
            print(f"{s['family']:<12}{s['setting']:<10}{float(s['spearman_abs']):>9.3f}{s['top20_overlap']:>6}"
                  f"{sign:>6}{float(s['brain_bg_ratio']):>9.2f}{float(s['brain_pct']):>10.2f}")


if __name__ == "__main__":
    sys.exit(main())