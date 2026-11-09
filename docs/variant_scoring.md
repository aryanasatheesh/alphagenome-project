# How the variant effect score is computed

*Written for Session D3 (Oct 2026). Answers Choo's 10/12 question about how the
score's summation works, and explains what `phase2/scoring_sensitivity.py` tests.*

## The score in one line

For every track, the score is

    score = log2(1 + Σ ALT) − log2(1 + Σ REF)  =  log2[(Σ ALT + 1) / (Σ REF + 1)]

where Σ is the sum of the model's predicted signal over the **501 bp centred on the
variant** (250 bp on each side plus the variant itself), at 1 bp resolution.

## Step by step

1. **Two predictions.** The model runs once on the 131,072 bp reference sequence and
   once on the same sequence with the single variant base swapped. Both return, for
   each DNase track, one predicted value per base pair, so 131,072 values per track.
2. **Window.** Keep only positions 65,286 to 65,786 (0-based, end inclusive). Position
   65,536 is the variant. That leaves 501 values per track per allele.
3. **Sum.** Add the 501 values. This number is the total predicted accessibility of
   the region around the variant, like the area under the curve in
   `results/sensitivity/profile_<rsid>.png`.
4. **Pseudocount.** Add 1 to each sum. This avoids log(0) when a track predicts
   nothing, and it shrinks fold changes in regions with very little signal.
5. **Log2 ratio.** Take log2(ALT sum / REF sum). The result is a fold change:
   +1 means the alternate allele doubles predicted accessibility, −1 halves it, and
   0 means no change.

Worked examples:

| Σ REF | Σ ALT | score | Meaning |
|---|---|---|---|
| 200 | 300 | log2(301/201) = **+0.58** | ALT about 1.5× more open |
| 300 | 200 | **−0.58** | the same effect in the other direction (the score is antisymmetric) |
| 0.5 | 1.0 | log2(2.0/1.5) = **+0.42** | without the +1 it would be +1.00: a closed region gets damped |
| 0 | 0 | **0** | no signal, no score |

The sign refers to hg38 ALT vs hg38 REF. The risk-oriented score
(`docs/allele_conventions.md`) turns it into risk vs non-risk allele.

## This is AlphaGenome's own recommended scorer

Checked against the source code (Oct 2026):

- `google-deepmind/alphagenome`, `src/alphagenome/models/variant_scorers.py`:
  `RECOMMENDED_VARIANT_SCORERS['DNASE'] = CenterMaskScorer(requested_output=DNASE,
  width=501, aggregation_type=DIFF_LOG2_SUM)`. ATAC uses the same; histone ChIP uses
  width 2001. Allowed widths are 501, 2001, 10001, 100001, 200001 or the whole input.
- `alphagenome-pytorch` 0.3.1, `variant_scoring/aggregations.py`: `DIFF_LOG2_SUM` is
  `log2(1 + alt_sum) − log2(1 + ref_sum)`, and `create_center_mask` uses the same window
  arithmetic as `run_ism_pilot.compute_variant_effect`.
- `scoring_sensitivity.py` recomputes the score three ways (our pipeline, the library
  function, our own float64 re-implementation) and stops if they disagree.

**Correction to the July notebook:** the July notes call this aggregation
`DIFF_SUM_LOG2`. That name belongs to a different scorer, Σ log2(x+1) over positions,
which is a sum of per-base log changes. Ours is `DIFF_LOG2_SUM`: the log of the sums.

## Why these choices, and what could change the answer

| Choice | Why this value | What could go wrong | Tested by |
|---|---|---|---|
| 501 bp window | About the width of one DNase peak (a regulatory element). It is the library default for DNase/ATAC. | A narrower window only sees the variant's own base. A much wider one adds distant peaks that the variant doesn't touch, which dilutes the effect. | `family = window`: 1 bp to the whole input |
| Sum, then log | Measures change in the element's total accessibility | A per-base log sum (`sumlog2`) counts many small shifts. Max (`log2max`) only sees the peak summit. | `family = aggregation` |
| +1 pseudocount | Library default | If predicted values are mostly below 1, the +1 compresses differences between tracks. That would make strongly accessible tissues look more affected. | `family = pseudocount`: 0.01 to 100 |
| 1 bp resolution | Most precise output | Histone ChIP exists only at 128 bp, so later work will use bins. | `family = resolution` |

**How to read `sensitivity_summary.tsv`.** For each SNP and setting:

- `spearman_abs` and `top20_overlap` ask whether the same tracks are the most affected.
- `brain_bg_ratio` and `brain_pct` ask whether the brain-vs-background conclusion holds.
- `sign_agree` asks whether the direction holds.

A conclusion is robust if it holds across the nearby settings: windows 101 to 2001,
pseudocounts 0.1 to 10, and the three sum-based aggregations. It doesn't have to hold
at extremes like the 1 bp window or the whole input.

---

## D3 results: how much does the scoring choice matter? (run Nov 2026)

`phase2/scoring_sensitivity.py`, DNase, 305 real tracks (16 brain, 255 background).
Full table: `results/sensitivity/sensitivity_summary.tsv`. Figure: `sensitivity_windows.png`.

### Checks passed

- Re-implementation vs `run_ism_pilot.compute_variant_effect`: max |diff| 2.6–3.7 × 10⁻⁷.
  The only difference is float64 vs float32 arithmetic.
- Today's model run vs the committed score files: max |diff| 4.9–5.0 × 10⁻⁷. That is the
  6-decimal rounding of the files. **The October scores are reproducible from the
  weights today.**
- Library scorer: *skipped*. The installed `alphagenome-pytorch` has no
  `variant_scoring` module. The equivalence to `DIFF_LOG2_SUM` therefore rests on
  reading the 0.3.1 source code (above).
- `diff_sum` and `diff_mean` agree in sign with our score on 100% of tracks for every
  SNP, as they must: log2((A+1)/(R+1)) has the sign of A − R.

### Brain ÷ background (ratio) and brain tracks' mean percentile within background (pct)

| Setting | rs1198588 | rs2238057 | rs2710323 | control rs72730413 |
|---|---|---|---|---|
| **501 bp (ours)** | **4.62 / 0.72** | **1.73 / 0.73** | **0.27 / 0.23** | **0.57 / 0.40** |
| window 101 | 2.32 / 0.72 | 1.33 / 0.69 | 0.32 / 0.23 | 0.54 / 0.33 |
| window 2001 | 6.32 / 0.77 | 2.05 / 0.74 | 0.27 / 0.25 | 0.58 / 0.45 |
| window 10001 | 4.86 / 0.81 | 1.87 / 0.79 | 0.35 / 0.33 | 0.54 / 0.45 |
| window 1 | 2.63 / 0.65 | 0.85 / 0.48 | 0.24 / 0.20 | 0.95 / 0.62 |
| whole input | 2.57 / 0.74 | 1.03 / 0.50 | 0.51 / 0.42 | 0.82 / 0.48 |
| pseudocount 0.01–10 | 4.43–4.91 / 0.71–0.78 | 1.72–1.73 / 0.73 | 0.25–0.27 / 0.22–0.23 | 0.47–0.59 / 0.39–0.42 |
| sum of logs | 3.17 / 0.81 | 1.31 / 0.71 | 0.19 / 0.22 | 0.36 / 0.45 |
| difference of sums | 2.46 / 0.79 | 1.11 / 0.70 | 0.17 / 0.22 | 0.34 / 0.46 |
| L2 of log1p | 2.72 / 0.84 | 0.94 / 0.58 | 0.27 / 0.21 | 0.49 / 0.45 |
| log2 of max | 6.91 / 0.85 | 0.98 / 0.63 | 0.15 / 0.19 | 0.53 / 0.40 |
| 501 bp at 128 bp bins | 5.54 / 0.81 | 1.65 / 0.74 | 0.22 / 0.17 | 0.52 / 0.33 |
| **range, all 21 settings** | **1.29–6.97 / 0.61–0.85** | **0.85–2.05 / 0.46–0.79** | **0.15–0.51 / 0.17–0.42** | **0.34–0.95 / 0.33–0.62** |

### Conclusions

1. **The pseudocount does not matter.** From 0.01 to 10, track rankings barely move
   (Spearman ≥ 0.93 for all four SNPs), signs agree on 100% of tracks, and the
   brain ratio changes by less than 10% for the three pilots (control: 0.47–0.59).
   Only 100 starts to compress scores.
2. **Every SNP's classification survives the reasonable settings.** Here "reasonable"
   means windows 101–10001, pseudocount 0.1–10 and the sum-based aggregations:
   - **rs1198588** is brain-enriched at **all 21 settings**, including the extremes
     (ratio ≥ 1.29, percentile ≥ 0.61).
   - **rs2710323** is brain-depleted (immune-dominated) at all 21 (ratio ≤ 0.51, percentile ≤ 0.42).
   - The **control** is never brain-enriched by ratio (≤ 0.95). By percentile it is
     ≤ 0.48 everywhere except the 1 bp and 11 bp windows.
3. **rs2238057 is the fragile one.** It is enriched across windows 101–10001 and the
   fold-change scores (ratio 1.31–2.05, percentile 0.69–0.79). With L2, peak-only or
   absolute-difference scoring the ratio falls to 0.94–1.11, though its brain tracks
   still rank above the background median (percentile 0.58–0.70). At the extreme
   windows (1, 11 bp; ≥ 100 kb) the enrichment disappears (ratio 0.85–1.06,
   percentile 0.46–0.50). Call it brain-leaning, not brain-specific.
4. **The size of the ratio depends on the window; the direction of the conclusion does
   not.** For rs1198588 the ratio runs from 2.32 (101 bp) to 6.32 (2001 bp). Report
   the 501 bp value with this range, and don't compare ratios across window choices.
5. **Exact track rankings are window-dependent for rs1198588.** Its Spearman with 501 bp
   is 0.70 at both 101 and 2001 bp, versus 0.85–0.99 for the other SNPs. Its effect has
   spatial structure beyond the 501 bp element (see `profile_rs1198588.png`). Report
   **sets** of tracks (e.g. "adult basal ganglia and cortex"), not a single top track.
6. **The percentile is more stable than the ratio.** Across all 21 settings
   rs1198588's percentile stays between 0.61 and 0.85, while its ratio varies more than
   5-fold (1.29–6.97). rs2710323's percentile never rises above 0.42. A mean ratio is pulled by a few large tracks;
   a percentile is not. **This favours a rank/percentile-based specificity statistic.**
7. **Absolute-change scoring down-weights brain.** Difference of sums and difference of
   means shrink brain enrichment (rs1198588 2.46, rs2238057 1.11). Absolute change grows
   with a track's baseline accessibility, so tracks that are highly open at the locus
   dominate. Fold change (ours) asks about relative change in each tissue, which is the
   right question for comparing tissues.
8. **128 bp resolution keeps the set-level conclusions** (rs1198588 5.54 / 6.97;
   rs2238057 1.65 / 2.03), but track rankings agree less (Spearman 0.40–0.93) and
   directions less often (63–92% of tracks). Histone ChIP (128 bp only) is usable for
   tissue-set questions. Track-level direction claims from it need care.

### Decision

**Keep the pre-specified 501 bp DIFF_LOG2_SUM score with pseudocount 1.** It is
AlphaGenome's recommended scorer, and the conclusions hold around it. The 2001 bp window
gives larger brain ratios. Switching to it because of that would mean choosing the
analysis after seeing the result. 2001 bp stays reserved for histone marks, as
AlphaGenome recommends.

### Still open

- How 128 bp bins relate to 1 bp values (sums or means): the `units` line in
  `sensitivity.log.txt`.
- Whether the tracks that change sign between windows are small-effect tracks:
  `phase2/sensitivity_followup.py`.