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