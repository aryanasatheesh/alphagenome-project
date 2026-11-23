# How we turn AlphaGenome output into tissue profiles, and what should change before the big runs

*Nov 18, 2026. Based on the code at commit `5ca2b7c` (`run_ism_pilot.py`, `phase2/tissue_profiles.py`)
and the AlphaGenome documentation.*

## 1. What the raw output is
- **Input:** 131,072 bp of hg38 sequence, one-hot encoded, centred on the variant.
- **Two predictions per variant:**
  - REF: the reference sequence
  - ALT: the same sequence with one base changed

  For each output head (DNase 384 channels, ATAC 256; 305 and 167 real), the model returns a
  **predicted signal profile at 1 bp resolution across the whole window**. This is its
  estimate of the coverage that the experiment behind each track would produce for this
  sequence in that cell type.
- **What a profile is and isn't:**
  - It is a prediction of a measurement, not a measurement.
  - It is specific to one cell type's experiment: the track head carries that cell type's
    regulatory environment.
  - Its scale is the scale of that experiment's normalised coverage, which differs between
    tracks.
- **Comparing ref and alt:** two curves per track, almost identical except near the variant,
  where a change in TF binding or nucleosome positioning raises or lowers the predicted
  accessibility of the surrounding element.

## 2. What we do with it now
1. **Score per track** (`compute_variant_effect`):
   - sum the REF profile and the ALT profile over the **501 bp** centred on the variant
   - score = log2(1 + Σ ALT) − log2(1 + Σ REF)
   - one number per track: the fold change in the total predicted accessibility of the local
     element (+1 doubles it, −1 halves it)
   - this is AlphaGenome's recommended DNase/ATAC scorer (CenterMaskScorer, width 501,
     DIFF_LOG2_SUM), checked against the source in D3
2. **Orient to risk** (`run_ism_pilot.py`):
   - ref/alt are set from the hg38 base, not the GWAS columns
   - the sign is flipped so that + means the GWAS **risk** allele opens chromatin
     (`risk_oriented_score`)
   - palindromic SNPs keep their size; their direction assumes the + strand
3. **Make tracks comparable** (`tissue_profiles.py`): divide each track's score by a
   reference scale for that track. Currently that's the mean |score| of 20 T1D leads.
4. **Summarise by tissue:** per organ-system category (equal weight, cell lines out)
   - size = mean |normalised score|
   - direction = mean signed normalised score
   - concordance = share of tracks agreeing in sign
5. **Assays stay separate.** DNase and ATAC are never pooled. ATAC has no brain tracks, so
   brain questions use DNase only.

## 3. Why the core score is the right one (keep it)
- **It answers the right question for tissue comparisons.** A fold change measures the
  *relative* change in each cell type. A tissue where the element is barely open and one
  where it is wide open can both register "the risk allele halves it". An absolute difference
  would favour tissues with high baseline coverage. D3 showed exactly that: absolute-difference
  scores shrank brain enrichment by up-weighting tracks that are highly accessible at baseline.
- **The window matches the biology.** 501 bp is the size of a regulatory element (a DNase
  hypersensitive site is a few hundred bp). Accessibility effects of a SNP are mostly local.
- **It's the developers' choice.** It is AlphaGenome's own recommended accessibility scorer,
  so our numbers are comparable with other uses of the model.
- **It's robust (D3).** Conclusions held across windows of 101–2001 bp, pseudocounts 0.1–10
  and the sum-based aggregations. Direction is stable wherever the effect is real: the 100
  most affected tracks never flip sign.
- **It's reproducible.** It matches the committed scores to within 5 × 10⁻⁷.

## 4. Where the current approach falls short, and what to do
In order of importance. Items A–C must be settled **before** the GPU runs (173 leads + matched
controls), because changing them later means rerunning the model.

### A. Model weights: we use one cross-validation fold, not the variant-scoring model
- **Now:** `WEIGHTS_PATH = ~/alphagenome-weights/model_fold_0.safetensors`. That model was
  trained with one quarter of the genome held out, for benchmarking track prediction.
- **Recommended:**
  - the PyTorch port recommends `model_all_folds.safetensors` (trained on all data)
  - the paper's variant-scoring model is a single distilled model that reproduces the all-fold
    teachers
  - a lead SNP that falls in fold 0's held-out region was never seen in training by our model
- **Do:**
  1. Confirm which released file is the distilled/all-folds model.
  2. Rescore the 40 SNPs already scored (10 pilots, 10 random, 20 T1D) on the Mac.
  3. Compare per-track scores (Spearman, sign agreement, the pilots' category profiles)
     between fold 0 and the recommended weights.
  4. Switch to the recommended weights for all future runs. If the scores change materially,
     report the October–November results as fold-0 results.

### B. Input length: 131 kb instead of the 1 Mb the model was trained and scored on
- **Now:** 131,072 bp, chosen because it fits on the Mac.
- **Official:** the scoring examples resize the interval to 1 Mb. Shorter inputs remove
  distal context; for local accessibility the effect is probably small, but it is untested.
- **Do:** on a GPU, score a handful of SNPs at both 131 kb and 1 Mb and compare as in A. If they
  agree closely, keep 131 kb (cheaper) and record the check. If not, run the big jobs at 1 Mb.

### C. We throw away the baseline accessibility
- **Now:** Σ REF and Σ ALT are computed, then only their log ratio is saved.
- **Why it matters:**
  - **Closed-chromatin effects:** a fold change can be large where the element is barely
    accessible (Σ 0.5 → 1.0 gives +0.42). That is a change in a region that is closed in that
    cell type, which matters less biologically than the same fold change at an open element.
  - **Activity gating and absolute effects** both need the sums.
  - **Changing the pseudocount later** needs the sums.
- **Do:** save Σ REF and Σ ALT per track in the same run. Zero extra compute: they are already
  computed. Then:
  - report whether the element is accessible in each tissue (e.g. Σ REF above that track's
    median over the matched controls)
  - mark effects in closed elements, instead of silently mixing them in

### D. The per-track scale must be robust, and come from a matched null
- **Now:** mean |score| of 20 T1D leads.
- **What the T1D run showed (Session 17):**
  - split-half reliability ≈ 0.53
  - the T1D and random-SNP scales disagree (Spearman −0.18)
  - one T1D SNP supplies a median 39% of the immune-track scale
  - the scale carries T1D's own tissue biology into every body comparison
  - "1" separates nothing from random SNPs
- **Two separate problems:**
  - **Too few reference variants, and the wrong kind:** fix with matched non-associated SNPs,
    pooled across leads (~865 at full scale).
  - **The mean is not a robust scale for heavy-tailed scores:** a few large effects dominate it.
- **Do:**
  - use the pooled matched controls as the reference
  - make the scale robust: the median |score| per track, or, with hundreds of controls,
    AlphaGenome's own approach of a per-track quantile against the reference
  - keep the ratio-to-mean as a sensitivity analysis

### E. Category summaries: the mean answers one question, and we need two
- The **category mean** asks "does this tissue as a whole respond?"
- For targets, the relevant question is often "is there **any** cell type in this tissue where
  the variant acts?" For rs1198588, a few basal-ganglia tracks carry the brain signal, and
  the mean dilutes it.
- A maximum on its own is biased towards categories with many tracks: more tracks, more
  chances.
- **Do:** report mean and max (or 90th percentile) per category, each compared with the same
  statistic computed on the matched controls. Comparing a statistic with its own null
  distribution removes the track-count bias.

### F. Checks worth adding (cheap)
- **Reverse complement:** score a few SNPs with the reverse-complemented sequence as well. The
  model was trained with reverse-complement augmentation, so the two should agree. If they
  differ noticeably, average them.
- **Effects that move rather than change accessibility:** a sum over 501 bp misses a variant
  that redistributes accessibility without changing the total. That's acceptable for
  accessibility, and it is the developers' choice; list it as a limitation.

## 5. What a score means, and what it doesn't
**A normalised score answers one question:** how much the model predicts the risk allele changes
the accessibility of the local element in the cell type of each track, relative to comparable
non-associated variants.

**It does not answer:**
- whether the variant is causal: lead SNPs tag a locus; fine-mapping comes later
- whether gene expression changes: that needs the RNA/CAGE heads and gene-level scorers
- whether the prediction is right: that needs comparison with measured allele-specific
  accessibility, planned for 198B

Direction is reliable only where the effect is large (D3) and for non-palindromic SNPs.

## 6. Verdict
**Keep the per-track score.** It is the right measure for comparing tissues, recommended by the
developers, and robust. A different way of comparing the ref and alt profiles isn't needed.

**Change what surrounds it, before any large run:**
- (A) switch to the recommended weights and measure the difference on the 40 scored SNPs
- (B) test input length on a GPU
- (C) save Σ REF / Σ ALT
- (D) normalise against a pooled matched null with a robust scale
- (E) add a null-calibrated "strongest track" statistic next to the category mean

A and C are code changes to `run_ism_pilot.py` that cost nothing at run time. A also needs
one rerun of the 40 SNPs on the Mac.