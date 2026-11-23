# Weights decision: fold 0 vs all-folds (Nov 23)

**Runs:** 40 SNPs (10 pilot, 10 random control, 20 T1D) rescored with `model_all_folds.safetensors` (md5 f26622eb9042b251789b048119d6d4ab), commit ecd96da, Mac/MPS, 131 kb. Commit 451a115. Fold 0 md5 cb18d831d374e6dab5758a39021a719b.

**Rule fixed before looking:** median Spearman ≥ 0.95, sign agreement ≥ 0.95, pilot conclusions unchanged → fold-0 results stand. Otherwise switch and redo the Session 17 analyses.

**Result: rule not met.** Median per-SNP Spearman across tracks is 0.75–0.90 in every set and modality. Sign agreement medians are 0.96–1.00, but the minimum is 0.41.
- Overall scale is unchanged: median ratio of summed |score| (all-folds/fold 0) is 1.01 for ATAC and 0.98 for DNase. Slopes below 1 come from regression dilution (both runs are noisy), not from a scale shift.
- The disagreement is real, not just ranking of near-zero tracks. Restricted to tracks with |score| ≥ 0.05, Spearman is still about 0.8. Strong single-track effects can change a lot: rs1198588 globus pallidus is 0.53 → 0.08, and rs61759532 ATAC has ρ = −0.12 over 157 tracks.
- DNase denominators D_t: Spearman 0.92 between the weight sets. The immune ruler is 1.54× the median track (1.57× with fold 0).

**Session 17 analyses redone with all-folds**
- Ruler split-half reliability: 0.34 (10 SNPs) → 0.51 (20 SNPs). Fold 0 was 0.35 → 0.52, so this replicates.
- T1D vs random-control ruler: Spearman −0.37 (fold 0 −0.18). T1D is still not a neutral ruler.
- Size-1 threshold: still not discriminating. Random controls still land in "shared" and "other_tissue".
- Specificity ranking across the 20 pilot + control SNPs: Spearman 0.96 between weight sets. Profile class is the same for 17 of 20; top body category for 13 of 20.
- **SNP-level change:** rs1198588 is no longer brain-leaning. Specificity was 2.97 (CI 1.80–6.05) and is now 0.42 (CI 0.22–0.99), and its class went from shared to none. Raw CNS mean |score| fell from 0.219 to 0.049, so the change is in the model output, not the ruler. rs2238057 (CACNA1C) rose from 0.48 to 1.29 (CI 0.65–2.89).

**Decision**
- Switch the default to all-folds (`run_ism_pilot.py` md5 67fb0fa39a54060dc5e282d18677892f).
- The methodological conclusions from Session 17 stand.
- SNP-level fold-0 results are superseded and stay labelled "fold 0".

**Implication for the pipeline:** model choice is a substantial source of track-level variation. Category averaging reduces it, but single-SNP calls can flip. A model-robustness check is needed, either for flagged SNPs or for all of them: rescore with folds 0–3 and report agreement, or use the fold ensemble. Question for Choo.

**Runtime:** 13 s/SNP on the Mac with cached sequences at 131 kb. That is about 40 min for all 173 leads, or about 2 h for 500 pool SNPs per model. A GPU is not needed at this scale; it is still needed for the 1 Mb test.