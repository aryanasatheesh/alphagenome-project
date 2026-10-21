# Allele conventions: GWAS alleles, the hg38 reference, and risk direction

Every variant effect score is a comparison between two sequences that differ at
one base. To interpret the sign of a score you need to know which base is in
which sequence and which base raises disease risk. These are three separate
questions, and GWAS files answer only one of them.

## Three different labels for the same two bases

| Label | Defined by | Meaning |
|---|---|---|
| **A1 / A2** | the GWAS file | A1 is the **effect allele**: BETA is the log odds ratio for carrying A1. A2 is the other allele. |
| **hg38 ref / hg38 alt** | the reference genome | ref is the base printed in hg38 at that position; alt is the other GWAS allele, which we substitute in. |
| **risk / non-risk** | the sign of BETA | risk is A1 if BETA > 0, A2 if BETA < 0. |

A1 is a statistical role, not a genomic one. Nothing in the definition of A1
says whether it is the hg38 reference base. Assuming "A2 = reference" (as the
project's first scripts did) is a common and silent mistake: it produces
plausible-looking scores with the wrong sign.

## What we observed in PGC3

In all 20 SNPs checked (10 pilot, 10 control), **the hg38 reference base was
A1**. Seventeen of the 20 are non-palindromic, and a non-palindromic SNP
reported on the minus strand would have matched neither allele. So PGC3 reports
alleles on the + strand, and in this file A1 appears to be the reference allele.

The code does not rely on this. `run_ism_pilot.py` reads the hg38 base at every
position and records which GWAS allele it is (`hg38_ref_is`). A future SNP where
hg38 carries A2 is handled correctly and is easy to spot.

## Columns, by file

**Selection outputs** (`phase1/data/scz_pgc3_*.tsv` and their `_hg38` versions)

| Column | Meaning |
|---|---|
| `a1` | GWAS effect allele (PGC3 `A1`) |
| `a2` | GWAS other allele (PGC3 `A2`) |
| `a1_freq_controls` | frequency of A1 in GWAS controls (PGC3 `FCON`) |
| `beta` | log odds ratio for A1 (PGC3 `BETA`) |

**ISM outputs** (`results/ism/*_ism_scores.tsv`, one row per SNP × track)

| Column | Meaning |
|---|---|
| `hg38_ref` | base in hg38 at the variant; present in the reference sequence |
| `hg38_alt` | the other GWAS allele; substituted into the alternate sequence |
| `hg38_ref_is` | `A1` or `A2`: which GWAS allele hg38 carries |
| `palindromic` | `True` for A/T and C/G SNPs (see below) |
| `variant_effect_score` | log2(signal with hg38_alt ÷ signal with hg38_ref) |
| `risk_oriented_score` | log2(signal with risk allele ÷ signal with non-risk allele); blank if BETA is 0 or missing |

`results/ism/*_ism_scores.variants.tsv` has one row per SNP with the same allele
columns plus `risk_allele`, and `status`/`note` saying whether the SNP was scored
or skipped and why.

## The two scores

Both scores are computed over a 501 bp window centred on the variant, from
predicted signal summed within that window with a pseudocount of 1:

```
variant_effect_score = log2[(Σ signal with hg38_alt + 1) / (Σ signal with hg38_ref + 1)]
```

`variant_effect_score` answers: *what happens to accessibility if the hg38 base
is replaced by the other allele?* Positive means the substitution opens chromatin.

`risk_oriented_score` answers the biological question: *does the schizophrenia
risk allele open or close chromatin?* Positive means the risk allele opens it.

```
risk_oriented_score = +variant_effect_score   if hg38_alt is the risk allele
                      −variant_effect_score   if hg38_ref is the risk allele
```

Equivalently, `sign(BETA) × variant_effect_score × (+1 if hg38_alt = A1, −1 if hg38_alt = A2)`.

**Worked example: rs58120505** A1 = T, A2 = C, BETA = +0.0896, so T is the
risk allele. hg38 has T at this position, so hg38_ref = T (A1) and
hg38_alt = C (A2). The risk allele is the reference base, so
`risk_oriented_score = −variant_effect_score`. If substituting C *raised*
predicted accessibility (positive `variant_effect_score`), the risk allele T
*lowers* it (negative `risk_oriented_score`).

**Which score to use.** Magnitude analyses (is this variant's effect larger in
brain than elsewhere?) use `|variant_effect_score|`; the two scores have the
same absolute value. Any statement about direction ("the risk allele reduces
accessibility in neurons") must use `risk_oriented_score`.

## Palindromic SNPs

For A/T and C/G SNPs the two alleles are each other's complements, so matching
against hg38 cannot tell whether the GWAS reported the + or − strand: either
way, one allele matches.

- **Magnitude is unaffected.** The substituted base is the one hg38 does not
  carry, whatever the strand, so the model sees the same two sequences.
- **Direction assumes the + strand.** If a palindromic SNP had been reported on
  the − strand, `hg38_ref_is`, `risk_allele` and the sign of
  `risk_oriented_score` would all be reversed.

How common: 15.1% of all PGC3 SNPs and 14.8% of genome-wide significant SNPs
are palindromic. In the current sets: rs1198588 (pilot) and rs2067987,
rs7326762 (controls).

Why we accept the + strand assumption: all 17 non-palindromic SNPs checked
matched on the + strand, and a mixed-strand file would be unusual.

How to verify it: compare `a1_freq_controls` with the frequency of the hg38
reference allele in a European reference panel (e.g. gnomAD non-Finnish
European). On the + strand the two should agree; on the − strand A1 would be
the other allele, with frequency near 1 − (reference frequency). The test
cannot decide SNPs with allele frequency near 0.5. Not yet done.

## Interaction with liftover

`liftover_to_hg38.py` rejects SNPs whose region maps to the − strand in hg38
(8,304 PGC3 SNPs, 0.108%). For a palindromic SNP in such a region the hg38 base
is the complement of the hg19 base, which is the *other* GWAS allele, so the
allele check would pass and silently reverse the direction. Rejecting strand
flips at liftover closes that hole. 1,383 of those 8,304 SNPs (16.7%) are
palindromic, close to the 15.1% genome-wide rate.

## Where this is enforced

| Rule | Code |
|---|---|
| A1/A2 written as GWAS roles, never as ref/alt | `get_lead_snps_pgc3.py`, `get_control_snps.py` |
| every output column present (no silent blanks) | same two scripts, guard before writing |
| strand flips rejected | `liftover_to_hg38.py`, `lift_one()` |
| hg38 base read and matched to A1/A2 per SNP; skipped if neither | `run_ism_pilot.py`, `assign_alleles()` |
| risk direction and palindromic flag recorded | `run_ism_pilot.py`, `assign_alleles()` |
| old `ref`/`alt` input files refused | `run_ism_pilot.py`, before the model loads |

## History

Until 16 October 2026 the selection scripts wrote A2 to a column named `ref` and
A1 to a column named `alt`. `run_ism_pilot.py` corrected this at run time, so
every score was correct, but the files were misleading. Scores produced before
and after the change were compared track by track: all 12,800 pilot and control
scores are identical, and the old `gwas_alleles_swapped = True` flag corresponds
exactly to `hg38_ref_is = A1`.