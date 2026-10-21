# Lifting GWAS SNPs from hg19 to hg38

**Scripts:** `phase1/liftover_to_hg38.py` (the pipeline step) and
`phase1/survey_liftover_edge_cases.py` (a whole-genome audit of that step).

## Why this step exists

PGC3 schizophrenia summary statistics (Trubetskoy et al., 2022, *Nature* 604:502–508)
report positions on **hg19 (GRCh37)**. AlphaGenome was trained on **hg38 (GRCh38)**.
Every SNP must be moved to hg38 before its sequence is extracted, or the model
would be shown the wrong stretch of DNA.


## Where it sits in the pipeline

```
PGC3 summary statistics (hg19)
  └─ get_lead_snps_pgc3.py / get_control_snps.py   select SNPs, write hg19 TSV
       └─ liftover_to_hg38.py                       hg19 → hg38 (this step)
            └─ run_ism_pilot.py                     fetch hg38 sequence, score
```

Liftover is deliberately a separate step: every hg38 file in the project is
produced by this one script, with a report and provenance log beside it.

## Setup (once per machine)

```bash
pip install pyliftover            # validated with 0.4.1 on Mac and Hoffman2
curl -o phase1/data/hg19ToHg38.over.chain.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
```

The chain file (~222 KB) is gitignored. The script looks for it in
`phase1/data/` relative to its own location, so it works from any directory.

pyliftover is used instead of UCSC's `liftOver` binary because the binary
requires a newer GLIBC than Hoffman2 provides. pyliftover reimplements
liftOver's conversion of single positions and reads the same UCSC chain file.

## Usage

```bash
python phase1/liftover_to_hg38.py \
    --input  phase1/data/scz_pgc3_lead_snps.tsv \
    --output phase1/data/scz_pgc3_lead_snps_hg38.tsv
```

Options: `--chain PATH`, `--chr-col NAME` (default `chr`), `--pos-col NAME`
(default `pos`), `--id-col NAME` (default `rsid`), `--report PATH`,
`--allow-failures`.

**Input:** a tab-separated file with a header. Chromosome as `6` or `chr6`;
position 1-based (the convention of GWAS files, VCF, dbSNP and Ensembl).

**Outputs:**

| File | Contents |
|---|---|
| `<output>.tsv` | every input column, with `chr`/`pos` replaced by hg38 values, plus `chr_hg19` and `pos_hg19` |
| `<output>.liftover_report.tsv` | one row per input SNP: status, number of mappings, hg38 position, strand, note |
| `<output>.liftover_log.txt` | command, date, chain file name and md5, pyliftover and Python versions, counts |

Allele columns are carried through untouched. Liftover moves coordinates only;
allele orientation is handled downstream (see `docs/allele_conventions.md`).

## Strict by default

A SNP is written to the output only if it maps to **exactly one** place, on a
**standard chromosome** (1–22, X, Y), on the **same chromosome** it came from,
on the **+ strand**.

If any SNP fails, the script writes the report and log but NOT the output
file, and exits with status 1. `--allow-failures` writes the
lifted rows anyway; failures stay listed in the report.

Two further guards:

- It refuses a file that already has `chr_hg19`/`pos_hg19` columns. Lifting
  twice would shift every SNP silently.
- It refuses to overwrite its own input.

## Coordinate convention

pyliftover works in 0-based coordinates. The script subtracts 1 before calling
it and adds 1 to the result, so input and output are both 1-based. This
off-by-one is the most likely way to corrupt coordinates silently, so it lives
in one function (`lift_one`) and nowhere else.

## Edge cases, and how often they occur in PGC3

The table counts every SNP in `PGC3_SCZ_european.vcf.tsv.gz` (7,659,767 SNPs),
produced by `survey_liftover_edge_cases.py` using the same `lift_one` function
as the pipeline (`results/liftover_survey/pgc3_european.summary.tsv`).

| Status | What happened | Why the SNP is not passed on | All SNPs | p < 5×10⁻⁸ | p > 0.5 |
|---|---|---|---|---|---|
| `lifted` | one mapping, same chromosome, + strand | — | 7,646,680 (99.829%) | 20,434 (99.888%) | 3,099,610 (99.849%) |
| `strand_flip` | maps to the − strand | the hg38 base is the complement of the hg19 base, so GWAS alleles no longer match the reference strand | 8,304 (0.108%) | 12 | 2,712 |
| `unmapped` | region has no hg38 equivalent | there is no coordinate to use | 3,809 (0.050%) | 10 | 1,559 |
| `alt_contig` | lands on an alt/random/unplaced contig | not part of the primary assembly AlphaGenome was trained on | 948 (0.012%) | 1 | 418 |
| `chrom_change` | lands on a different chromosome | almost always an assembly artefact for a single SNP | 26 (<0.001%) | 0 | 11 |
| `multi` | maps to several places | no way to choose; the old code silently took the first | 0 | 0 | 0 |
| `unknown_chrom`, `invalid_input` | chromosome not in chain / unparseable input | — | 0 | 0 | 0 |
| **failures** | | | **13,087 (0.171%)** | **23 (0.112%)** | **4,700 (0.151%)** |

What the table shows:

- **Strand flips are the most common failure** (63% of failures), more common
  than unmapped regions. For non-palindromic SNPs a strand flip would make the
  hg38 base match neither GWAS allele, and `run_ism_pilot.py` would stop. For
  palindromic (A/T, C/G) SNPs it would instead match the *wrong* allele and
  silently reverse the risk direction. Rejecting strand flips here prevents that.
  1,383 of the 8,304 strand-flip SNPs (16.7%) are palindromic, close to the
  15.1% rate across all PGC3 SNPs, so roughly 1 in 6 strand flips is a SNP
  whose direction would have been silently reversed.
- **No SNP maps to multiple places** with this chain file. The check costs
  nothing and protects other datasets or chain files.
- **None of the 173 lead SNPs, 10 pilot SNPs or 10 control SNPs hit any edge
  case.** Of the 23 genome-wide significant SNPs that fail:
  - **12 lie in the excluded MHC**, all strand flips, clustered at
    chr6:26.74–26.76 Mb (hg19).
  - **11 lie within 1 Mb of a lead SNP that lifted**, in three loci: 8 near
    rs5751191 (chr22:42.53–42.55 Mb; 7 unmapped, 1 alt contig), 2 near
    rs12489270 (chr3:16.86 Mb) and 1 near rs12138231 (chr1:150.15 Mb).
  - **None lies outside the MHC and more than 1 Mb from a lead**, so no
    schizophrenia locus is lost. The 11 would matter only if we later scored
    every SNP in a locus instead of its lead SNP.
- **Planning for more controls:** with 0.151% of null SNPs failing, a random
  draw of 50 controls has about a 7% chance of containing at least one SNP that
  cannot be lifted, and 100 controls about 14%. Strict mode will then stop and
  name it; the fix is to replace that control, not to relax the check.

The survey takes about 70 seconds on a Hoffman2 compute node:

```bash
qrsh -l h_data=4G,h_rt=1:00:00
python phase1/survey_liftover_edge_cases.py \
    --gwas data/PGC3_SCZ_european.vcf.tsv.gz \
    --out-prefix results/liftover_survey/pgc3_european
```

## What liftover does *not* check

**Whether the reference base itself changed between builds.** At some positions
hg38 corrected the hg19 reference base. Liftover moves the coordinate
correctly but cannot see bases. `run_ism_pilot.py` catches this when it reads
the hg38 base: a non-palindromic SNP whose hg38 base matches neither allele is
skipped and listed; one that matches A2 instead of the usual A1 shows up in the
`hg38_ref_is` column.

**Rescue of strand-flipped SNPs.** They could be recovered by
reverse-complementing their alleles. We exclude them instead because none
affects the current SNP sets, and a rescue would add a correction step that
would itself need validating.

## Validation

1. **Synthetic tests.** A hand-built chain file and 13 test SNPs covering each
   status, both edges of a mapped block, the first base of an unmapped gap, and
   invalid input (chromosome 23, non-numeric position, position 0). Every
   expected status and hg38 position was produced. The guards (re-lift
   refusal, overwrite refusal, missing column, missing chain file) were each
   triggered and each stopped the script.
2. **Reproduction of the summer 2026 files.** Re-lifting the hg19 lead, pilot
   and control files gave identical hg38 positions for 173/173, 10/10 and 10/10
   SNPs. The August lead-SNP liftover had no committed script; it is now
   reproducible.
3. **External spot checks (summer 2026).** rs58120505 → chr7:1,990,232 and
   rs4129585 → chr8:142,231,572 match Ensembl GRCh38.
4. **Whole-genome survey** (table above): 7,659,767 SNPs read, matching the
   row count of the source file; 20,457 genome-wide significant SNPs, matching
   the count from lead-SNP selection.
5. **Real-chain self-test** (`tests/liftover/check_liftover.py`): lifts six real
   PGC3 SNPs with known outcomes (two lifted to known hg38 positions; one each
   unmapped, alt contig, strand flip and invalid) and triggers every guard.
   Anyone can run it after downloading the chain file; it should print
   `ALL CHECKS PASSED`.

## Testing it yourself

Step-by-step instructions for a first-time user, including how to send
feedback, are in [`docs/TESTING_liftover.md`](TESTING_liftover.md).