# Methods Draft
**AlphaGenome Variant Effect Project · Aryana Satheesh · Luo Lab**

This document contains thesis-ready methods prose, built incrementally from the lab notebook after each completed phase. It is not meant to be complete until the project is finished — sections are added as work is done.

---

## Completed sections

### 1. GWAS data acquisition and processing

Schizophrenia GWAS summary statistics were obtained from the CLOZUK+PGC2 meta-analysis (scz2018clozuk; N ≈ 105,000 cases and controls), downloaded from the Psychiatric Genomics Consortium (PGC). The file contained 8,064,800 SNPs, of which 18,088 reached genome-wide significance (p < 5 × 10⁻⁸). Coordinates in the summary statistics were in GRCh37/hg19, confirmed by cross-referencing rsID positions against the Ensembl GRCh37 REST API for all five pilot variants.

To obtain a set of independent risk loci, we applied the following filtering procedure in Python. First, SNPs were restricted to those reaching genome-wide significance (p < 5 × 10⁻⁸). Insertion-deletion variants were removed, retaining only single-nucleotide substitutions, as AlphaGenome accepts only single-base substitutions as input. SNPs were then sorted by p-value in ascending order and a 1 Mb clumping window was applied per chromosome: the most significant SNP was selected as the lead SNP for each locus, and all remaining SNPs within 1 Mb on the same chromosome were excluded before selecting the next lead. The major histocompatibility complex region (chromosome 6, ~28–32 Mb) was excluded due to its complex linkage disequilibrium structure. This procedure yielded 141 independent lead SNPs across the genome. A subset of five pilot SNPs was selected for initial pipeline development, corresponding to the five most significant lead SNPs outside the MHC region.

### 2. Coordinate liftover (hg19 → hg38)

AlphaGenome was trained on and requires GRCh38/hg38 sequences. Since the GWAS summary statistics used hg19 coordinates, all 141 lead SNP positions were lifted over to hg38 using pyliftover (v0.4.1), a Python implementation of the UCSC liftOver algorithm, with the UCSC hg19ToHg38.over.chain.gz chain file. Lifted coordinates were validated by cross-referencing the five pilot SNP positions against the Ensembl GRCh38 REST API; all five matched exactly. All 141 SNPs were successfully lifted with zero failures.

---

## Sections to be written

- Sequence extraction (131,072 bp hg38 windows, ref/alt construction)
- ISM inference pipeline (variant scoring, mask, aggregation)
- Track selection and brain vs non-brain comparison
- Visualization approach
- *(fine-tuning sections, if completed)*

---

*Last updated: July 21, 2026*
