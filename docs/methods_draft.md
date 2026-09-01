# Methods Draft
**AlphaGenome Variant Effect Project · Aryana Satheesh · Luo Lab**

This document contains thesis-ready methods prose, built incrementally from the lab notebook after each completed phase. It is not meant to be complete until the project is finished — sections are added as work is done. Each section includes the biological reasoning behind methodological choices, so they can be adapted directly into a thesis or paper methods section.

---

## Completed sections

### 1. GWAS data acquisition

#### 1.1 Choice of GWAS dataset

We used summary statistics from the Psychiatric Genomics Consortium wave 3 schizophrenia GWAS (PGC3; Trubetskoy et al., 2022, Nature; PubMed 35396580), the largest schizophrenia GWAS conducted to date. The PGC3 meta-analysis combined data from multiple international cohorts and identified 270 distinct risk loci across ancestries. Summary statistics were downloaded from the PGC data repository (https://doi.org/10.6084/m9.figshare.19426775).

We used the European ancestry subset (`PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz`), which includes approximately 53,386 cases and 77,258 controls (effective sample size ≈ 58,749). This dataset contains 7,659,767 SNPs, of which 20,457 reach genome-wide significance (p < 5 × 10⁻⁸).

An earlier version of this analysis used the 2018 CLOZUK+PGC2 GWAS (scz2018clozuk; 8,064,800 SNPs, 18,088 genome-wide significant), but was updated to PGC3 to take advantage of the larger sample size and correspondingly greater number of discovered risk loci (173 independent loci in European after our filtering, compared to 141 from the 2018 data).

#### 1.2 Choice of ancestry subset

The European ancestry subset was chosen for several reasons. First, linkage disequilibrium (LD) patterns — which determine how we define independent loci through clumping — are population-specific, and using a single-ancestry analysis avoids confounding LD structures from different populations in the clumping procedure. Second, downstream reference panels and functional annotations are most complete for European ancestry data. Third, the European subset is directly comparable to our earlier analysis using the 2018 CLOZUK+PGC2 GWAS (also European ancestry), enabling cross-study validation of the pipeline.

Importantly, the ISM analysis itself is ancestry-agnostic: AlphaGenome is a sequence-to-function model that predicts chromatin accessibility from DNA sequence alone, without any population-genetic information. The ancestry consideration affects only which variants we select as inputs (through LD-dependent clumping), not the model's predictions for those variants. The multi-ancestry primary analysis (7,585,077 SNPs; 21,723 genome-wide significant) is available for future expansion.

#### 1.3 File format

The PGC3 summary statistics are provided in PGCsumstatsVCFv1.0 format, a tab-separated file with VCF-style `##` comment headers. Key columns are: CHROM (chromosome), ID (rsID), POS (base pair position), A1 (effect/alternate allele), A2 (other/reference allele), BETA (log odds ratio), SE (standard error), PVAL (p-value), NCAS (number of cases), NCON (number of controls), and NEFF (effective sample size). This differs from the 2018 CLOZUK file format in column naming conventions (e.g., BETA vs OR, PVAL vs P) but encodes equivalent information.

#### 1.4 Coordinate system verification

Genomic coordinates in the PGC3 summary statistics are in GRCh37/hg19. This was confirmed by cross-referencing the position of a known SNP (rs2007044) between the GWAS file and the Ensembl coREST API. The GWAS reports rs2007044 at chr12:2,344,960, which matches the Ensembl GRCh37 position exactly. The GRCh38/hg38 position for this SNP is chr12:2,235,794 — a discrepancy of ~109 kb, characteristic of build differences. The same verification was performed in our earlier analysis of the 2018 GWAS and the hg19 coordinate system is consistent across PGC releases.

### 2. MHC dominance and the need for structured SNP selection

A naive approach to selecting variants of interest — simply taking the top N SNPs ranked by p-value — fails dramatically for schizophrenia GWAS due to the overwhelming signal from the Major Histocompatibility Complex (MHC) region on chromosome 6. We systematically documented this before applying any filtering, as follows.

When all 20,457 genome-wide significant SNPs in the European PGC3 analysis are ranked by p-value, all 50 of the top 50 SNPs fall on chromosome 6 within a ~1.3 Mb window (positions ~27.5–28.8 Mb, hg19). The first non-chromosome-6 SNP does not appear until rank #1,169 (rs58120505, chr7:2,029,867, p = 2.235 × 10⁻²⁴). Of the 20,457 genome-wide significant SNPs, 6,104 (29.8%) are on chromosome 6, and 5,994 of those fall specifically within the extended MHC region (chr6:25–34 Mb).

This concentration arises from several properties of the MHC. The region is one of the most gene-dense and polymorphic in the human genome, encoding the HLA genes critical for immune function. LD blocks in the MHC can extend soacross megabases — far longer than the typical ~100–200 kb in European populations — meaning thousands of SNPs are correlated with each other, all reflecting a small number of underlying causal signals. The MHC does harbor genuine schizophrenia risk variants (the complement component C4 locus is a well-characterized example), but the region's extreme LD structure makes it nearly impossible to fine-map causal variants using standard approaches. For an ISM-based analysis, running AlphaGenome on thousands of tightly correlated MHC SNPs would be computationally wasteful and scientifically uninformative — the model would predict nearly identical effects for SNPs that are in near-perfect LD.

We therefore excluded the extended MHC region (chr6:25,000,000–34,000,000 in hg19) from lead SNP selection. This window is deliberately conservative, extending beyond the classical MHC boundaries (~28–33 Mb) to capture the full extent of long-range LD in the region. The top 50 SNPs (pre-exclusion) are preserved in a separate output file for documentation.

### 3. Lead SNP selection

After excluding the MHC, 14,463 genome-wide significant single-nucleotide variants remained. Insertion-deletion variants were removed prior to this step, as AlphaGenome accepts only single-nucleotide substitutions (the model predicts the effect of swapping one base for another at a specific position in a 131,072 bp sequence window; multi-base insertions or deletions cannot be represented in this framework). In the PGC3 European dataset, no indels were present among the genome-wide significant variants, as the VCF format restricts to biallelic SNPs.

To define independent risk loci, we applied greedy distance-based clumping. SNPs were sorted by p-value in ascending order (most significant first). The most significant SNP was selected as the lead SNP for the first locus. All remaining SNPs within 1 Mb on the same chromosome were then removed from consideration, and the process was repeated until no unassigned SNPs remained. The 1 Mb clumping window is a widely used convention in post-GWAS analysis (employed by PGC, PLINK, and standard GWAS pipelines) that approximately captures the extent of LD in European populations: SNPs separated by more than 1 Mb are unlikely to be in strong linkage disequilibrium and can be treated as independent signals.

This procedure yielded **173 independent lead SNPs** distributed across the autosomal genome. For initial pipeline development and testing, the 10 most significant lead SNPs (by p-value) were selected as a pilot set. These span 9 chromosomes (chr1, chr2, chr3, chr4, chr7, chr8, chr10, chr12, chr15), confirming that the MHC exclusion and clumping successfully produce a genome-wide distribution rather than a regionally clustered set.

### 4. Coordinate liftover (hg19 → hg38)


AlphaGenome was trained on GRCh38/hg38 genome sequences and requires hg38 coordinates as input. Because the PGC3 summary statistics use hg19 coordinates (Section 1.4), all 173 lead SNP positions were converted from hg19 to hg38 before use with the model.

Coordinate liftover was performed using pyliftover (v0.4.1), a Python implementation of the UCSC liftOver algorithm. The UCSC hg19ToHg38.over.chain.gz chain file encodes the mapping between genome builds: for each contiguous block of sequence in hg19, the chain file specifies the corresponding position in hg38. Positions shift between builds because hg38 incorporated improved assemblies of centromeric and telomeric regions, closed sequence gaps, fixed errors, and added alternate locus representations. These changes cause position offsets that vary by genomic region — from a few base pairs to several megabases — which is why a systematic liftover is necessary rather than applying a uniform offset.

A technical note on coordinate conventions: GWAS summary statistics and most genomics text files use 1-based coordinates (the first base of a chromosome is position 1). Pyliftover, following the BED format convention, uses 0-based half-open intervals internally. Our liftover script converts from 1-based to 0-based before querying pyliftover, then converts the result back to 1-based for the output files.

All 173 lead SNPs were successfully lifted with zero failures. Lifted coordinates were validated by two independent methods:

1. **Cross-study validation.** rs4129585 (chr8) appears in both the PGC3 pilot set and our earlier 2018 CLOZUK pilot set. Its hg38 position from both independent liftovers is identical (chr8:142,231,572), confirming consistency.

2. **Direct database verification.** The hg38 position of rs58120505 (the most significant non-MHC SNP, not present in any prior analysis) was checked against the Ensembl GRCh38 REST API. The API reports chr7:1,990,232, matching our lifted position exactly. Alleles (T/A/C/G in dbSNP, with our ref=C and alt=T among them) are also consistent.

# Methods Draft
**AlphaGenome Variant Effect Project · Aryana Satheesh · Luo Lab**
 
This document contains thesis-ready methods prose, built incrementally from the lab notebook after each completed phase. It is not meant to be complete until the project is finished — sections are added as work is done. Each section includes the biological reasoning behind methodological choices, so they can be adapted directly into a thesis or paper methods section.
 
---
 
## Completed sections
 
### 1. GWAS data acquisition
 
#### 1.1 Choice of GWAS dataset
 
We used summary statistics from the Psychiatric Genomics Consortium wave 3 schizophrenia GWAS (PGC3; Trubetskoy et al., 2022, Nature; PubMed 35396580), the largest schizophrenia GWAS conducted to date. The PGC3 meta-analysis combined data from multiple international cohorts and identified 270 distinct risk loci across ancestries. Summary statistics were downloaded from the PGC data repository (https://doi.org/10.6084/m9.figshare.19426775).
 
We used the European ancestry subset (`PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz`), which includes approximately 53,386 cases and 77,258 controls (effective sample size ≈ 58,749). This dataset contains 7,659,767 SNPs, of which 20,457 reach genome-wide significance (p < 5 × 10⁻⁸).
 
An earlier version of this analysis used the 2018 CLOZUK+PGC2 GWAS (scz2018clozuk; 8,064,800 SNPs, 18,088 genome-wide significant), but was updated to PGC3 to take advantage of the larger sample size and correspondingly greater number of discovered risk loci (173 independent loci in European after our filtering, compared to 141 from the 2018 data).
 
#### 1.2 Choice of ancestry subset
 
The European ancestry subset was chosen for several reasons. First, linkage disequilibrium (LD) patterns — which determine how we define independent loci through clumping — are population-specific, and using a single-ancestry analysis avoids confounding LD structures from different populations in the clumping procedure. Second, downstream reference panels and functional annotations are most complete for European ancestry data. Third, the European subset is directly comparable to our earlier analysis using the 2018 CLOZUK+PGC2 GWAS (also European ancestry), enabling cross-study validation of the pipeline.
 
Importantly, the ISM analysis itself is ancestry-agnostic: AlphaGenome is a sequence-to-function model that predicts chromatin accessibility from DNA sequence alone, without any population-genetic information. The ancestry consideration affects only which variants we select as inputs (through LD-dependent clumping), not the model's predictions for those variants. The multi-ancestry primary analysis (7,585,077 SNPs; 21,723 genome-wide significant) is available for future expansion.
 
#### 1.3 File format
 
The PGC3 summary statistics are provided in PGCsumstatsVCFv1.0 format, a tab-separated file with VCF-style `##` comment headers. Key columns are: CHROM (chromosome), ID (rsID), POS (base pair position), A1 (effect/alternate allele), A2 (other/reference allele), BETA (log odds ratio), SE (standard error), PVAL (p-value), NCAS (number of cases), NCON (number of controls), and NEFF (effective sample size). This differs from the 2018 CLOZUK file format in column naming conventions (e.g., BETA vs OR, PVAL vs P) but encodes equivalent information.
 
#### 1.4 Coordinate system verification
 
Genomic coordinates in the PGC3 summary statistics are in GRCh37/hg19. This was confirmed by cross-referencing the position of a known SNP (rs2007044) between the GWAS file and the Ensembl REST API. The GWAS reports rs2007044 at chr12:2,344,960, which matches the Ensembl GRCh37 position exactly. The GRCh38/hg38 position for this SNP is chr12:2,235,794 — a discrepancy of ~109 kb, characteristic of build differences. The same verification was performed in our earlier analysis of the 2018 GWAS and the hg19 coordinate system is consistent across PGC releases.
 
### 2. MHC dominance and the need for structured SNP selection
 
A naive approach to selecting variants of interest — simply taking the top N SNPs ranked by p-value — fails dramatically for schizophrenia GWAS due to the overwhelming signal from the Major Histocompatibility Complex (MHC) region on chromosome 6. We systematically documented this before applying any filtering, as follows.
 
When all 20,457 genome-wide significant SNPs in the European PGC3 analysis are ranked by p-value, all 50 of the top 50 SNPs fall on chromosome 6 within a ~1.3 Mb window (positions ~27.5–28.8 Mb, hg19). The first non-chromosome-6 SNP does not appear until rank #1,169 (rs58120505, chr7:2,029,867, p = 2.235 × 10⁻²⁴). Of the 20,457 genome-wide significant SNPs, 6,104 (29.8%) are on chromosome 6, and 5,994 of those fall specifically within the extended MHC region (chr6:25–34 Mb).
 
This concentration arises from several properties of the MHC. The region is one of the most gene-dense and polymorphic in the human genome, encoding the HLA genes critical for immune function. LD blocks in the MHC can extend across megabases — far longer than the typical ~100–200 kb in European populations — meaning thousands of SNPs are correlated with each other, all reflecting a small number of underlying causal signals. The MHC does harbor genuine schizophrenia risk variants (the complement component C4 locus is a well-characterized example), but the region's extreme LD structure makes it nearly impossible to fine-map causal variants using standard approaches. For an ISM-based analysis, running AlphaGenome on thousands of tightly correlated MHC SNPs would be computationally wasteful and scientifically uninformative — the model would predict nearly identical effects for SNPs that are in near-perfect LD.
 
We therefore excluded the extended MHC region (chr6:25,000,000–34,000,000 in hg19) from lead SNP selection. This window is deliberately conservative, extending beyond the classical MHC boundaries (~28–33 Mb) to capture the full extent of long-range LD in the region. The top 50 SNPs (pre-exclusion) are preserved in a separate output file for documentation.
 
### 3. Lead SNP selection
 
After excluding the MHC, 14,463 genome-wide significant single-nucleotide variants remained. Insertion-deletion variants were removed prior to this step, as AlphaGenome accepts only single-nucleotide substitutions (the model predicts the effect of swapping one base for another at a specific position in a 131,072 bp sequence window; multi-base insertions or deletions cannot be represented in this framework). In the PGC3 European dataset, no indels were present among the genome-wide significant variants, as the VCF format restricts to biallelic SNPs.
 
To define independent risk loci, we applied greedy distance-based clumping. SNPs were sorted by p-value in ascending order (most significant first). The most significant SNP was selected as the lead SNP for the first locus. All remaining SNPs within 1 Mb on the same chromosome were then removed from consideration, and the process was repeated until no unassigned SNPs remained. The 1 Mb clumping window is a widely used convention in post-GWAS analysis (employed by PGC, PLINK, and standard GWAS pipelines) that approximately captures the extent of LD in European populations: SNPs separated by more than 1 Mb are unlikely to be in strong linkage disequilibrium and can be treated as independent signals.
 
This procedure yielded **173 independent lead SNPs** distributed across the autosomal genome. For initial pipeline development and testing, the 10 most significant lead SNPs (by p-value) were selected as a pilot set. These span 9 chromosomes (chr1, chr2, chr3, chr4, chr7, chr8, chr10, chr12, chr15), confirming that the MHC exclusion and clumping successfully produce a genome-wide distribution rather than a regionally clustered set.
 
### 4. Coordinate liftover (hg19 → hg38)
 
AlphaGenome was trained on GRCh38/hg38 genome sequences and requires hg38 coordinates as input. Because the PGC3 summary statistics use hg19 coordinates (Section 1.4), all 173 lead SNP positions were converted from hg19 to hg38 before use with the model.
 
Coordinate liftover was performed using pyliftover (v0.4.1), a Python implementation of the UCSC liftOver algorithm. The UCSC hg19ToHg38.over.chain.gz chain file encodes the mapping between genome builds: for each contiguous block of sequence in hg19, the chain file specifies the corresponding position in hg38. Positions shift between builds because hg38 incorporated improved assemblies of centromeric and telomeric regions, closed sequence gaps, fixed errors, and added alternate locus representations. These changes cause position offsets that vary by genomic region — from a few base pairs to several megabases — which is why a systematic liftover is necessary rather than applying a uniform offset.
 
A technical note on coordinate conventions: GWAS summary statistics and most genomics text files use 1-based coordinates (the first base of a chromosome is position 1). Pyliftover, following the BED format convention, uses 0-based half-open intervals internally. Our liftover script converts from 1-based to 0-based before querying pyliftover, then converts the result back to 1-based for the output files.
 
All 173 lead SNPs were successfully lifted with zero failures. Lifted coordinates were validated by two independent methods:
 
1. **Cross-study validation.** rs4129585 (chr8) appears in both the PGC3 pilot set and our earlier 2018 CLOZUK pilot set. Its hg38 position from both independent liftovers is identical (chr8:142,231,572), confirming consistency.
2. **Direct database verification.** The hg38 position of rs58120505 (the most significant non-MHC SNP, not present in any prior analysis) was checked against the Ensembl GRCh38 REST API. The API reports chr7:1,990,232, matching our lifted position exactly. Alleles (T/A/C/G in dbSNP, with our ref=C and alt=T among them) are also consistent.
### 5. Sequence extraction and allele assignment
 
AlphaGenome accepts a fixed input of 131,072 bp of one-hot encoded DNA sequence. For each lead SNP, we extracted the hg38 sequence window centered on the variant position (65,536 bp on either side) via the Ensembl REST API. Sequences were one-hot encoded with the channel order A, C, G, T; ambiguous or unknown bases were encoded as all-zero vectors.
 
A methodological point that materially affects the analysis concerns allele orientation. In GWAS summary statistics, the two allele columns (here A1 and A2) denote *statistical roles* rather than positions relative to the reference assembly: A1 is the effect allele, whose association with disease risk is quantified by the effect size estimate, and A2 is the other allele. Neither is defined with reference to the genome build, and either may correspond to the reference genome base at that position. In the PGC3 European pilot set, the hg38 reference base matched A1 — the effect allele — at all ten positions.
 
For in silico mutagenesis this distinction is not cosmetic. The reference sequence supplied to the model must be the actual genomic sequence, and the alternate sequence must differ from it at exactly the variant position; supplying the two in the wrong orientation inverts the sign of every variant effect score. We therefore determined allele identity empirically rather than by column convention: for each variant, the reference allele was defined as whichever of A1 or A2 matched the hg38 base at that position, and the alternate allele as the other. Variants for which the reference base matched neither reported allele were flagged and excluded, as this would indicate a coordinate or strand problem rather than a labeling ambiguity; no pilot variants were excluded on this basis.
 
### 6. In silico mutagenesis and variant effect scoring
 
For each variant, the reference and alternate sequences were passed independently through the pretrained AlphaGenome model, yielding predicted signal at every position across the 131,072 bp window for each output track. Predictions were obtained at 1 bp resolution for chromatin accessibility outputs (ATAC and DNase); the histone ChIP output head is available only at 128 bp binned resolution.
 
Variant effect scores were computed following the aggregation recommended for chromatin accessibility outputs. A spatial mask of 501 bp centered on the variant was applied, restricting scoring to the local sequence neighborhood on the rationale that a variant disrupting a transcription factor binding site alters accessibility locally rather than across the full input window. Within this mask, predicted signal was summed separately for the reference and alternate sequences, and the effect score for each track computed as the log2 ratio of the two sums with a pseudocount of 1:
 
  score = log2[(Σ alt + 1) / (Σ ref + 1)]
 
The pseudocount avoids undefined values for tracks with no predicted signal at the locus. Positive scores indicate that the alternate allele is predicted to increase accessibility; negative scores indicate a predicted decrease. Because the score is on a log2 scale it is directly interpretable as fold change.
 
### 7. Output track annotation and padding
 
Interpreting per-track scores requires mapping each output channel to the biological sample it represents. Track metadata was obtained from the AlphaGenome research repository (`OutputMetadataResponse_ORGANISM_HOMO_SAPIENS.textproto`), which records, for each track, the assay, ontology term, biosample name, biosample type, developmental stage, and data source.
 
Two features of the output representation required attention before this metadata could be used.
 
First, **the output heads are zero-padded to round tensor dimensions.** The ATAC head contains 167 real tracks padded to 256 channels; DNase contains 305 real tracks padded to 384; the histone ChIP head contains 1,116 real tracks padded to 1,152. Padding entries are explicitly labeled as such in the metadata, carry empty biosample and assay fields, and occupy a contiguous block at the end of each head, so that real tracks occupy indices 0 through N−1 in metadata order. Padding channels carry no biological meaning and were excluded before computing any summary statistic; retaining them dilutes per-track means and inflates the apparent number of tracks. This also reconciles two track counts that appear inconsistent — the number of real tracks and the width of the output tensor — which differ only by the padding.
 
Second, **assigning tissue labels by substring matching on biosample names is error-prone** and was rejected in favor of a manually curated mapping. Naive matching produced false positives including "gastrocnemius medialis" and "gastroesophageal sphincter" (matching "astro"), and several renal samples described as "renal cortex" or "renal cortical" (matching "cortex"). Peripheral nerve samples such as sciatic and tibial nerve, which are not central nervous system tissue, were likewise excluded. The final classification was constructed by reading the complete set of biosample names present in each output head and assigning central nervous system, eye, and other categories by exact match.
 
Under this classification, central nervous system coverage differs markedly across output heads: 1 of 167 real ATAC tracks, 23 of 305 DNase tracks, 70 of 1,116 histone ChIP tracks, and 35 of 733 RNA-seq tracks. The ATAC head therefore provides essentially no basis for brain-specific inference in the pretrained model, and analyses of tissue specificity rely on the DNase and histone ChIP outputs. This limitation also motivates the fine-tuning component of the project, since fine-tuning on cell-type-resolved fetal brain ATAC-seq would supply precisely the coverage the pretrained ATAC head lacks.
 
Verification that metadata row order corresponds to output channel order was performed by inspecting the model implementation and the structure of the metadata file. Padding entries appear contiguously at the end of each head with real tracks preceding them in order, and the count of padding channels matches the count of output channels with no predicted signal variance at a test locus (79 for DNase), confirming that real tracks are indexed 0 through N−1 as recorded.
 
---
 
## Sections to be written
 
- Negative controls (external: non-brain trait GWAS; internal: matched non-significant variants)
- Tissue-specific scoring: central nervous system versus non-neural tracks
- Visualization approach (delta tracks, heatmaps)
- Analysis and interpretation of variant effects across cell types
- Discussion: in silico mutagenesis as a complement to statistical fine mapping under linkage disequilibrium
- *(Fine-tuning on Ziffra fetal brain ATAC-seq, if completed)*
---
 
*Last updated: August 25, 2026*
