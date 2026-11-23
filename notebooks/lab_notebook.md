# AlphaGenome Project — Lab Notebook
**Researcher:** Aryana Satheesh  
**PI:** Prof. Chongyuan Luo  
**Mentor:** Cuining (Choo) Liu

This notebook documents every work session in chronological order: what was run, what the output was, decisions made, and any blockers. After each completed phase, these notes get converted into thesis-ready methods prose. 

---
 
## Session 1
**Goal:** Get AlphaGenome running locally on Mac (Task 1 of summer plan)

### Context
AlphaGenome was previously installed and working on Hoffman2 (alphagenome_env, Python 3.12, CUDA 11.8). This session sets up a local Mac environment for faster development and testing. The Mac workflow is: develop + test locally → push to GitHub → run at scale on Hoffman2.

### Environment setup

**System:** MacBook Pro, Apple Silicon (arm64)
**Python available before setup:** 3.9.6 (system) — too old, alphagenome-pytorch requires 3.12+

**Step 1: Install Homebrew**
```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
# then run the three PATH commands printed by the installer:
echo >> /Users/aryana/.zprofile
echo 'eval "$(/opt/homebrew/bin/brew shellenv zsh)"' >> /Users/aryana/.zprofile
eval "$(/opt/homebrew/bin/brew shellenv zsh)"
```

**Step 2: Install Python 3.12**
```bash
brew install python@3.12
# installed: Python 3.12.13 at /opt/homebrew/bin/python3.12
```

**Step 3: Create virtual environment and install alphagenome-pytorch**
```bash
python3.12 -m venv ~/alphagenome-env
source ~/alphagenome-env/bin/activate
pip install alphagenome-pytorch
```
Installed cleanly. Key packages: alphagenome-pytorch==0.3.1, torch==2.13.0, numpy==2.5.1

**To reactivate this environment in future sessions (locally):**
```bash
source ~/alphagenome-env/bin/activate
```

**Step 4: Download model weights**
```bash
pip install huggingface_hub
hf download gtca/alphagenome_pytorch model_fold_0.safetensors --local-dir ~/alphagenome-weights
# downloaded 921 MB to: /Users/aryana/alphagenome-weights/model_fold_0.safetensors
```
Note: there are 4 folds total (model_fold_0 through model_fold_3). We download fold 0 only for now — sufficient for development and pilot runs. Ensemble across all 4 folds would be used for final results.

**Step 5: Smoke test — forward pass**
```bash
PYTORCH_ENABLE_MPS_FALLBACK=1 python3 - <<'EOF'
import torch
import numpy as np
from alphagenome_pytorch import AlphaGenome

device = "mps" if torch.backends.mps.is_available() else "cpu"
print(f"using device: {device}")

model = AlphaGenome.from_pretrained(
    "/Users/aryana/alphagenome-weights/model_fold_0.safetensors",
    device=device
)
model.eval()
print("model loaded successfully")

seq = np.random.randint(0, 4, size=(1, 131072))
dna = torch.tensor(np.eye(4)[seq], dtype=torch.float32).to(device)

with torch.inference_mode():
    outputs = model.predict(dna, organism_index=0)  # 0 = human, 1 = mouse

print("forward pass OK")
print("output keys:", list(outputs.keys()))
print("atac shape:", outputs['atac'][1].shape)
EOF
```

### Output
```
using device: mps
model loaded successfully
[UserWarning: aten::logspace.out not supported on MPS, falling back to CPU — harmless, handled by MPS_FALLBACK]
forward pass OK
output keys: ['atac', 'dnase', 'procap', 'cage', 'rna_seq', 'chip_tf', 'chip_histone', 'contact_maps', 'splice_sites', 'splice_site_usage', 'splice_junctions']
atac shape: torch.Size([1, 131072, 256])
```

### Notes on output format
- Output is a dict keyed by modality
- Each value is a tuple; index [1] gives the prediction tensor
- ATAC shape `[1, 131072, 256]` = 1 sequence × 131,072 positions × 256 cell-type tracks
- `organism_index=0` = human (1 = mouse)
- The MPS fallback warning is expected and harmless — `PYTORCH_ENABLE_MPS_FALLBACK=1` handles it silently

### Track metadata (from Choo, July 3)
Source: `google-deepmind/alphagenome_research` →
`/src/alphagenome_research/model/metadata/OutputMetadataResponse_ORGANISM_HOMO_SAPIENS.textproto`

| Head | Adult Postmortem | Pre/Perinatal Postmortem | Cell Line | Other/Unknown | Total Brain | Head Total |
|---|---|---|---|---|---|---|
| CHIP_HISTONE | 22 | 16 | 47 | 15 | 100 | 1116 |
| SPLICE_SITE_USAGE | 22 | 14 | 16 | 2 | 54 | 734 |
| CHIP_TF | 1 | 0 | 47 | 4 | 52 | 1617 |
| RNA_SEQ | 13 | 15 | 16 | 2 | 46 | 667 |
| CAGE | 0 | 0 | 2 | 32 | 34 | 546 |
| SPLICE_JUNCTIONS | 11 | 7 | 8 | 1 | 27 | 367 |
| DNASE | 4 | 2 | 8 | 7 | 21 | 305 |
| ATAC | 0 | 0 | 1 | 0 | 1 | 167 |
| **TOTAL** | **73** | **54** | **145** | **63** | **335** | **5519** |

**Important:** ATAC has only 1 brain track out of 167 total. DNase has 21 brain tracks and CHIP_HISTONE has 100. For brain-specific ISM analysis, DNase and CHIP_HISTONE are more informative than ATAC in the pretrained model.

### Status
Task 1 complete. AlphaGenome runs locally on Mac with MPS acceleration.


**Goal:** Task 2 — understand variant_scoring docs (mask, aggregation) before writing ISM script
 
**Source:** https://www.alphagenomedocs.com/variant_scoring.html
 
### How variant scoring works — pipeline overview
 
Variant scoring has four steps:
 
**Step 1: Make REF and ALT predictions**
Feed two sequences into AlphaGenome — one with the reference allele, one with the alternate allele — and get back full prediction tensors for both. Output shape: `[1, 131072, N_tracks]` for each. We already have the ref/alt `.npz` files ready for 5 pilot SNPs.
 
**Step 2: Apply a spatial mask**
Rather than scoring the entire 131,072 bp window, define a region of interest (the "mask") and discard everything outside it. For chromatin accessibility (ATAC/DNase), the recommended mask is a **501 bp window centered on the variant**. This makes biological sense — a SNP disrupting a TF binding site changes local chromatin accessibility right around where it sits, not across the whole window.
 
The mask collapses the `[131072]` position axis down to the 501 positions around the SNP.
 
**Step 3: Aggregate spatially and compute ALT − REF**
Reduce those 501 positions to a single scalar per track. For chromatin accessibility specifically, the recommended aggregation is:
 
```
log2[(sum(ALT) + 1) / (sum(REF) + 1)]
```
 
The +1 is a pseudocount to avoid log(0). Result: one number per track.
- Positive = alt allele opens chromatin
- Negative = alt allele closes chromatin
This is the variant effect score.
 
**Step 4: Optionally aggregate across tracks**
After one score per track, optionally average/max across subsets of tracks — e.g., mean across all 21 brain DNase tracks → single "brain DNase effect" per SNP. This is how we produce the brain vs non-brain comparison.
 
### AggregationType naming convention
Options like `DIFF_LOG2_SUM` are read right to left (innermost operation first):
- `DIFF_LOG2_SUM` → apply log2 first, then sum, then compute ALT − REF
- This is the chromatin accessibility formula above

### Key concepts 
 
**What is a mask?**
A spatial filter that restricts scoring to a local window around the SNP (501 bp for ATAC/DNase), so irrelevant signal from the rest of the 131 kb window is excluded. Can also be gene-body or exon-based for splicing/expression modalities.
 
**What is aggregation?**
The mathematical operation that collapses masked positions into a single scalar per track. For ATAC/DNase: log2-ratio of summed signals. For RNA-seq: log-fold change of mean signal over exons. The result is always one number per track, enabling comparison across many tracks.
 
**How do you measure impact on one specific track?**
After computing the variant effect score for all tracks, index into the result by track index to get the score for a single tissue/cell type. The track metadata file (`OutputMetadataResponse_ORGANISM_HOMO_SAPIENS.textproto`) maps track indices to tissue names.
 
### Important for our script
The `variant_scoring` module in `alphagenome_pytorch` implements mask creation, aggregation, and recommended scorers per modality. Our pilot script just needs to: load sequences from `.npz` → call the variant scorer → collect scores per track → subset to brain vs non-brain tracks.
 
### Status
Task 2 complete. Understand mask and aggregation; know which scorer to use for ATAC/DNase.

 
**Goal:** Task 3 — confirm scz2018clozuk GWAS coordinates are hg38 by verifying ref allele at pilot SNP positions
 
### Method
Queried Ensembl REST API for hg38 and hg19 positions of all 5 pilot SNPs and compared to GWAS file coordinates.
 
```bash
# hg38 query
for rsid in rs2007044 rs12668848 rs2660304 rs12416331 rs4129585; do
    echo -n "$rsid | "
    curl -s "https://rest.ensembl.org/variation/human/${rsid}?content-type=application/json" | python3 -c "
import json, sys
d = json.load(sys.stdin)
for m in d.get('mappings', []):
    loc = str(m.get('location',''))
    if ':' in loc and 'PATCH' not in loc:
        print(loc, '| alleles:', m.get('allele_string'), '| strand:', m.get('strand'))
"
done
# hg19: same command using grch37.rest.ensembl.org instead
```
 
### Results
 
**hg38 (Ensembl GRCh38) — none match GWAS positions:**
| rsid | GWAS pos | hg38 pos (dbSNP) | offset |
|---|---|---|---|
| rs2007044 | chr12:2,344,960 | chr12:2,235,794 | ~109 kb |
| rs12668848 | chr7:2,020,995 | chr7:1,981,360 | ~40 kb |
| rs2660304 | chr1:98,512,127 | chr1:98,046,571 | ~466 kb |
| rs12416331 | chr10:104,928,914 | chr10:103,169,157 | ~1.76 Mb |
| rs4129585 | chr8:143,312,933 | chr8:142,231,572 | ~1.08 Mb |
 
**hg19 (Ensembl GRCh37) — all match GWAS positions exactly:**
| rsid | GWAS pos | hg19 pos (dbSNP) | match? |
|---|---|---|---|
| rs2007044 | chr12:2,344,960 | chr12:2,344,960 | ✅ |
| rs12668848 | chr7:2,020,995 | chr7:2,020,995 | ✅ |
| rs2660304 | chr1:98,512,127 | chr1:98,512,127 | ✅ |
| rs12416331 | chr10:104,928,914 | chr10:104,928,914 | ✅ |
| rs4129585 | chr8:143,312,933 | chr8:143,312,933 | ✅ |
 
### Conclusion
**scz2018clozuk uses hg19 (GRCh37) coordinates.** All 5 pilot SNP positions match hg19 exactly and are off from hg38 by 40 kb – 1.76 Mb, which is the classic signature of a genome build mismatch. Alleles are consistent. This confirms Choo's suspicion that most GWAS use hg19.
 
### Action required: liftover hg19 → hg38
AlphaGenome uses hg38, so all coordinates in `scz_lead_snps.tsv` (141 SNPs) must be lifted over before running ISM. Plan:
- Use UCSC liftOver on Hoffman2 (where the TSV files live)
- Chain file needed: `hg19ToHg38.over.chain.gz`
- Input: BED-format coordinates derived from `scz_lead_snps.tsv`
- Output: hg38 coordinates, rejoined with rsID/allele columns
- Verify: spot-check a few lifted positions against the dbSNP hg38 positions confirmed above
Note: `extract_sequences.py` was written assuming hg38 input — after liftover it can be used as-is, only the input coordinates change.
 
### Status
Task 3 complete. GWAS confirmed hg19. Liftover is the immediate next coding step before ISM can run.
 
---

## Session 2 — July 23-24, 2026
**Goal:** Reorganize project file structure, set up GitHub repo, and configure VS Code
 
### GitHub setup
Created private repo and pushed initial commit:
```bash
# on Hoffman2
git init
git add .gitignore README.md phase1/get_lead_snps.py phase2/extract_and_mutate.py docs/ notebooks/
git commit -m "Initial commit: project structure, get_lead_snps.py, extract_and_mutate.py, lab notebook"
git remote add origin https://github.com/aryanasatheesh/alphagenome-project.git
git push -u origin master
```
Data files (`.gz`, `.fa`, `.safetensors`, `.npz`, `.bw`) are excluded from git via `.gitignore` — they stay on disk only and are too large for version control.

Authentication uses a GitHub Personal Access Token (classic) with `repo` scope, entered as the password when pushing.
 
### Git workflow going forward
```bash
# standard workflow after a work session:
git add -A
git commit -m "description of changes"
git push
 
# when switching machines:
git pull    # on the machine you're starting work on
# ... do work ...
git add -A && git commit -m "message" && git push
# on the other machine later:
git pull
```
 
GitHub is the bridge between Mac and Hoffman2. Code lives on GitHub; large data files live on each machine separately.
 
### VS Code setup
- Installed VS Code on Mac
- Signed in with GitHub
- Installed **Remote - SSH** extension (by Microsoft)
- Attempted to connect directly to Hoffman2 via Remote-SSH — **failed** due to:
  - `noexec` restriction on Hoffman2 home directory filesystem (VS Code server binary can't run)
  - Symlink workaround to `/u/project/cluo/aryasath/.vscode-server` partially worked but server still crashed with "channel closed" error
- **Workaround adopted:** Clone repo locally on Mac, edit in VS Code, push to GitHub, pull on Hoffman2. SSH to Hoffman2 from a VS Code terminal tab for remote work.
```bash
# on Mac
cd ~
git clone https://github.com/aryanasatheesh/alphagenome-project.git
# then File → Open Folder → ~/alphagenome-project in VS Code
```
 
### Status
GitHub repo live, VS Code configured, git sync verified (push from Mac → pull on Hoffman2 works).
 

**Goal:** Complete environment setup on Hoffman2, verify both machines ready
 
### Hoffman2 environment setup
 
**Step 1: Create conda environment**
```bash
module load anaconda3
conda create -n alphagenome_env python=3.12 -y
conda activate alphagenome_env
```
Created with Python 3.12.11.
 
**To activate in future sessions:**
```bash
module load anaconda3
conda activate alphagenome_env
```
 
**Step 2: Install alphagenome-pytorch**
Direct `pip install alphagenome-pytorch` failed because numpy couldn't compile from source (Hoffman2's GCC 4.8.5 is too old for the meson build system). Fix:
```bash
conda install numpy -y                      # prebuilt binary, avoids GCC issue
pip install alphagenome-pytorch --no-deps    # install without re-pulling numpy
pip install torch safetensors einx           # install remaining dependencies
```
Installed: alphagenome-pytorch==0.3.0, torch==2.6.0+cu124, numpy==2.0.1, einx==0.4.3, safetensors==0.8.0
 
**Step 3: Download model weights**
```bash
mkdir -p ~/alphagenome-weights
cd ~/alphagenome-weights
pip install huggingface_hub
hf download gtca/alphagenome_pytorch model_fold_0.safetensors --local-dir .
# downloaded 921 MB to ~/alphagenome-weights/model_fold_0.safetensors
```
 
**Step 4: Verification**
```bash
python3 -c "
import torch
from alphagenome_pytorch import AlphaGenome
print('torch:', torch.__version__)
print('cuda available:', torch.cuda.is_available())
print('alphagenome imported OK')
"
```
Output:
```
torch: 2.6.0+cu124
cuda available: False    # expected — login nodes have no GPU
alphagenome imported OK
```
 
Model load test on login node failed with OOM (`can't allocate memory`) — this is expected because login nodes have limited RAM. The model (~920 MB parameters) needs a compute node. Imports and weight download are verified; full forward pass will be confirmed when we submit the first GPU job.
 
### Environment summary — both machines
 
| Component | Mac (local dev) | Hoffman2 (batch compute) |
|---|---|---|
| Python | 3.12.13 | 3.12.11 |
| torch | 2.13.0 | 2.6.0+cu124 |
| alphagenome-pytorch | 0.3.1 | 0.3.0 |
| numpy | 2.5.1 | 2.0.1 |
| GPU type | MPS (Apple Silicon) | CUDA (on compute nodes) |
| Model weights | ~/alphagenome-weights/ | ~/alphagenome-weights/ |
| Activate env | `source ~/alphagenome-env/bin/activate` | `module load anaconda3 && conda activate alphagenome_env` |
| Run scripts | `PYTORCH_ENABLE_MPS_FALLBACK=1 python3 script.py` | submit via SLURM/qsub |
 
### How switching between machines works
- **Code** (scripts, notebooks, docs): lives on GitHub. Push from one machine, pull on the other.
- **Large data** (GWAS files, hg38.fa, model weights, bigwigs): stays on each machine separately. Too large for git.
- **Rule of thumb:** develop and test on Mac with small inputs (1–5 SNPs), push to GitHub, pull on Hoffman2, run at scale (141+ SNPs, batch GPU jobs).
### Status
Task 2 complete. Both environments verified and ready.
 
**Next:** Task 4 — Download PGC3 SCZ GWAS, rerun data acquisition pipeline, document liftover logic thoroughly.

---

## Session 3 - August 3, 2026
**Goal:** Download the latest PGC3 SCZ GWAS, build a robust lead SNP pipeline with documented MHC analysis, liftover to hg38, and validate
 
### Context and motivation for switching to PGC3
In our July 22 meeting, Choo recommended switching from the 2018 CLOZUK+PGC2 GWAS (scz2018clozuk) to the latest available SCZ GWAS from PGC. The reasoning is straightforward: larger sample sizes yield more statistical power, which means more discovered loci and more precise effect size estimates. Since we are using GWAS hits as inputs to AlphaGenome (not re-analyzing the GWAS itself), using the most comprehensive available set of risk loci gives us the richest set of variants to characterize computationally.
 
### PGC3 SCZ GWAS (Trubetskoy et al., 2022, Nature)
- **Citation:** Trubetskoy et al., "Mapping genomic loci implicates genes and synaptic biology in schizophrenia," Nature, 2022. PubMed: 35396580
- **Download:** https://doi.org/10.6084/m9.figshare.19426775
- **Sample:** European ancestry subset: ~53,386 cases, ~77,258 controls (effective N ≈ 58,749)
- **Multi-ancestry (primary):** ~69,369 cases, ~236,642 controls across European, East Asian, African American, and Latino ancestries; identified 270 distinct risk loci (headline result)
### Choosing the European ancestry subset
The figshare download contains multiple files for different ancestry groups and analysis configurations. We downloaded two:
- `PGC3_SCZ_european.vcf.tsv.gz` (229 MB) — European ancestry only
- `PGC3_SCZ_primary.vcf.tsv.gz` (226 MB) — multi-ancestry primary analysis
We chose **European ancestry** for the initial analysis. The reasoning requires nuance:
 
For ISM with AlphaGenome specifically, ancestry does not affect the prediction itself — AlphaGenome is purely sequence-based, takes hg38 DNA as input, and does not know or care about population genetics. The reference genome is the same for everyone, and the ISM result ("does this SNP change predicted chromatin?") is ancestry-agnostic.
 
However, ancestry matters for **which SNPs we select as inputs**, because:
1. **Allele frequencies differ across populations.** A variant common enough to detect in Europeans may be rare in East Asians, and vice versa. The GWAS p-value depends on allele frequency in the study population.
2. **Linkage disequilibrium (LD) patterns differ.** When we apply a 1 Mb clumping window to select independent lead SNPs, the result depends on the LD structure of the study population. In populations with longer LD blocks, more SNPs get clumped together.
3. **The lead SNP from clumping may not be the causal SNP.** It is the most statistically significant SNP within a correlated block, and which SNP that is depends on LD structure.
4. **Downstream tools are better calibrated for European populations.** LD reference panels, fine-mapping methods, and most existing functional annotations are most complete for European ancestry data.
5. **Comparability with prior work.** The 2018 CLOZUK+PGC2 GWAS was also European ancestry.
Either European or primary (multi-ancestry) could be justified. We use European for now and have the primary file available if we want to expand later. This is a decision worth revisiting with Choo — the multi-ancestry analysis has more power (270 loci vs. ~230 in European-only) but introduces complexity in LD-based SNP selection.
 
### File format inspection
The PGC3 files use a VCF-style tab-separated format (`PGCsumstatsVCFv1.0`) with `##` comment headers. Key columns:
 
| Column | Meaning | Notes |
|---|---|---|
| CHROM | Chromosome | No 'chr' prefix |
| ID | rsID | |
| POS | Position | **hg19 (GRCh37)** — confirmed below |
| A1 | Effect allele | This is the **alt** allele |
| A2 | Other allele | This is the **ref** allele |
| BETA | Log odds ratio | Different from 2018 file which used OR |
| SE | Standard error | |
| PVAL | P-value | Column name differs from 2018 (was "P") |
| NCAS | N cases | |
| NCON | N controls | |
| NEFF | Effective N | |
 
Compared to the 2018 CLOZUK file: column names differ (CHROM vs CHR, ID vs SNP, PVAL vs P, BETA vs OR) but the data structure is equivalent. A1 is still the effect/alt allele, A2 is still the ref allele.
 
### Coordinate system verification
Applied the same sanity check method from Session 3: looked up rs2007044 in the European file.
```bash
zcat PGC3_SCZ_european.vcf.tsv.gz | grep -v "^##" | grep "rs2007044"
# Output: chr12  rs2007044  2344960  A  G  ...
```
POS = 2,344,960, which we previously confirmed is the hg19 position for this SNP (hg38 would be 2,235,794). **Conclusion: PGC3 also uses hg19 coordinates.** This is consistent across PGC releases.
 
### Summary statistics overview
 
| Metric | European | Primary |
|---|---|---|
| Total SNPs | 7,659,767 | 7,585,077 |
| Genome-wide significant (p < 5e-8) | 20,457 | 21,723 |
 
Compared to 2018 CLOZUK: 8,064,800 total SNPs, 18,088 GW-significant. PGC3 has slightly fewer total SNPs (different imputation panel) but more genome-wide significant hits (larger sample = more power).
 
### MHC dominance analysis — why naive "top N SNPs" fails
 
Choo asked us to first pull the top 10–50 SNPs by p-value and note where they fall. This is a pedagogically important exercise: it demonstrates a fundamental challenge in schizophrenia genetics.
 
**Result: all 50 of the top 50 SNPs by p-value are on chromosome 6.**
 
They span positions ~27.5–28.8 Mb (hg19), which falls squarely within the Major Histocompatibility Complex (MHC) region (~25–34 Mb on chr6). The first non-chr6 SNP does not appear until **rank #1,169** (rs58120505, chr7:2,029,867, p = 2.235e-24).
 
**Why does this happen?** The MHC region is one of the most gene-dense and polymorphic regions in the human genome. It encodes proteins critical for immune function (HLA genes). Several properties make it dominate GWAS results:
1. **Extreme genetic diversity.** The MHC harbors more common variants than almost any other region, giving GWAS more "chances" to detect associations.
2. **Complex and extended LD structure.** LD blocks in the MHC can stretch megabases — much longer than typical genomic LD (~100–200 kb in Europeans). This means thousands of SNPs are correlated with each other, all reflecting a small number of underlying causal signals.
3. **Strong biological effect.** The MHC genuinely harbors schizophrenia risk variants (complement component C4 is a well-characterized example), but the signal is amplified by the region's unusual genetic architecture.
4. **Difficult to interpret mechanistically.** The extreme LD makes it nearly impossible to fine-map causal variants in the MHC using standard methods. For ISM, running AlphaGenome on thousands of correlated MHC SNPs would be computationally wasteful and scientifically uninformative.
**Quantitative breakdown:**
- GW-sig SNVs on chr6: 6,104 / 20,457 (29.8%)
- GW-sig SNVs in the MHC region specifically (chr6:25–34 Mb): 5,994
- GW-sig SNVs outside chr6: 14,353
So nearly 30% of all genome-wide significant SNPs come from one region that spans ~9 Mb (0.3% of the genome). This is why we exclude the MHC from lead SNP selection — not because it's unimportant biologically, but because its unusual genetic structure makes it unsuitable for the per-variant ISM approach.
 
### Lead SNP pipeline (get_lead_snps_pgc3.py)
 
**Memory optimization.** The initial version of the script loaded all 7.6M SNPs into memory before filtering, which caused a MemoryError on Hoffman2 login nodes (limited RAM). Fixed by filtering on the fly: only genome-wide significant SNVs are retained during reading, reducing memory usage from ~7.6M rows to ~20k.
 
**Pipeline steps:**
1. **Read and filter to GW-significance (p < 5e-8):** Stream through the gzipped file, parsing each line and keeping only SNPs with PVAL < 5e-8. This reduces 7,659,767 SNPs to 20,457.
2. **Remove indels:** AlphaGenome accepts only single-nucleotide substitutions as input (it predicts the effect of swapping one base for another at a specific position). Multi-base variants (insertions, deletions) cannot be represented in this framework. In practice, 0 indels were found among the GW-significant SNPs in this dataset (likely because the PGC3 VCF format already restricts to biallelic SNPs).
3. **Document MHC dominance:** Save the top 50 SNPs to a separate file (`scz_pgc3_mhc_analysis.tsv`) for reference, and report the MHC statistics described above.
4. **Exclude MHC region:** Remove all SNPs in chr6:25,000,000–34,000,000 (hg19). This is a deliberately conservative window that extends beyond the classical MHC boundaries to capture the full extent of the extended MHC LD. After exclusion: 14,463 SNVs remain.
5. **1 Mb clumping:** Apply greedy distance-based clumping to define independent loci. The algorithm sorts SNPs by p-value (most significant first), selects the top SNP as a lead, removes all SNPs within 1 Mb on the same chromosome, and repeats. The 1 Mb window is standard in GWAS (used by PGC, PLINK, and most post-GWAS tools) because it roughly captures the extent of LD in European populations — SNPs more than 1 Mb apart are unlikely to be in strong LD and can be considered independent signals. After clumping: **173 independent lead SNPs**.
6. **Select pilot SNPs:** Take the top 10 lead SNPs by p-value for initial pipeline development and testing.
**Output:**
- `scz_pgc3_lead_snps.tsv` — 173 independent lead SNPs (hg19)
- `scz_pgc3_pilot_10snps.tsv` — 10 pilot SNPs (hg19)
- `scz_pgc3_mhc_analysis.tsv` — top 50 SNPs documenting MHC dominance
### Pilot SNPs (hg19)
 
| rsid | chr | pos | ref | alt | pval |
|---|---|---|---|---|---|
| rs58120505 | 7 | 2,029,867 | C | T | 2.235e-24 |
| rs2238057 | 12 | 2,384,005 | G | T | 8.502e-22 |
| rs1198588 | 1 | 98,552,832 | T | A | 1.731e-21 |
| rs4702 | 15 | 91,426,560 | A | G | 2.794e-21 |
| rs13107325 | 4 | 103,188,709 | T | C | 2.900e-21 |
| rs2710323 | 3 | 52,815,905 | C | T | 1.229e-19 |
| rs12129573 | 1 | 73,768,366 | A | C | 2.282e-18 |
| rs4129585 | 8 | 143,312,933 | C | A | 5.109e-18 |
| rs778371 | 2 | 233,743,109 | G | A | 1.495e-17 |
| rs11191580 | 10 | 104,906,211 | C | T | 1.772e-17 |
 
Note: rs4129585 (chr8) was also in our original 2018 pilot set — provides cross-study continuity.
 
### Liftover hg19 → hg38
 
Applied the same liftover procedure as Session 4 using pyliftover with the UCSC hg19ToHg38.over.chain.gz chain file.
 
**How liftover works (conceptual):** The human reference genome has been assembled multiple times as sequencing technology improved. hg19 (GRCh37, released 2009) and hg38 (GRCh38, released 2013) differ because hg38 incorporated better assemblies of centromeric regions, added alternate loci, fixed sequencing errors, and closed gaps. These changes shift the coordinates of most genomic positions — a gene that starts at position X in hg19 may start at position X+Δ in hg38, where Δ varies by region (from a few bp to several megabases). The UCSC "chain file" encodes these coordinate mappings: it specifies, for each contiguous block in hg19, where that block maps in hg38. Pyliftover reads this chain file and performs the lookup.
 
**Coordinate convention detail:** GWAS summary statistics and most genomics text files use 1-based coordinates (position 1 = first base). Pyliftover internally uses 0-based half-open intervals (consistent with BED format and Python indexing). Our script converts: subtract 1 before calling pyliftover, add 1 to the result.
 
**Results:**
```
scz_pgc3_lead_snps.tsv: 173 lifted, 0 failed
scz_pgc3_pilot_10snps.tsv: 10 lifted, 0 failed
```
 
### Liftover validation
 
**Method 1: Cross-reference with previous validation.** rs4129585 appears in both the 2018 and PGC3 pilot sets. Its hg38 position from this liftover (chr8:142,231,572) matches exactly what we validated against dbSNP in Session 3.
 
**Method 2: Spot-check new SNP against Ensembl.** Queried rs58120505 (the #1 non-MHC SNP, not in our previous set):
```bash
curl -s "https://rest.ensembl.org/variation/human/rs58120505?content-type=application/json"
# Output: 7:1990232-1990232 | alleles: T/A/C/G
```
Our liftover gives chr7:1,990,232 — exact match. Alleles T/A/C/G include our ref=C and alt=T.
 
### Pilot SNPs (hg38, final)
 
| rsid | chr | pos (hg38) | ref | alt | pval |
|---|---|---|---|---|---|
| rs58120505 | 7 | 1,990,232 | C | T | 2.235e-24 |
| rs2238057 | 12 | 2,274,839 | G | T | 8.502e-22 |
| rs1198588 | 1 | 98,087,276 | T | A | 1.731e-21 |
| rs4702 | 15 | 90,883,330 | A | G | 2.794e-21 |
| rs13107325 | 4 | 102,267,552 | T | C | 2.900e-21 |
| rs2710323 | 3 | 52,781,889 | C | T | 1.229e-19 |
| rs12129573 | 1 | 73,302,683 | A | C | 2.282e-18 |
| rs4129585 | 8 | 142,231,572 | C | A | 5.109e-18 |
| rs778371 | 2 | 232,878,399 | G | A | 1.495e-17 |
| rs11191580 | 10 | 103,146,454 | C | T | 1.772e-17 |
 
### Problems encountered and solutions
 
1. **figshare download URL:** The initial wget to `figshare.com/ndownloader/articles/19426775/versions/2` returned a 0-byte file (HTTP 202 Accepted without data). Fixed by querying the figshare API for individual file download URLs and using those directly.
2. **Memory error on login node:** The initial `get_lead_snps_pgc3.py` loaded all 7.6M SNPs into a list before filtering, exceeding login node RAM. Fixed by rewriting to filter on the fly during reading — only ~20k significant SNPs are kept in memory.
3. **Coordinate system (hg19 vs hg38):** Already known from Session 3 that PGC GWAS uses hg19, but verified independently for PGC3 using the same rs2007044 check. This confirms the pattern holds across PGC releases.
### Output files on Hoffman2
All in `/u/project/cluo/aryasath/alphagenome_project/phase1/data/`:
- `PGC3_SCZ_european.vcf.tsv.gz` — raw GWAS summary statistics (229 MB)
- `PGC3_SCZ_primary.vcf.tsv.gz` — multi-ancestry version (226 MB, downloaded but not used yet)
- `scz_pgc3_mhc_analysis.tsv` — top 50 SNPs (MHC documentation)
- `scz_pgc3_lead_snps.tsv` — 173 lead SNPs, hg19
- `scz_pgc3_lead_snps_hg38.tsv` — 173 lead SNPs, hg38
- `scz_pgc3_pilot_10snps.tsv` — 10 pilot SNPs, hg19
- `scz_pgc3_pilot_10snps_hg38.tsv` — 10 pilot SNPs, hg38
### Status
✅ Task 4 complete. PGC3 GWAS downloaded, processed, MHC dominance documented, 173 lead SNPs identified, 10 pilot SNPs selected, all lifted to hg38 and validated.
 
**Immediate next steps:**
1. Write ISM inference script using variant_scoring module
2. Run ISM on 10 pilot SNPs (locally on Mac first, then Hoffman2 for batch)
3. Visualize delta ATAC tracks, compare brain vs non-brain

---

## Session 4 — August 4, 2026
**Goal:** Write and run the ISM inference script on the 10 PGC3 pilot SNPs
 
### Approach
Rather than pre-extracting sequences to `.npz` files (the spring approach), the ISM script fetches sequence windows on demand from the Ensembl REST API. This avoids needing the 3 GB hg38 FASTA on the Mac and keeps the pipeline self-contained for small pilot runs. For batch runs on Hoffman2 we will switch to reading from the local `hg38.fa` to avoid API rate limits.
 
### Pipeline steps implemented in `run_ism_pilot.py`
1. Fetch a 131,072 bp hg38 window centered on the variant
2. Verify which GWAS allele matches the hg38 reference base; assign true ref/alt
3. One-hot encode ref and alt sequences (A=0, C=1, G=2, T=3)
4. Run both through the pretrained AlphaGenome model
5. Apply a 501 bp spatial mask centered on the variant
6. Compute per-track variant effect score: `log2[(sum(alt)+1) / (sum(ref)+1)]`
7. Save one row per SNP × modality × track
### Key finding: GWAS A1/A2 are not genome ref/alt
Before running ISM, we spot-checked rs58120505: the hg38 base at chr7:1,990,232 is **T**, but our pipeline had labeled `ref=C, alt=T`. Checking all 10 pilot SNPs showed the same pattern — **in all 10 cases the hg38 reference base matched the GWAS A1 ("effect allele"), not A2.**
 
This is not a liftover error. In GWAS summary statistics, A1 and A2 denote *statistical roles*:
- **A1 = effect allele** — the allele whose effect on risk is measured by BETA
- **A2 = other allele** — the non-effect allele
Neither is defined relative to the reference assembly. Either can be the reference genome allele. Our `get_lead_snps_pgc3.py` had assumed A2 = ref, which is not reliable.
 
For ISM this distinction is critical: the "ref" sequence must be the actual hg38 sequence (otherwise we would be mutating a base that isn't there), and the "alt" sequence must contain the other allele. Feeding the sequences in the wrong orientation would flip the sign of every variant effect score.
 
**Fix:** the ISM script now checks the hg38 base at the variant position and assigns:
- `true_ref` = whichever of A1/A2 matches hg38
- `true_alt` = the other one
- Records a `gwas_alleles_swapped` flag per SNP for transparency
- Errors out if hg38 matches *neither* allele (would indicate a genuine coordinate problem)
This does **not** invalidate the lead SNP selection, clumping, or liftover — none of those depend on allele orientation.
 
### Run details
- Ran on Mac (MPS), all 10 pilot SNPs
- Two forward passes per SNP (ref and alt), 131,072 bp each
- Modalities scored: ATAC and DNase
- Output: 6,400 rows in `results/ism/pilot_ism_scores.tsv`
### Results summary (all channels, pre-padding-filter — see Session 9)
 
| SNP | chr | ATAC max | ATAC mean | DNase max | DNase mean | DNase tracks>0.1 |
|---|---|---|---|---|---|---|
| rs2710323 | 3 | 5.19 | 0.53 | 4.42 | 0.41 | 243 |
| rs2238057 | 12 | 0.32 | 0.02 | 0.37 | 0.03 | 34 |
| rs1198588 | 1 | 0.15 | 0.02 | 0.70 | 0.06 | 45 |
| rs4702 | 15 | 0.17 | 0.02 | 0.24 | 0.04 | 59 |
| rs13107325 | 4 | 0.04 | 0.01 | 0.21 | 0.03 | 60 |
| rs4129585 | 8 | 0.34 | 0.03 | 0.42 | 0.06 | 88 |
| rs11191580 | 10 | 0.14 | 0.03 | 0.19 | 0.06 | 70 |
| rs12129573 | 1 | 0.06 | 0.02 | 0.27 | 0.04 | 48 |
| rs778371 | 2 | 0.08 | 0.02 | 0.12 | 0.01 | 2 |
| rs58120505 | 7 | 0.04 | 0.00 | 0.05 | 0.01 | 0 |
 
**Note:** these summary statistics were computed across all 256 ATAC / 384 DNase output channels. Session 9 established that a substantial fraction of those channels are zero-padding, so these means are diluted and must be recomputed on real tracks only. The raw per-track scores are unaffected.
 
### Preliminary observations
- **rs2710323** is a large outlier, with effects across most channels. A variant that changes accessibility in nearly every tissue is more consistent with a ubiquitous mechanism than a brain-specific one — per Choo's earlier guidance, this pattern argues *against* SCZ specificity.
- **rs58120505** has the strongest non-MHC GWAS p-value but essentially no predicted chromatin effect. GWAS significance and predicted regulatory impact are clearly not the same thing. Possible explanations: the causal mechanism is not chromatin-mediated; or the lead SNP is not the causal variant but is in LD with it.
- Several SNPs (rs1198588, rs4129585) show larger DNase than ATAC effects, which is consistent with DNase having far broader tissue coverage in this model.
### Status
✅ First real ISM results. Pipeline works end to end.

---

## Session 5 — August 21–25, 2026
**Goal:** Map output channels to tissues so we can compute brain vs non-brain scores (Step 1 of post-meeting plan)
 
### Context — feedback from Choo (Aug 21 meeting)
Main points from the meeting:
 
1. **Sequence models sidestep the LD problem.** In real data, LD means a whole correlated block lights up and it is very hard to say which SNP is causal. Fine mapping tries to solve this statistically. AlphaGenome has no notion of LD or inheritance — it evaluates each variant independently from sequence alone, so we can perturb one SNP without its neighbors. This is a genuine argument for sequence-to-function models as a complement to fine mapping, and belongs in the thesis introduction.
2. **The European-ancestry constraint is less relevant for our purposes** than it would be for an empirical study, for the same reason — we are not relying on LD structure for the ISM itself.
3. **We need negative controls to interpret score magnitude.** Without a baseline we cannot say whether 0.34 is large or small. Two options: an external control (a non-brain trait GWAS such as height, via the GWAS Catalog), and an internal control (non-significant SNPs drawn from the same GWAS, e.g. lowest p-value within a 1 Mb region).
4. **Brain-specific scoring is itself a control.** We would not expect SCZ variants to affect liver as much as brain. Scoring brain tracks separately, rather than averaging over all tissues, is both the biological question and a sanity check.
5. Clumping could use predefined LD/recombination regions instead of a fixed ±1 Mb window.
6. Look up W.E.I.R.D. science (bias in which populations are studied).
Agreed plan: (1) brain vs non-brain tracks, (2) negative controls, (3) write up the fine-mapping argument.
 
### Obtaining track metadata
Downloaded the official metadata from the `google-deepmind/alphagenome_research` repo:
```bash
curl -sL -o human_output_metadata.textproto \
  "https://raw.githubusercontent.com/google-deepmind/alphagenome_research/main/src/alphagenome_research/model/metadata/OutputMetadataResponse_ORGANISM_HOMO_SAPIENS.textproto"
```
2.5 MB textproto. Each track record contains `name`, `strand`, `ontology_term`, `biosample {type, name, stage}`, `assay`, `data_source`, and `nonzero_mean`. Records are grouped under `output_metadata { output_type: ... }` blocks.
 
Parsed into `phase2/metadata/track_metadata.tsv` with one row per track.
 
### Track counts per output type
| output_type | rows |
|---|---|
| ATAC | 256 |
| DNASE | 384 |
| CAGE | 640 |
| RNA_SEQ | 768 |
| CHIP_HISTONE | 1152 |
| CHIP_TF | 1664 |
| PROCAP | 128 |
| SPLICE_SITES | 5 |
| SPLICE_SITE_USAGE | 734 |
| SPLICE_JUNCTIONS | 367 |
| CONTACT_MAPS | 28 |
| **TOTAL** | **6126** |
 
These match the model's output tensor widths exactly (ATAC 256, DNase 384), which was the first reassurance that the metadata corresponds to this checkpoint.
 
### Reconciling with Choo's track table
Choo's table listed ATAC = 167 and DNase = 305 total tracks, which did not match the 256/384 output width. Initially unclear whether this indicated a version mismatch. Counting *unique biosample names* gave ATAC 168 and DNase 306 — one more than Choo's numbers in each case. That "+1" turned out to be the empty string shared by padding rows (see below). Both numbers were correct; they were counting different things.
 
### False positives in naive keyword matching
A first pass classified tracks as brain-related by substring matching. This produced clear errors:
- "g**astro**cnemius medialis", "g**astro**esophageal sphincter" → matched "astro"
- "renal **cortex** interstitium", "renal **cortic**al epithelial cell" → matched "cortex"/"cortic"
- "sciatic nerve", "tibial nerve" → peripheral nervous system, not CNS
**Lesson:** substring matching on biosample names is unsafe. Replaced with a manually curated exact-match list built by reading the actual names present in the file. Recorded in `track_metadata_classified.tsv` with a `tissue_class` column (CNS / eye / other).
 
Eye tissue (`retina`, `retinal pigment epithelial cell`) was classified separately because Choo suggested eye as a candidate negative-control tissue.
 
### Tissue classification results
| Modality | CNS | eye | other | total |
|---|---|---|---|---|
| ATAC | **1** | 0 | 255 | 256 |
| DNASE | 23 | 2 | 359 | 384 |
| CHIP_HISTONE | 70 | 1 | 1081 | 1152 |
| RNA_SEQ | 35 | 0 | 733 | 768 |
 
**Major implication: ATAC cannot support brain-specificity claims in the pretrained model.** One CNS track (motor neuron) out of 167 real tracks. DNase (23 CNS) and especially CHIP_HISTONE (70 CNS) are where brain signal lives. This is also a strong motivation for the eventual fine-tuning on Ziffra fetal-brain ATAC — it would add the brain ATAC coverage the pretrained model lacks.
 
DNase CNS track indices: 0, 3, 22, 69, 70, 71, 72, 75, 77, 91, 92, 148, 153, 174, 201, 221, 246, 247, 251, 263, 268, 274, 292.
 
### Verifying that metadata row *i* = model channel *i*
Matching *counts* does not prove matching *order*. If the PyTorch port permuted output channels, every brain label assigned by index would be silently wrong — and wrong in a way that still produces plausible-looking numbers. We tested this before using the labels.
 
**Test 1 — correlate predicted signal with `nonzero_mean` (failed to resolve).**
Compared per-channel mean prediction against the metadata's `nonzero_mean`. Results: ATAC r = −0.437, DNase r = 0.018 (shuffled control r = 0.016). Inconclusive at best. Two candidate explanations: the ordering is wrong, or the proxy is bad (`nonzero_mean` is a genome-wide training statistic; we compared it to mean predicted signal in a single 131 kb window, likely on a different scale).
 
**Test 2 — do same-tissue tracks correlate with each other? (ambiguous).**
If ordering is correct, tracks labeled as the same tissue should have more similar prediction profiles than random tracks. Using 128 bp binned DNase output:
- 4 astrocyte tracks (indices 3, 70, 71, 72): mean pairwise r = **0.849** — high, and notable because those indices are not adjacent
- 3 SK-N-* neuroblastoma lines (148, 153, 174): r = 0.514 — at or below the null
- All 23 CNS tracks: r = 0.566 vs null 0.535, z = 1.06, p = 0.126 — not significant
The astrocyte result argued for correct ordering; the SK-N-* result argued against. Also noticed 79 DNase tracks had zero variance at this locus, which broke the first null computation (NaN) and required nan-safe handling.
 
**Test 3 — inspect the library's own code (decisive).**
Reading `alphagenome_pytorch/named_outputs.py` revealed a `_placeholder_tracks()` helper that creates tracks named `"Padding"`, plus a public `NamedOutputs.strip_padding()` method. This pointed directly at the answer: **the output heads are zero-padded to round dimensions.**
 
Checking the parsed metadata confirmed it:
 
| Head | Real tracks | Padding | Real index range | Padding index range |
|---|---|---|---|---|
| ATAC | **167** | 89 | 0–166 | 167–255 (contiguous) |
| DNASE | **305** | 79 | 0–304 | 305–383 (contiguous) |
| CHIP_HISTONE | **1116** | 36 | 0–1115 | 1116–1151 (contiguous) |
 
Padding rows are literally `name: "Padding"` with empty `assay` and empty `biosample`, and they sit contiguously at the **end** of each head. The real-track counts (167, 305, 1116) are exactly Choo's table.
 
### Conclusion: the index mapping was correct all along
Real tracks occupy indices 0…N−1 in metadata order, so metadata row *i* is model channel *i* for all real tracks. The earlier ambiguous tests were **diluted by padding channels**:
- 384 − 305 = 79 padding tracks, matching exactly the 79 zero-variance DNase channels observed in Test 2
- The random null in Test 2 was drawn partly from dead padding channels, depressing the comparison
- The astrocyte r = 0.849 was genuine signal
### Consequences for existing results
1. **The 6,400 ISM scores are valid and do not need rerunning** — they were computed at 1 bp resolution on human predictions (see note below). They do need **filtering**: drop `track_index >= 167` for ATAC and `>= 305` for DNase before any summary statistic. All max/mean values in the Session 8 table were computed over full width and are therefore diluted.
2. **Any future per-track analysis must filter padding first.** This is easy to miss and would silently corrupt results.
### Secondary finding: model output structure
While debugging, discovered that model outputs are **dicts keyed by resolution in base pairs**, not tuples:
- `out['atac'][1]` → shape `[1, 131072, 256]` — 1 bp resolution
- `out['atac'][128]` → shape `[1, 1024, 256]` — 128 bp binned (131072/128 = 1024)
- `out['chip_histone']` has **only** key `128` — this head is not available at 1 bp resolution, which is why indexing it with `[1]` raised `KeyError: 1`
- `out['contact_maps']` is a bare tensor `[1, 64, 64, 28]`, not a dict
Our ISM script used `[1]`, i.e. 1 bp resolution, which was the intended behaviour. To add CHIP_HISTONE we must use the 128 bp output and adjust the spatial mask accordingly (501 bp ≈ 4 bins at 128 bp resolution).
 
### Library metadata API (for future use)
`alphagenome_pytorch.named_outputs` provides `TrackMetadataCatalog` and `NamedOutputs`, which support filtering tracks by metadata and stripping padding automatically:
- `TrackMetadataCatalog.from_file(path)` accepts `.csv`/`.tsv`/`.parquet`; expects an `output_type` column and `track_name` or `name`; re-indexes by sorted `track_index`
- `TrackMetadataCatalog.load_builtin()` is **not** usable — no metadata ships with the package (`data/` contains only an empty `__init__.py`); it raises `FileNotFoundError` and points to a `scripts/extract_track_metadata.py` we do not have
- `NamedOutputs.from_raw(...)` and `.strip_padding()` / `.select(...)` would let us do padding removal and tissue filtering through the library rather than by hand
Worth converting `track_metadata_classified.tsv` into the catalog's expected schema so we use the library's own code path instead of manual index arithmetic.
 
### Status
✅ Track-to-tissue mapping established and verified. Padding identified and characterised.
⬜ Brain vs non-brain scoring — next session, now a short task.
 
### Files added
- `phase2/metadata/human_output_metadata.textproto` — official metadata (2.5 MB)
- `phase2/metadata/track_metadata.tsv` — parsed, one row per track
- `phase2/metadata/track_metadata_classified.tsv` — adds `tissue_class` (CNS / eye / other)
### Next steps
1. Filter padding from `pilot_ism_scores.tsv`; recompute summary statistics on real tracks only
2. Compute per-SNP CNS vs non-CNS scores for DNase; use eye tracks as a tissue control
3. Add CHIP_HISTONE to the ISM run (70 CNS tracks — best brain coverage available)
4. Build negative controls: height GWAS (external) and matched non-significant SNPs (internal)
5. Write up the fine-mapping / LD argument for the thesis introduction
---

## Session 6 — August 25, 2026
**Goal:** Filter padding from ISM scores and compute CNS vs non-CNS variant effects (Step 1 complete)
 
### Padding filter
Dropped `track_index >= 167` (ATAC) and `>= 305` (DNase) from `pilot_ism_scores.tsv`:
- rows read: 6,400
- padding dropped: 1,680
- kept: 4,720
### Scoring approach
Computed three per-SNP summaries per modality, because they answer different questions:
- **CNS mean |score|** — how broadly the variant perturbs brain tissue
- **CNS max |score|** — strongest single brain effect (outlier-sensitive)
- **CNS mean − non-CNS mean** — tissue *specificity*, the quantity relevant to Choo's question
Eye tracks (retina, retinal pigment epithelial cell) reported separately as a candidate control tissue.
 
### DNase results (305 real tracks, 23 CNS), sorted by specificity
 
| SNP | CNS mean | CNS max | other mean | CNS−other | eye mean |
|---|---|---|---|---|---|
| rs1198588 | 0.1600 | 0.5061 | 0.0587 | **+0.1013** | 0.5163 |
| rs4702 | 0.0596 | 0.1315 | 0.0476 | +0.0119 | 0.0815 |
| rs2238057 | 0.0534 | 0.2278 | 0.0427 | +0.0107 | 0.0317 |
| rs11191580 | 0.0817 | 0.1511 | 0.0718 | +0.0099 | 0.0942 |
| rs58120505 | 0.0161 | 0.0331 | 0.0161 | −0.0001 | 0.0279 |
| rs778371 | 0.0142 | 0.0419 | 0.0187 | −0.0045 | 0.0382 |
| rs12129573 | 0.0377 | 0.1908 | 0.0449 | −0.0073 | 0.0457 |
| rs13107325 | 0.0325 | 0.1560 | 0.0445 | −0.0120 | 0.0059 |
| rs4129585 | 0.0623 | 0.1926 | 0.0785 | −0.0162 | 0.0885 |
| rs2710323 | 0.1701 | 0.5702 | 0.5470 | **−0.3769** | 0.2764 |
 
### ATAC results are not interpretable
With only 1 CNS track (motor neuron) out of 167 real tracks, "CNS mean" and "CNS max" are the same single number and have no spread. The ATAC panel is reported in the notebook for completeness but should **not** be presented as evidence of tissue specificity.
 
### Per-track breakdown — which CNS tracks respond
To check whether CNS means are broad effects or single-track artifacts, examined the top CNS tracks for the two extreme SNPs.
 
**rs1198588** (all-track mean |score| = 0.0693, 95th pct = 0.3071):
```
0.5061  putamen
0.4574  neural progenitor cell
0.4229  head of caudate nucleus
0.3324  dorsolateral prefrontal cortex
0.3251  neuronal stem cell
0.2577  ecto neural progenitor cell
0.2219  SK-N-MC
0.1923  spinal cord
```
Broad and **neuronally coherent**: striatum (putamen, caudate), cortex (DLPFC), and neural progenitor/stem populations. Six of the top eight are neurons or neuronal precursors. The top three exceed this variant's own 95th percentile across all tracks, so this is not one outlier track carrying an average. Striatum and DLPFC are both well-established schizophrenia-relevant regions.
 
**rs2710323** (all-track mean |score| = 0.5168, 95th pct = 1.8744):
```
0.5702  brain microvascular endothelial cell
0.4782  astrocyte of the spinal cord
0.4614  astrocyte of the hippocampus
0.4393  smooth muscle cell of the brain vasculature
0.3512  brain pericyte
0.2977  astrocyte of the cerebellum
0.2715  astrocyte
0.1574  neural progenitor cell
```
Top CNS hits are **vascular and glial**, not neuronal — endothelium, pericytes, smooth muscle, astrocytes. Only one neuronal progenitor track appears, at 0.157. Every one of these scores sits far below this variant's own typical track (mean 0.517, 95th pct 1.874), so the brain response is a weak subset of a genome-wide-scale perturbation.
 
### Interpretation (preliminary)
The two extremes separate on cell-type composition, not just magnitude:
- **rs1198588** — focused neuronal effect standing out against its own low background. Best candidate for a brain-specific regulatory variant in this pilot set.
- **rs2710323** — broadly active variant where brain tracks respond *less* than average, and the brain tracks that do respond are vascular/glial rather than neuronal. Consistent with Choo's point that a variant affecting all tissues is likely not disease-specific.
### Caveats — what this does not establish
1. **No null distribution.** CNS−other values span −0.38 to +0.10, with most within ±0.02 of zero. Without negative controls we cannot say whether +0.10 is notable or ordinary. This is the gap the controls are designed to fill.
2. **Tracks are not independent.** Accessibility profiles are correlated across related tissues, so eight neuronal tracks agreeing is not eight independent observations. Effective sample size is smaller than 23.
3. **Eye is not a clean control tissue.** rs1198588's eye mean (0.516) exceeds its CNS mean (0.160). Retina is neural-derived developmentally, so this may be biologically sensible rather than contradictory — but it means retina cannot serve as a non-neural negative control. A genuinely unrelated tissue (liver, muscle) would be a better contrast.
4. **Single model fold.** All results use `model_fold_0`. Ensembling across the four folds would indicate how stable these scores are.
### Status
✅ Step 1 complete — padding filtered, CNS vs non-CNS computed, per-track composition examined.
⬜ Step 2 — negative controls.
 
### Next steps
1. **Internal control:** sample non-significant SNPs from PGC3 matched on genomic context (lowest p-value within 1 Mb windows), run the same ISM, build a null for CNS−other
2. **External control:** height GWAS (GWAS Catalog) lead SNPs — a trait with no expected brain-chromatin mechanism
3. Add CHIP_HISTONE (70 CNS tracks) at 128 bp resolution; mask becomes ~4 bins
4. Consider a non-neural control tissue set (liver, muscle) rather than eye

---

## Session 7 — August 31, 2026
**Goal:** Build a first internal negative-control SNP set and run ISM on it (Step 2 of the post-meeting plan)
 
### Why controls
Session 6 left one gap: no reference distribution. A CNS − non-CNS value of +0.10 (rs1198588) cannot be called large or small without knowing what the same statistic looks like for variants with no disease association.
 
### Design decisions
- **Random first, matched later.** Matching on allele frequency and GC content (Li & Ernst, 2025: GC content alone explains much apparent deep-learning variant signal) is the better design but more work. A random set validates the machinery first.
- **n = 10**, the same size as the pilot set. This is a machinery check, not a null: a mean and spread from 10 points is very uncertain. 50+ needed before drawing conclusions.
### Selection criteria (`phase1/get_control_snps.py`)
1. **p > 0.5.** Clearly non-associated. Deliberately far above 5×10⁻⁸: SNPs at p ≈ 10⁻⁵ may be real but underpowered signal.
2. Single-nucleotide variants only.
3. Outside the extended MHC (chr6:25–34 Mb, hg19), as for lead SNPs.
4. ≥ 1 Mb from every one of the 173 lead SNPs, so no control sits in a risk locus's LD block.
5. At most one per chromosome, to avoid clustering.
6. Reservoir sampling (k = 1 per chromosome) while streaming the 7.6M-row file, so memory stays small on the login node. Seed 42.
### Controls selected
Alleles are written here as GWAS roles (A1 = effect allele).
 
| rsid | chr | pos hg19 | pos hg38 | A1 | A2 | p |
|---|---|---|---|---|---|---|
| rs2067987 | 2 | 242,027,926 | 241,088,511 | T | A | 0.584 |
| rs3777506 | 6 | 75,884,136 | 75,174,420 | G | A | 0.832 |
| rs12807017 | 11 | 9,679,145 | 9,657,598 | T | C | 0.957 |
| rs7137768 | 12 | 77,724,773 | 77,330,993 | T | C | 0.536 |
| rs7326762 | 13 | 86,578,698 | 86,004,563 | T | A | 0.539 |
| rs72730413 | 14 | 67,384,310 | 66,917,593 | G | A | 0.529 |
| rs78807154 | 16 | 8,200,245 | 8,150,243 | C | T | 0.698 |
| rs6033101 | 20 | 11,297,353 | 11,316,705 | A | G | 0.606 |
| rs2898160 | 21 | 30,487,496 | 29,115,175 | A | G | 0.931 |
| rs62225795 | 22 | 48,554,947 | 48,159,130 | A | G | 0.659 |
 
The chr6 control (75.9 Mb) is well outside the MHC window. Liftover: 10/10.
 
### ISM attempt failed silently
Ran `run_ism_pilot.py` on the controls. The run completed but **re-scored the 10 pilot SNPs** and overwrote `pilot_ism_scores.tsv` with identical values. Two causes, the first hiding the second:
1. The edits making `run_ism_pilot.py` read its input path from the command line had not been saved, so the hardcoded pilot path was used.
2. The control TSVs existed only on Hoffman2 (never committed and pushed), so the Mac did not have them. Because of (1) the script never tried to open them, so there was no "file not found" error.
**Lesson:** a 15-minute run can produce a plausible-looking wrong file. Scripts should print the input and output paths they actually resolved, before doing any work.
 
### Status
✅ Control selection script written; 10 controls selected and lifted.
⬜ Control ISM run (redone in Session 8).
 
---

## Session 8 — October 10–11, 2026
**Goal:** Re-enter the project after a month away (start of C&S Bio 198A), sync both machines, run ISM on the controls
 
### Inventory of both machines
- Both on the same commit (`34cf0f2`); nothing had diverged.
- Mac: the `run_ism_pilot.py` argument edits **had** saved but were uncommitted.
- Hoffman2: control TSVs untracked; `data/` (PGC3 files, 456 MB) and `phase2/data/` untracked. The `.gz` files are correctly gitignored. Did **not** use `git add -A` on Hoffman2, to avoid pushing large files.
### Fixes
- Pushed the control TSVs from Hoffman2 (`19f15ed`).
- Committed the argument edits (`b6bae16`). The first push was rejected because the remote had the control commit; resolved with `git pull --rebase`.
- Added a printout of the resolved input and output paths at the start of `main()` (committed later, `6393041`).
- **Environment:** `python3` on the Mac now resolves to Homebrew's Python, which has no torch (`ModuleNotFoundError`). Fix: always `source ~/alphagenome-env/bin/activate` first. `conda` and `module` are Hoffman2 commands and do not exist on the Mac.
### Control ISM run
```bash
PYTORCH_ENABLE_MPS_FALLBACK=1 python3 run_ism_pilot.py \
    phase1/data/scz_pgc3_control_10snps_hg38.tsv control_ism_scores.tsv
```
6,400 scores to `results/ism/control_ism_scores.tsv` (committed in Session 10).
 
Summary over all output channels (padding not yet removed, so means are diluted, as in Session 4):
 
| SNP | ATAC max | ATAC mean | DNase max | DNase mean |
|---|---|---|---|---|
| rs2067987 | 0.3079 | 0.0457 | 0.3412 | 0.0376 |
| rs3777506 | 0.0925 | 0.0051 | 0.1518 | 0.0209 |
| rs12807017 | 0.0380 | 0.0065 | 0.0490 | 0.0185 |
| rs7137768 | 0.1073 | 0.0107 | 0.2431 | 0.0488 |
| rs7326762 | 0.0222 | 0.0057 | 0.0418 | 0.0031 |
| rs72730413 | 0.7608 | 0.0303 | 0.9421 | 0.1074 |
| rs78807154 | 0.1396 | 0.0095 | 0.2071 | 0.0247 |
| rs6033101 | 0.1195 | 0.0128 | 0.3394 | 0.0261 |
| rs2898160 | 0.1499 | 0.0189 | 0.1552 | 0.0386 |
| rs62225795 | 0.0767 | 0.0261 | 0.1614 | 0.0575 |
 
### Observations
- **hg38 reference = A1 for all 10 controls**, as for all 10 pilots: 20/20. This is not specific to lead SNPs; it is a property of how PGC3 reports alleles.
- **The null has a tail.** rs72730413 (p = 0.53) has DNase max 0.94 and mean 0.107 (all channels), both larger than 9 of the 10 pilot SNPs; only rs2710323 is higher. Magnitude alone will not separate disease variants from controls; tissue specificity has to.
### Review before the Choo meeting (October 11–12)
Re-derived the whole project logic step by step (written up as the "Getting Recalibrated" slides, week of 10/12). Points worth keeping:
- **MHC window.** chr6:25–34 Mb (hg19) is a rounded, conservative version of the extended MHC (xMHC, Horton et al., 2004, *Nat Rev Genet* 5:889–899; ~25.7–33.4 Mb). The rounded window is the default in GWASLab and is used in SCZ work (e.g. *Transl Psychiatry* 2022, doi:10.1038/s41398-022-01854-9). Other studies use different windows (e.g. chr6:26–34 Mb), so it is a stated choice. 20,457 − 5,994 MHC SNPs = 14,463 exactly, so no other filter removed SNPs.
- **Padding.** Every padded head is padded to the next multiple of 128 (ATAC 167→256, DNase 305→384, histone 1,116→1,152). Future code should drop padding by the metadata name `Padding`, not by index thresholds, and assert the real-track counts.
- **Unresolved:** our CNS counts differ from Choo's table (DNase 23 vs 21; histone ChIP 70 vs 100, probably cell lines). To be settled as a lab group.
- **Candidate specificity statistics:** CNS − non-CNS (scales with overall loudness), log2(CNS mean / non-CNS mean) (scale-free), and mean CNS percentile rank within the variant (exact null 0.5). Choose empirically against controls.
- **What each control answers:** random internal controls give the scale of an ordinary variant but sit in random genomic context; an external trait (height) controls for "GWAS hits sit in regulatory DNA"; the within-variant CNS vs non-CNS contrast cancels locus context. The height comparison is the one that supports a disease-specific claim.
### Status
✅ Both machines in sync; control ISM scores produced.
⬜ Control analysis (needs `score_tracks.py`).
 
---

## Session 9 — October 12–13, 2026
**Goal:** Present the recalibration to Choo; plan the week
 
### Feedback from Choo (October 12 meeting)
1. **Record the A1/A2 vs hg38 ref/alt mismatch per SNP.** Direction matters: if hg38 has A at a position and the GWAS says A raises risk, then substituting the alternate allele *reduces* risk.
2. **Understand what the variant-scoring summation does and means.**
    - summing signal in the 501 bp window to get an average of signal around the variant
3. **Brain vs non-brain tissue classification will be decided as a lab group** at the next group meeting.
4. **Plot before summarising.** Show variant score (y) against CNS vs non-CNS (x); compare means *and variances* as an intermediate step; see whether a high-variance score is driving things and which statistic could represent it.
5. **Work out the edge cases of the liftover.**
6. **Build visualisations** that give users intuition about the range and variability of scores.
Also: send Choo the hg19 → hg38 pipeline; hand the pipeline to another undergraduate in about two weeks to test.
 
### Decisions (October 13)
- **Build risk direction into the pipeline now**, not later, documenting exactly how GWAS alleles are substituted.
- **This week:** Session C = liftover package (strict script, edge cases, documentation); Session D = `score_tracks.py` data layer, distribution plots, scoring sensitivity.
- **Deferred:** height GWAS control and a 50-SNP matched control set, until the score itself is understood.
---

## Session 10 — October 13–14, 2026
**Goal:** Session C, steps C0–C2: audit the existing liftover code and replace it with one strict, documented script
 
### C0 — Environments
- Committed the outstanding work: path printout (`6393041`), control ISM scores (`4aaf666`), notebook/methods edits (`aa0afbc`). Hoffman2 pulled to match.
- Mac had no chain file and no pyliftover. Installed `pyliftover 0.4.1` in `alphagenome-env`; downloaded `hg19ToHg38.over.chain.gz` (222 KB) to `phase1/data/` (gitignored).
- Hoffman2: chain file present; `pyliftover` imports. Note: `which python` gives the base Anaconda Python even with `alphagenome_env` active.
- Do not run `find ~` on the Mac: it walks iCloud folders and hangs.
### C1 — Audit of the old code
Liftover logic lived in two places: `phase1/liftover_to_hg38.py` and the end of `phase1/get_control_snps.py`.
 
| Problem | Consequence |
|---|---|
| `convert_coordinate()` returns a **list**; code took `[0]` without checking its length | a SNP mapping to several places silently gets one of them |
| strand returned but ignored | a SNP in a region inverted in hg38 passes; its hg38 base is the complement of the hg19 base |
| `replace('chr', '')` on contig names | `chr6_GL000256v2_alt` would become chromosome `6_GL000256v2_alt` |
| failed SNPs dropped, exit status 0 | a short output file looks as valid as a complete one |
| chain path differed between the two files and matched neither real location | worked only from one particular folder |
| `liftover_to_hg38.py` hardcoded to the 2018 CLOZUK files | **the PGC3 hg38 coordinates used since August had no committed script** |
 
### C2 — New `phase1/liftover_to_hg38.py`
Read pyliftover's source to confirm its behaviour: it returns `None` for an unknown chromosome, `[]` for an unmapped position, and otherwise a list sorted by chain score.
 
- **Command-line tool:** `--input`, `--output`, optional `--chain`, `--chr-col`, `--pos-col`, `--id-col`, `--report`, `--allow-failures`.
- **Strict by default.** Writes a SNP only if it maps uniquely, to a standard chromosome, to the same chromosome, on the + strand. Otherwise it records one of 7 failure statuses: `unmapped`, `multi`, `strand_flip`, `alt_contig`, `chrom_change`, `unknown_chrom`, `invalid_input`. Any failure means no output file and exit status 1; `--allow-failures` overrides.
- **Guards:** refuses to re-lift a lifted file; refuses to overwrite its input; names missing columns; prints the chain download command.
- **Outputs:** lifted TSV (input columns, hg38 `chr`/`pos`, plus `chr_hg19`/`pos_hg19`); per-SNP report; provenance log (command, chain md5, versions). The log is a separate file because comment lines at the top of a TSV would be read as the header.
- 1-based ↔ 0-based conversion lives only in `lift_one()`.
- Allele columns are not touched: liftover moves coordinates only.
### Testing
- **Synthetic chain + 13 SNPs:** every status, the last base of a mapped block and the first base of a gap, and invalid input (chromosome 23, text position, position 0). All correct. Each guard triggered and stopped the script.
- **Real data:** re-lifted the 173 PGC3 lead SNPs (hg19) to a temporary file: 173 lifted, 0 failed, **0 positions differ** from the August file. The August PGC3 liftover is now reproducible from committed code.
### Status
✅ Strict liftover script written and validated.
 
---
 
## Session 11 — October 16, 2026
**Goal:** Make every hg38 file come from the new script, remove the second liftover implementation, and make allele handling explicit
 
### Regenerated hg38 files
- Committed the script (`05a120a`).
- Re-lifted the lead, pilot and control files: **173/173, 10/10, 10/10 positions identical** to the committed versions (differ 0, missing 0, extra 0). Committed with their reports and logs (`692bf7d`). Every line shows as changed in git because of the new `chr_hg19`/`pos_hg19` columns and the switch from `\r\n` to `\n` line endings.
### One liftover implementation (C2b)
- Removed the embedded liftover from `get_control_snps.py`, which now writes hg19 only and prints the liftover command to run next (`79a5261`). It no longer needs pyliftover.
- **Regression test:** re-ran control selection on Hoffman2 with seed 42; the output was **byte-identical** to the committed control file, so only the liftover code was removed.
### Allele refactor
**Problem:** both selection scripts wrote A2 to a column named `ref` and A1 to `alt`. `run_ism_pilot.py` corrected this at run time, so every score was right, but the files said something false.
 
**Direction formula, corrected.** The earlier rule `risk_oriented_score = −sign(beta) × score` only holds when hg38 ref = A1. General form:
```
risk_oriented_score = sign(beta) × variant_effect_score × (+1 if hg38_alt = A1, −1 if hg38_alt = A2)
```
i.e. log2(signal with risk allele ÷ signal with non-risk allele). Blank when beta = 0.
 
**Palindromic SNPs.** For A/T and C/G SNPs, matching to hg38 cannot reveal which strand the GWAS used. Magnitude is unaffected; direction assumes the + strand. Supported by 17/17 non-palindromic SNPs matching on the + strand (a minus-strand non-palindromic SNP would match neither allele). **rs1198588 (best pilot candidate) is palindromic**, as are controls rs2067987 and rs7326762.
 
**A1 — selection scripts.** Write `a1`, `a2`, and `a1_freq_controls` (PGC3 `FCON`, frequency of A1 in controls; needed later for allele-frequency-matched controls and for checking palindromic strand). New column order: `rsid chr pos a1 a2 a1_freq_controls beta pval`.
 
**A2 — regenerate and verify (Hoffman2).** Check against the committed files: rsid, position, beta and p unchanged; new `a1` = old `alt`; new `a2` = old `ref`; frequency numeric in [0, 1].
- First run: lead, pilot and MHC files OK; **control file had a blank frequency in all 10 rows.** Cause: the frequency line was missing from the control script's row dictionary, and `csv.DictWriter` **silently writes '' for a column missing from a row** (whereas an *extra* key raises an error).
- Fix: added the line, plus a guard in both selection scripts that stops if any output column is missing from a row. Re-run: all four files OK. Re-lifted.
**A3 — `run_ism_pilot.py`.**
- Reads `a1`, `a2`, `beta`. Refuses old `ref`/`alt` files *before* loading the model.
- `assign_alleles()` (pure function) decides `hg38_ref`, `hg38_alt`, `hg38_ref_is` (A1/A2), `palindromic`, the risk allele and the orientation.
- New output columns: `a1 a2 beta hg38_ref hg38_alt hg38_ref_is palindromic … variant_effect_score risk_oriented_score`. New per-SNP file `*.variants.tsv` with status (`scored`/`skipped`) and reason. Exits 1 if any SNP is skipped, but still writes the rest.
- Sequence fetch, one-hot encoding and `compute_variant_effect()` unchanged.
- Fixed two latent bugs: an indel crashed the run before anything was saved; a `chr8`-style name was passed to Ensembl unchanged.
- Tested with a stand-in model on 13 synthetic SNPs (every A1/A2 × beta-sign combination, palindromic, beta = 0 and NA, lowercase alleles, and five SNPs that must be skipped). Scores and risk scores matched hand calculations, and matched the original script's output exactly; the old `gwas_alleles_swapped` flag equals `hg38_ref_is == A1` on every row.
**A4 — real-data regression.** Re-ran ISM on pilots and controls with the new script and compared track by track with the committed scores: **all 12,800 scores matched**, alleles matched, `hg38_ref_is` = A1 for all 20 SNPs, palindromic = rs1198588, rs2067987, rs7326762. Committed scores plus `*.variants.tsv`.
 
### Status
✅ One liftover implementation; all hg38 files reproducible.
✅ Every score records its alleles, hg38 assignment and risk direction.
 
---
 
## Session 12 — October 18, 2026
**Goal:** Measure how often each liftover edge case occurs in real data (C3) and document the package (C4)
 
### Why a separate survey script
`liftover_to_hg38.py` holds all rows in memory (right for hundreds of SNPs; would hit the August login-node memory error on 7.6M). `phase1/survey_liftover_edge_cases.py` streams the GWAS file and keeps only counts, using the **same `lift_one()`**, so it tests the pipeline's own code. Tested on a synthetic PGC3-format file (including p = 6×10⁻⁸ correctly *not* counted as significant and p = "NA" counted in neither p-value subset). Python 3.6+ compatible.
 
### Run
Hoffman2 compute node (`qrsh -l h_data=4G,h_rt=1:00:00`), pyliftover 0.4.1, **71 s**.
 
| Status | All SNPs | p < 5×10⁻⁸ | p > 0.5 |
|---|---|---|---|
| lifted | 7,646,680 (99.829%) | 20,434 (99.888%) | 3,099,610 (99.849%) |
| strand_flip | 8,304 (0.108%) | 12 | 2,712 |
| unmapped | 3,809 (0.050%) | 10 | 1,559 |
| alt_contig | 948 (0.012%) | 1 | 418 |
| chrom_change | 26 | 0 | 11 |
| multi | 0 | 0 | 0 |
| **total** | **7,659,767** | **20,457** | **3,104,310** |
| palindromic | 1,159,302 (15.1%) | 3,022 (14.8%) | 468,739 (15.1%) |
 
Sanity checks: 7,659,767 rows and 20,457 significant SNPs both match the counts from Session 3. No non-numeric p-values.
 
### Where the 23 significant failures are
- **12 in the excluded MHC**, all strand flips, at chr6:26.74–26.76 Mb (hg19).
- **11 within 1 Mb of a lead SNP that lifted**: 8 near rs5751191 (chr22:42.53–42.55 Mb; 7 unmapped, 1 alt contig), 2 near rs12489270 (chr3:16.86 Mb), 1 near rs12138231 (chr1:150.15 Mb).
- **0 elsewhere.** No schizophrenia locus is lost. These 11 would matter only if every SNP in a locus were scored, not just the lead.
### Palindromic share of each failure type
strand_flip 1,383/8,304 (16.7%); unmapped 606/3,809 (15.9%); alt_contig 159/948 (16.8%); chrom_change 5/26 (19.2%). All close to the 15.1% background.
 
### Interpretation
- **Strand flips are the most common failure (63%).** The old code would have passed all 8,304. About 1 in 6 is palindromic, where the hg38 base matches the *other* GWAS allele and the risk direction would be silently reversed. Non-palindromic ones would have made `run_ism_pilot.py` stop ("matches neither").
- **No multi-mappings** with this chain file.
- **Planning for more controls:** 0.151% of null SNPs fail, so a random 50 has ~7% chance (100: ~14%) of including one that cannot be lifted. Strict mode will name it; replace it rather than relax the check.
### Documentation and test package
- `docs/liftover_README.md`: why, setup, usage, strict rules, coordinate convention, edge-case table with these frequencies, what liftover does not check, validation.
- `docs/allele_conventions.md`: A1/A2 vs hg38 ref/alt vs risk allele, column glossary, both score formulas with a worked example (rs58120505), palindromic SNPs and how to verify their strand.
- `tests/liftover/check_liftover.py` + `edge_cases_hg19.tsv`: six real PGC3 SNPs with known outcomes (rs58120505 → chr7:1,990,232; rs1198588 → chr1:98,087,276; rs748439 unmapped; rs4993393 alt contig; rs13213200 strand flip; one invalid) plus every guard; 22 checks.
- `docs/TESTING_liftover.md`: step-by-step instructions and feedback questions for lab testers.
### Status
✅ Liftover package written, validated, surveyed genome-wide and documented.
⬜ Tag `liftover-v1.0` and share with lab testers.
⬜ Next steps: `score_tracks.py` (data layer), distribution plots, scoring-window sensitivity.
 
---

## Session 13 (D1) — October 14–21, 2026
**Goal:** Work out what every AlphaGenome output track is and label each one. This answers Choo's Slack question (Oct 14: "what alphagenome track is which? how do i tell which is brain or ipsc?") and her 10/12 note to classify tissues as a lab group.
 
**Decision before starting:** research the classification ourselves and build a preliminary label file plus the pipeline that uses it. The lab group then reviews the decisions rather than starting from scratch.
 
### Sources
- **Track metadata:** `google-deepmind/alphagenome_research`, `OutputMetadataResponse_ORGANISM_HOMO_SAPIENS.textproto` (md5 `474494fa3859c361d40a735b8d9a531a`). Metadata row *i* = output channel *i* of that head.
- **Ontologies:** Uberon anatomy (release 2026-10-01) and Cell Ontology (release 2026-06-08), `basic.obo` versions. Stored in `phase2/metadata/ontologies/` (gitignored, large).
- **Cellosaurus:** for cell-line identity and known problems.
- **ENCODE biosample pages:** for the parent line and donor of in-vitro samples.
### What a track is
- **Heads and channels.** Each head is one assay, and each channel is one training experiment.
- **Padding.** Heads are padded with all-zero channels named `Padding` up to a multiple of 128, always as one block at the end:
  | Head | Real tracks | Channels |
  |---|---|---|
  | ATAC | 167 | 256 |
  | DNase | 305 | 384 |
  | CAGE | 546 | 640 |
  | RNA-seq | 667 | 768 |
  | Histone ChIP | 1,116 | 1,152 |
  | TF ChIP | 1,617 | 1,664 |
  | PRO-cap | 12 | 128 |
  6,126 channels in total, 563 of them padding.
- **Fields.** Each real track has an ontology term, a biosample type, a biosample name, a life stage (set **per track**), the assay, a histone mark or TF, and a data source.
- **Ontology prefix = kind of sample:**
  - UBERON: tissue
  - CL: cell type
  - EFO/CLO: cell line
  - NTR: provisional ENCODE term
- **Biosample type = where the cells came from:**
  - TISSUE
  - PRIMARY_CELL
  - CELL_LINE
  - IN_VITRO_DIFFERENTIATED_CELLS
  - ORGANOID
### How to tell brain from iPSC (the answer to Choo's question)
- **Real brain** = TISSUE or PRIMARY_CELL whose ontology term is part of the central nervous system.
- **Stem-cell-derived neural** = IN_VITRO_DIFFERENTIATED_CELLS. These share names with real cell types (e.g. "neural progenitor cell"), so only the type separates them.
- **Pluripotent lines** (H1, H9, WTC11, iPS lines) = CELL_LINE.
- **ESC- vs iPSC-derived:** the stage of an in-vitro sample is the stage of the **parent line's donor**. `embryonic` → ESC-derived; `adult`/`child` → iPSC-derived.
  - Confirmed on ENCODE: motor neurons are iPSC-derived from a 49-year-old donor (ENCBS745APJ); glutamatergic neurons "originated from WTC11" (ENCBS229NUR).
- **Fetal vs adult:** read the stage per track, because the same sample name can be embryonic for one assay and adult for another. For example, "brain" DNase is embryonic, while some "brain" histone tracks are adult.
### Method: `phase2/build_track_classes.py`
1. Parse the metadata and the two ontologies.
2. **Find candidates with the ontology, not names.** Candidates are tracks whose term is part of the CNS (UBERON:0001017), PNS (UBERON:0000010) or eye (UBERON:0000970), or whose cell type descends from a neural cell type (CL:0002319, CL:0000540, CL:0000125, CL:0011020, CL:0000047, CL:0000133, CL:0011012). Names are used only to flag extra candidates for review.
3. **Every candidate needs a reviewed row in `phase2/metadata/biosample_curation.tsv`.** That file has 132 rows, each with a group, cell class, brain region and evidence. Rows match on (ontology ID, biosample type, name), with `*` as a wildcard; the most specific row wins.
4. **The build stops** if:
   - an ontology candidate has no reviewed row,
   - a curated row matches no track, or
   - two rows match the same sample equally well.
   So new metadata can't change the brain set unnoticed.
5. Everything else is `non_neural`.
6. **Output:** `phase2/metadata/track_classes.tsv`, one row per channel for all heads. Columns: `output_type, track_index, track_name, is_padding, ontology_curie, biosample_type, biosample_name, stage, life_stage_group, group, analysis_set, origin, cell_class, brain_region, assay, histone_mark, tf, data_source, gtex_tissue, strand, evidence`. Plus `.summary.tsv` and `.log.txt` (input md5s and ontology versions).
7. **To change a decision:** edit the row's `group` and rerun. No code changes are needed.
### Problems found and fixed while building
- **Ontology parser bug.** My first OBO parser kept qualifiers such as `{is_inferred="true"}` inside the parent ID, so some `is_a` links were silently lost. After the fix, the build found 7 candidates the curation had missed and stopped, as designed:
  - temporal lobe, parietal lobe and corpus callosum → added as CNS_tissue
  - four extraocular muscles (caught by the eye search) → reviewed as non_neural
- **Cell-line IDs.** Some EFO IDs I had first written into the curation for cell lines could not be verified. I replaced them with exact name matching (`*` ID). Only the IDs I checked are kept: EFO:0002860 SK-N-MC, EFO:0003072 SK-N-SH, EFO:0005721 SK-N-DZ, EFO:0002782 H54.
### Traps found
- **Name matching gives false hits:** "g**astro**cnemius", "renal **cortex**", "P**eye**r's patch".
- **Retina** is filed under eye in Uberon, not CNS.
- **Located in the brain is not the same as neural:** brain microvascular endothelium, pericytes and vascular smooth muscle are blood-vessel cells.
- **Neuroblastoma lines** (SK-N-SH, SK-N-DZ, BE2C, SH-SY5Y) are peripheral nervous system / neural crest, not brain.
- **SK-N-MC** is an Askin tumour, not a neuroblastoma (Cellosaurus CVCL_0530).
- **H54** is an A-172 derivative (CVCL_7185).
- **UCSF-4** is an ESC line (CVCL_A022) and **L1-S8** is an iPSC line (CVCL_AV88).
- **Tumour lines in the curation:** D721Med (medulloblastoma), PFSK-1 (PNET), WERI-Rb-1 (retinoblastoma).
- **Cell-type names carry no location:** in the Cell Ontology, "neuron" and "astrocyte" are not linked to the CNS, and "neural progenitor cell" is not even classed under neural cell.
### Groups → analysis sets
- **brain:** CNS_tissue, CNS_primary_cell
- **neural_model:** neural_in_vitro, CNS_organoid
- **excluded:** neural_crest_in_vitro, CNS_non_neural, CNS_tumor_line, PNS, PNS_tumor_line, eye, pluripotent
- **background:** non_neural
"Excluded" samples are neither clearly brain nor clearly unrelated to brain, so they sit out of **both** sides of the comparison. They can still be analysed as their own groups.
 
### Counts (real tracks)
 
| | ATAC | DNase | Histone | RNA-seq | CAGE | TF ChIP |
|---|---|---|---|---|---|---|
| brain | **0** | 16 | 70 | 41 | 66 | 4 |
| neural_model | 1 | 4 | 32 | 9 | 0 | 14 |
| excluded | 3 | 30 | 128 | 43 | 36 | 134 |
| background | 163 | 255 | 886 | 574 | 444 | 1,465 |
 
- **ATAC has no brain tracks.** Its one neural track is iPSC-derived motor neurons. This is the strongest argument for fine-tuning on fetal brain ATAC.
- **Brain by life stage:**
  - DNase: 2 prenatal, 10 postnatal, 4 unknown (cultured astrocytes)
  - Histone: 21 prenatal, 37 postnatal, 12 unknown
- **Brain histone tracks span 11 marks** (e.g. 10 H3K27ac, 15 H3K4me3, 6 H3K27me3), so histone must be analysed **per mark**. Active and repressive marks must not be averaged together.
### Changes from the August DNase "CNS" list (23 → 16 brain)
- **Kept (12):**
  - astrocyte; astrocytes of the cerebellum, hippocampus and spinal cord
  - brain (embryonic), spinal cord (embryonic)
  - frontal cortex, dorsolateral prefrontal cortex
  - putamen, head of caudate nucleus
  - cerebellum, cerebellar cortex
- **→ neural_model (4):** neuronal stem cell, neural progenitor cell, ecto neural progenitor (ESC-derived); glutamatergic neuron (iPSC-derived).
- **→ excluded (7):** 3 brain blood-vessel cell types; neural crest cell; SK-N-SH; SK-N-DZ; SK-N-MC.
- **Added (4), missed in August:** pons, globus pallidus, Ammon's horn, posterior cingulate gyrus.
- **Session 1 table comparison:** the Session 1 brain table (21 DNase / 100 histone) is close to our brain + neural_model (20 / 102). The difference is probably mostly whether stem-cell-derived neural cells count as brain. To check track by track with Choo.
### Run on the Mac
`python phase2/build_track_classes.py` → all curated rows matched, no unreviewed candidates, no ambiguous matches. Output counts are identical to the sandbox build.
 
### Documentation
`docs/alphagenome_tracks_guide.md` covers:
- what a track is
- the metadata fields
- how to tell brain from iPSC
- the traps
- the label scheme and counts
- how to change a decision
- changes from August
- five open questions for the lab group
- sources
Written so other lab members (e.g. Nadia) can use it.
 
### Open questions for the lab group
1. Do stem-cell-derived neural cells count as brain?
2. Brain tumour lines: exclude (current) or keep as their own group?
3. FANTOM primary neurons and neuronal stem cells with no location (CAGE only): brain?
4. Split prenatal and postnatal brain? DNase has only 2 prenatal tracks.
5. Is pituitary correctly excluded?
### Status
✅ Every output channel labelled; brain vs stem-cell-derived vs cell line answered and documented.
⬜ Lab-group review of the 5 questions; then edit the curation file and rebuild.
 
---
 
## Session 14 (D2) — October 22–29, 2026
**Goal:** Join the ISM scores to the new labels in one place, then follow Choo's 10/12 notes:
- plot score (y) by CNS / non-CNS (x)
- compare means **and variances**
- plot the distributions before summarizing
- visualize the range and variability of scores
### `phase2/score_tracks.py`: the single place where scores meet labels
- **Input:** any number of ISM score files (`--scores pilot=… control=…`) and `track_classes.tsv`.
- **Checks (stops if any fails):**
  - every score row has a label, so the channel numbering is intact
  - padding is one block at the end of each head
  - after dropping padding, each SNP has exactly the head's real-track count (305 DNase, 167 ATAC), so a truncated or duplicated run can't pass
  - no duplicate rows
  - each SNP has one set of alleles across all rows
  - old-format files (`true_ref`/`true_alt`) are refused with a hint
  - set names must be unique
- **Outputs:**
  - `results/summaries/track_scores_long.tsv`: one row per set × SNP × real track, with SNP fields, label fields, score, |score| and risk-oriented score
  - `track_scores_by_set.tsv`: per SNP × head × analysis set, with n, mean, variance, median and max of |score|, and mean risk-oriented score
  - a log with input md5s
- **Run:** pilot and control files each have 6,400 rows. 1,680 padding rows were dropped from each, leaving 4,720 per set and 9,440 in total (6,100 DNase + 3,340 ATAC).
- **Regression check:** recomputing the August CNS vs non-CNS summary from the new table with the **August** CNS list reproduced the August numbers for **10/10 pilots**. So any change in results comes from the new labels, not from the code.
### `phase2/plot_track_scores.py`: look before summarizing
Reads only the tidy table and draws DNase by default. |score| is on a log axis shared across panels; values below 10⁻⁴ are drawn at 10⁻⁴. Colours: brain blue, neural_model orange, excluded green, background grey.
1. `dnase_1_per_snp_strips.png`: one panel per SNP. Every track's |score| is a dot, split by analysis set, with each set's mean marked. Pilots come first, then controls, each sorted by brain ÷ background.
2. `dnase_2_means_variances.png`: brain vs background mean |score| (left) and variance (right), one point per SNP, pilots vs controls.
3. `dnase_3_score_ranges.png`: each SNP's full spread across all tracks, with brain tracks overlaid.
4. `dnase_4_top_tracks_<rsid>.png`: the tracks with the largest |score|, coloured by set. Drawn for rs1198588, rs2710323, rs2238057 and rs72730413.
Plot fixes made along the way: overlapping per-row y labels were replaced by one shared label, and a colliding label was moved.
 
### Results (DNase; 16 brain vs 255 background tracks; mean |score| with variance in brackets)
 
| Set | SNP | Brain | Background | Ratio |
|---|---|---|---|---|
| pilot | rs1198588 | 0.219 (0.0369) | 0.047 (0.0030) | **4.62** |
| pilot | rs2238057 | 0.075 (0.0035) | 0.044 (0.0036) | **1.73** |
| pilot | rs58120505 | 0.018 (0.0001) | 0.016 (0.0001) | 1.17 |
| pilot | rs11191580 | 0.084 (0.0006) | 0.072 (0.0016) | 1.16 |
| pilot | rs4702 | 0.044 (0.0023) | 0.049 (0.0026) | 0.90 |
| pilot | rs12129573 | 0.035 (0.0032) | 0.045 (0.0040) | 0.78 |
| pilot | rs13107325 | 0.031 (0.0023) | 0.045 (0.0028) | 0.69 |
| pilot | rs4129585 | 0.045 (0.0013) | 0.077 (0.0042) | 0.58 |
| pilot | rs778371 | 0.009 (0.0001) | 0.018 (0.0003) | 0.49 |
| pilot | rs2710323 | 0.154 (0.0216) | 0.576 (0.4203) | **0.27** |
| controls | 10 SNPs | | | 0.56–1.13 |
 
### Interpretation
- **Only rs1198588 and rs2238057 are clearly above every control.** rs58120505 (1.17) and rs11191580 (1.16) are only just above the top control (1.13). With 10 controls this is descriptive, not a test.
- **rs1198588's ratio rose from 2.73 (August labels) to 4.62.** This comes only from cleaning the labels: removing neuroblastoma and vascular tracks, and adding the basal ganglia tracks August had missed.
- **Variance shows what the mean hides.** rs1198588's brain variance is 12× its background variance (0.0369 vs 0.0030), so the brain tracks are not one population:
  - **High** in adult neuron-containing tissue: globus pallidus 0.53, putamen 0.51, pons 0.43, caudate 0.42, Ammon's horn 0.38, posterior cingulate 0.35, DLPFC 0.33.
  - **Low** in cultured astrocytes (≤ 0.09), cerebellum (0.02) and fetal brain (0.02).
  - High in ESC-derived neural progenitors (0.46) but with the **opposite sign**: adult tissue +, NPCs −.
  - Its top 30 tracks also include eye/RPE, retinoblastoma (WERI-Rb-1), medulloblastoma (D721Med) and Ewing-family lines (A673, SK-N-MC). "Neural/neuroectodermal" may describe it better than brain-only.
  - It is palindromic, so the magnitude is solid but the **direction assumes the + strand**.
- **rs2238057** is highest in cerebellum (cerebellar cortex 0.23, cerebellum 0.20), with all brain tissue in the same direction.
- **rs2710323** has the largest effects of any SNP (max 4.42, immature NK cells), but **all of its top 30 tracks are immune or blood** (T, NK and myeloid cells, thymus). This points to a non-brain mechanism.
- **Control rs72730413** has large effects (max 0.94, bigger than 9 of 10 pilots) but is not brain-enriched (0.57). Big effect is not the same as brain effect.
- **Cultured astrocytes are low for every pilot.** Question for the lab: are cultured astrocytes a good proxy for brain astrocytes?
### Status
✅ Scores joined to labels with checks; August numbers reproduced 10/10.
✅ Distributions, means, variances, ranges and top tracks plotted (Choo's 10/12 notes).
⬜ Choose a specificity statistic; build a larger matched control set.
⬜ Decide whether A673 and SJCRH30 belong in background.
 
---

## Session 15 (D3) — October 29 – November 7, 2026
**Goal:** Understand exactly what the variant effect score measures (Choo, 10/12: "understand how the variant-scoring summation works") and test whether our conclusions depend on how it is computed

### What the score is
- `score = log2(1 + Σ ALT) − log2(1 + Σ REF)`, Σ over the 501 bp centred on the variant (positions 65,286–65,786 of the 131,072 bp input, 0-based), per track, at 1 bp resolution.
- Meaning: fold change in the total predicted accessibility of the element around the variant. +1 = ALT doubles it, −1 = halves it.
- The +1 pseudocount avoids log(0) and damps fold changes where the signal is small (Σ 0.5 → 1.0 gives +0.42 rather than +1.00).
- **Verified against source code:**
  - AlphaGenome's `RECOMMENDED_VARIANT_SCORERS['DNASE']` is `CenterMaskScorer(width=501, aggregation_type=DIFF_LOG2_SUM)`.
  - `alphagenome-pytorch` 0.3.1 defines `DIFF_LOG2_SUM` as `log2(1 + alt_sum) − log2(1 + ref_sum)`, and its window arithmetic (`create_center_mask`) is identical to ours.
  - **Our score is AlphaGenome's own recommended DNase scorer.**
- **Correction to Session 1 notes (July):** the "AggregationType naming convention" paragraph calls this aggregation `DIFF_SUM_LOG2`. That is a different scorer (Σ of per-base log2(x+1) differences). Ours is `DIFF_LOG2_SUM`.
- Written up in `docs/variant_scoring.md`.

### Sensitivity analysis: `phase2/scoring_sensitivity.py`
- **SNPs:** rs1198588, rs2238057 (brain-leaning pilots), rs2710323 (immune pilot), rs72730413 (control with large effects). One model run (ref + alt) per SNP on the Mac (MPS).
- **21 ways of rescoring all 305 real DNase tracks:**
  - windows 1, 11, 101, 501, 2001, 10001, 100001 bp and the whole input
  - pseudocounts 0.01, 0.1, 1, 10, 100
  - six aggregations: log of sums (ours), sum of logs, difference of means, difference of sums, L2 of log1p, log ratio of peaks
  - 501 and 2001 bp windows on the 128 bp output
- **Metrics per setting:**
  - Spearman of |score| vs 501 bp
  - top-20 overlap
  - sign agreement (tracks with |score| > 0.001)
  - brain ÷ background mean |score|
  - brain tracks' mean percentile within background
- **Built-in checks.** The script stops unless its baseline equals the pipeline's score, the library scorer, and the committed scores.
- **Tested before running.** Used a stand-in model with the real track labels. The window arithmetic was compared with a port of the library mask at every width and resolution, each aggregation was checked by hand, and a tampered score file was correctly refused.
- **Pre-specified criterion (Oct 29, before running):** a conclusion is robust if it holds for windows 101–2001 bp, pseudocounts 0.1–10 and the sum-based aggregations.

### Checks
- Baseline vs `run_ism_pilot.compute_variant_effect`: max |diff| 2.6–3.7 × 10⁻⁷ (float64 vs float32).
- Today's run vs the committed October scores: max |diff| 4.9–5.0 × 10⁻⁷, which is the 6-decimal rounding. **Scores are reproducible from the weights.**
- Library check skipped: the installed `alphagenome-pytorch` has no `variant_scoring` module.

### Results (brain ÷ background / brain percentile)

| Setting | rs1198588 | rs2238057 | rs2710323 | rs72730413 (control) |
|---|---|---|---|---|
| 501 bp (ours) | 4.62 / 0.72 | 1.73 / 0.73 | 0.27 / 0.23 | 0.57 / 0.40 |
| 101 bp | 2.32 / 0.72 | 1.33 / 0.69 | 0.32 / 0.23 | 0.54 / 0.33 |
| 2001 bp | 6.32 / 0.77 | 2.05 / 0.74 | 0.27 / 0.25 | 0.58 / 0.45 |
| all 21 settings | 1.29–6.97 / 0.61–0.85 | 0.85–2.05 / 0.46–0.79 | 0.15–0.51 / 0.17–0.42 | 0.34–0.95 / 0.33–0.62 |

### Interpretation
1. **Pseudocount is irrelevant** between 0.01 and 10 (Spearman ≥ 0.93; signs 100%; pilot ratios change < 10%).
2. **Every classification holds within the pre-specified neighbourhood:**
   - rs1198588 is brain-enriched at all 21 settings.
   - rs2710323 is brain-depleted at all 21.
   - The control is never brain-enriched by ratio.
3. **rs2238057 is brain-leaning, not brain-specific.**
   - Enriched with fold-change scores at windows 101–10001.
   - With L2, peak or absolute-difference scoring the ratio drops to 0.94–1.11, though brain tracks still rank above the background median (percentile 0.58–0.70).
   - At the extreme windows (1 bp, 11 bp, ≥ 100 kb) the enrichment disappears.
4. **Ratio size depends on the window** (rs1198588: 2.32 → 6.32 from 101 to 2001 bp). The direction of the conclusion doesn't.
5. **rs1198588's exact track ranking is the most window-sensitive** (Spearman 0.70 at 101 and 2001 bp vs 0.85–0.99 for the other SNPs). Report track *sets*, not one top track.
6. **The percentile is steadier than the mean ratio** (rs1198588: 0.61–0.85 vs 1.29–6.97). This argues for a rank/percentile-based specificity statistic.
7. **Absolute-difference scores shrink brain enrichment** because they up-weight tracks that are highly accessible at baseline. Fold change compares relative change per tissue, which is the right question for tissue comparisons.
8. **128 bp bins keep set-level conclusions** but agree less on rankings (Spearman 0.40–0.93) and direction (63–92%). That is relevant for histone later.

### Decision
Keep the pre-specified 501 bp `DIFF_LOG2_SUM` score with pseudocount 1. Don't switch to 2001 bp because it gives larger ratios: that would be choosing the analysis after seeing the result. 2001 bp stays reserved for histone marks, as AlphaGenome recommends.

### Follow-up (Nov 8): `phase2/sensitivity_followup.py`
- **128 bp bins hold sums, not averages.** Over the same 512 bp, the bins add up to 0.88–0.98 of the summed 1 bp signal (median across 289–305 tracks per SNP). Averages would give 0.008. → For histone: sum the 128 bp bins over 2001 bp, pseudocount 1.
- **Direction is stable wherever the effect is real.** For all 4 SNPs, the 100 most affected tracks keep their sign at 101, 501 and 2001 bp (100/100). All sign changes fall in tracks ranked below 100, with scores near zero; rs1198588 had the most (142/199 and 148/199 agree). That fully explains its 0.81 / 0.83 sign agreement in the summary.
- **rs1198588's brain tracks:**
  - The 10 largest (globus pallidus, putamen, pons, caudate, Ammon's horn, posterior cingulate, DLPFC, spinal cord, frontal cortex, spinal-cord astrocytes) are positive at all three windows.
  - Sign changes occur only in tracks with |score| ≤ 0.02 (fetal brain, cerebellum, astrocyte cultures).
  - Globus pallidus, putamen and pons are near full size already at 101 bp, while posterior cingulate and DLPFC are ≈0.05 at 101 bp and ≈0.33–0.35 at 501 bp. The effect sits at different distances from the variant in different regions.
- **rs2238057:** all 12 brain tissues are positive at every window; astrocyte cultures are near zero.
- **Why rs1198588's ratio grows at 2001 bp:** its brain scores *shrink* (putamen 0.51 → 0.42), but the background shrinks faster. Background effects are more confined to the variant, while the brain effect is spread wider. This doesn't change the decision to keep 501 bp.
- Control rs72730413: strong but very local effects (posterior cingulate/DLPFC/caudate −0.38 to −0.42 at 101 bp, ≈ −0.05 at 2001 bp).

### Status
✅ Score definition verified and documented; Session 1 naming error corrected.
✅ Sensitivity analysis run; conclusions robust except rs2238057 (brain-leaning).
✅ 128 bp units known (sums); sign changes confined to near-zero tracks.
⬜ Optional: look at `profile_rs1198588.png` alongside the 101 bp vs 501 bp cortex pattern.

---

## Session 16 (E1) — November 9–17, 2026
**Goal:** Act on Choo's Nov 9 feedback:
- normalise scores before comparing tracks, starting with a control GWAS (type 1 diabetes)
- report effect size **and** effect direction per tissue
- group tracks by organ system, averaging within a category before comparing categories

Preliminary results by the next meeting use only the 10 pilots: the 173-lead run needs a GPU.

### Nov 9 meeting with Choo: decisions
- **Raw scores can't be compared between tracks** (AlphaGenome itself normalises against common variants per track).
- **Two normalisations discussed:**
  - **SNP-by-SNP matching** (more rigorous): 1–5 control SNPs per SCZ SNP, matched on allele frequency in a European reference and on distance to the nearest gene. Each SCZ score is divided by the mean of its controls' scores in the same track.
  - **Control GWAS** (simpler, recommended starting point): divide SCZ effects by the summed or averaged effects of a control trait's variants in each track.
  - **Plan:** start with the control GWAS, upgrade to matching later.
- **Outputs per SNP:** effect size and effect direction per tissue. A variant can have large effects in many tissues but in opposite directions.
- **Grouping:** organ-system categories replace "brain vs everything else". Flat averages are skewed by over-represented tissues: ATAC is 92/167 LCL tracks.

### State of the repo before this session
- The Nov 9 scripts were written but never committed: `annotate_snps.py`, `assign_tissue_categories.py`, `tissue_profiles.py`. They were found on the Mac Nov 17 and committed (`f3f76c6`).
- `run_ism_pilot.py` was briefly moved into `phase1/` and back (`31e68ff`, `6bfcf04`); it is byte-identical to the October version.
- The Hoffman2 Slurm migration deleted the July conda env from scratch.
  - New env: `/u/project/cluo/aryasath/.conda_envs/alphagenome_env` (alphagenome-pytorch 0.3.1, torch 2.15.0, numpy 2.5.4).
  - GPU: only an old P4 card was granted (unsupported), so all model runs this week were on the Mac.

### Choosing the control GWAS
**Trait (Choo):** type 1 diabetes. The intent is a trait with no expected brain mechanism, so its variants' effects in brain tracks represent "ordinary" variant effects there.

**Study:** Chiou et al. 2021, *Nature* 594:398–402, GWAS Catalog **GCST90014023**.
- 18,942 cases / 501,638 controls, European, 9 cohorts (incl. UK Biobank, FinnGen, T1DGC), TOPMed-imputed.
- **Why this one:**
  - genome-wide (Robertson 2021 is larger but ImmunoChip-based, so dense only in immune regions)
  - European, like our PGC3 file
  - the most cases of the genome-wide T1D GWASs found (Crouch 2021 preprint 15,573; Forgetta 2020 9,358)
  - full summary statistics are public
- **Known caveat from the start:** T1D risk variants are enriched in T-cell regulatory DNA (and acinar, ductal and beta cells), so the denominator might be larger in immune tracks.

### `phase1/select_leads.py`: lead selection for any GWAS
- **Same rules as `get_lead_snps_pgc3.py`:**
  - p < 5×10⁻⁸
  - single-base A/C/G/T alleles
  - autosomes
  - MHC excluded
  - greedy 1 Mb clumping
- **Additions:**
  - a column map, so any GWAS file works
  - an optional frequency filter (≥1%)
  - a seeded random sample
  - a log (input md5, count per filter, seed)
- **Gate before use:** on the PGC3 file (md5 `6ebe2376f5cda972d37efa0f214c4df0`) it reproduced the **173 leads exactly**, with and without the frequency filter (FCON). So the filter removes no SCZ lead, and applying it to T1D keeps the two traits on the same rules.
  - The only apparent mismatch was ID formatting of the one lead without an rsID (`8:4180090_T_A`). The script now keeps the file's own IDs.

### T1D download and build check
- **Download:** downloaded on the Hoffman2 data-transfer node (`dtn.hoffman2.idre.ucla.edu`; 5 GB, too slow for home internet). Stored outside the repo at `/u/project/cluo/aryasath/gwas_data/t1d_chiou2021/`. Uncompressed md5 `ea8d72e8261f10d70e62242c9dc8ae3c` = the GWAS Catalog's.
- **The deposited file is GRCh38, not hg19.** The paper describes its own analysis as hg19, but rs2007044 is at chr12:**2,235,794** (hg38; hg19 is 2,344,960), matching the file name and yaml.
- **Consequences:**
  - MHC limits were lifted with our chain file (`hg19ToHg38.over.chain.gz`, pyliftover) from hg19 chr6:25,000,000–34,000,000 to **hg38 chr6:24,999,772–34,032,223**, so both traits exclude the same region.
  - The T1D leads skip `liftover_to_hg38.py`.
  - `run_ism_pilot.py` checks every hg38 base against A1/A2, which also confirms the build.

### T1D lead selection (Hoffman2 CPU node)

| Filter | Count |
|---|---|
| rows | 62,116,689 |
| p < 5×10⁻⁸ | 46,732 |
| not single-base A/C/G/T | 3,130 |
| frequency < 1% (incl. the paper's rare large-effect variants) | 3,069 |
| in MHC (hg38 limits) | 33,601 (72% of significant SNPs) |
| candidates | 6,932 |
| **independent leads** | **79** (A1 risk-raising 33, protective 46) |

- **Random 20 for the preliminary run (seed 2026):** rs10751776, rs1493696, rs3024493, rs55893453, rs12464462, rs4490209, rs777029891, rs9385401, rs1947178, rs61839660, rs11066320, rs238265, rs1350275, rs55993634, rs9891059, rs61759532, rs35327136, rs1808094, rs34536443, rs601338.
- **Files:**
  - `phase1/data/t1d_chiou2021_lead_snps_hg38.tsv`
  - `t1d_chiou2021_control20_hg38.tsv`
  - `t1d_chiou2021_lead_snps_hg38.log.txt`

### ISM on the 20 T1D SNPs (Mac, same script and version as the pilots)
- 20/20 scored, 12,800 track scores (`results/ism/t1d_control20_ism_scores.tsv`).
- hg38 reference = A2 (the non-effect allele) at all 20, as expected for GWAS-SSF files.
- 3 palindromic SNPs: rs4490209, rs55993634, rs34536443.

### Tissue categories: `phase2/assign_tissue_categories.py`
- **Moved from `phase1/` to `phase2/`.** From `phase1/` it crashed (`No module named 'build_track_classes'`).
- **Method:** ontology system ancestors (Uberon 2026-10-01, CL 2026-06-08) plus reviewed overrides; reuses the D1 neural groups.
- **Result:** 5,563 real tracks → 22 categories (`phase2/metadata/track_categories.tsv`).
  - DNase: 198 main-analysis tracks (tissue, primary cell, stem-cell-derived) + 107 cell lines. ATAC: 62 + 105.
  - CNS = the 16 reviewed brain DNase tracks (2 prenatal).
  - Brain vasculature and pericytes → vascular.
  - Cell lines get a category but are excluded from the main comparison.
- **Decision (Aryana, Nov 17):** the 4 prenatal limb DNase tracks (301–304) → `embryonic_extraembryonic`. Whole embryonic limbs are mixed tissue; the Uberon limb terms have no organ-system anchor. DNase now has no unassigned tissue tracks.
- **Gap vs the D1 rule:** the script reports unassigned samples instead of stopping. Read the "other" list every run.

### Control-GWAS normalisation: `phase2/tissue_profiles.py --method ratio`
- **Moved to `phase2/`.** The Nov 5 quantile method is kept (byte-identical output) and is still the default.
- **Steps (DNase, 198 main-analysis tracks, padding dropped by label):**
  1. **Denominator per track:** D_t = mean |score| of the 20 T1D SNPs in track t.
  2. **Normalised score:** n(i, t) = risk-oriented score(i, t) ÷ D_t.
     - |n| = 1 means "as large as a typical T1D lead in this track"
     - sign: + = the SCZ risk allele opens chromatin
  3. **Per SNP × category** (categories with ≥2 tracks, weighted equally):
     - size = mean |n|
     - direction = mean signed n
     - concordance = share of tracks with the majority sign
  4. **Summary:**
     - brain size (CNS)
     - body mean (13 non-neural body categories)
     - top body category (≥3 tracks)
     - specificity = brain − body mean
     - class at size 1 (brain-specific / shared / other-tissue / none)
- **Guards:**
  - 1,000 bootstrap resamples of the 20 T1D SNPs (95% intervals for brain size and specificity)
  - tracks with D_t in the bottom 5% flagged; second run with them dropped (`results/profiles/drop_low_denominator/`)
  - palindromic SNPs' direction labelled strand-assumed
- **Built-in null check:** the 10 random non-significant SCZ SNPs (October controls) go through the same normalisation.
- **Tested before use** with the random SNPs as a stand-in denominator. The figures were rendered and fixed: capped colour scale, log axes for ratios, per-panel close-up scales.

### Status
✅ Control GWAS chosen, verified (md5, build) and selected with the SCZ rules (gate passed).
✅ 20 T1D SNPs scored; tissue categories built; ratio normalisation run (all tracks and low-denominator tracks dropped).
➡️ Results and evaluation: Session 17.

---

## Session 17 (E2) — November 17, 2026: what the control-GWAS normalisation shows
**Goal:** Read the T1D-normalised profiles of the 10 pilots and decide whether control-GWAS normalisation should carry the pipeline forward.

All numbers: DNase, commit `5ca2b7c`, outputs in `results/profiles/`. n = 10 pilots, 20 T1D SNPs, 10 random SCZ SNPs: **descriptive, not a test.**

### Figures: what each shows
1. **`dnase_ratio_heatmap.png`**
   - **Layout:** one row per pilot (most brain-specific at top); one column per tissue category, neural categories first (white rule after CNS and after the neural block).
   - **Colour:** the category's direction (mean signed ratio; blue = risk allele opens, orange = closes), capped at ±4.2 (95th percentile).
   - **Number in each cell:** size (mean |ratio|).
   - **Reads:**
     - rs1198588 is the only row whose CNS cell is strongly blue while body cells are pale, and its two neighbouring neural cells (neural stem-cell-derived 4.2, retina 4.5) are strongly **orange**
     - rs2238057's largest cell is heart (3.9), not CNS (1.5)
     - rs2710323's row is orange almost everywhere, with immune 17.6 the largest value in the figure
     - rows 4–9 are pale (≤ ~2.8 anywhere)
2. **`dnase_ratio_slope.png`**
   - **Layout:** one line per pilot from the body mean (left) to brain (right), log scale, dashed line at 1. Lines are coloured when brain ≥ 1.5× the body mean (blue) or ≤ 1/1.5 (orange).
   - **Reads:** **1 pilot higher in brain (rs1198588: 0.67 → 3.65), 3 lower** (rs2710323 5.79 → 2.51, rs12129573 and rs13107325 1.34 → ~0.5). The other six change little.
   - *Fix to make:* the title is clipped at the right edge.
3. **`dnase_ratio_scatter.png`**
   - **Layout:** brain size (y) against the strongest non-neural body category (x), log–log, lines at 1 and the brain = body diagonal. Pilots coloured by class; the 10 random SCZ SNPs as triangles.
   - **Reads:**
     - rs1198588 is the only point well above the diagonal
     - **the triangles sit among the pilots**, and two of them are above 1 in brain
   - This is the key figure for calibration (below).
4. **`dnase_ratio_closeups.png`**
   - **Layout:** every track's signed ratio, grouped by category, for the 4 most and 4 least brain-specific pilots (own y-scale per panel; bar = category mean).
   - **Reads:**
     - rs1198588's brain values are a spread, not one population: globus pallidus +10.3, pons +9.1, putamen +7.4, Ammon's horn +6.0, caudate +5.6, with a few CNS tracks near 0 (fetal brain, cerebellum, astrocyte cultures, as in D2/D3)
     - its neural-progenitor and RPE tracks are uniformly negative (−4.8 to −6.1)
     - rs11191580 is negative in almost every track of every category
     - rs2710323's immune tracks reach about −35
5. **`dnase_ratio_denominators.png`**
   - **Layout:** each track's T1D denominator D_t (log scale) by category, median bar, dashed line = median of all tracks, flagged bottom-5% tracks orange.
   - **Reads:**
     - immune tracks sit as a tight cluster at the top (median 1.57× all tracks)
     - CNS is slightly above the median (1.14×; 1 flagged track)
     - liver, vascular and PNS are lowest
     - the spread within categories is about 3-fold

The `drop_low_denominator/` copies are the same figures without the 10 flagged tracks.

### Results (all tracks; 95% bootstrap interval)

| SNP (gene) | Brain size | Body mean | Top body category | Brain − body | Brain direction (concordance) | Class |
|---|---|---|---|---|---|---|
| rs1198588 (ENSG00000297419; nearest coding gene DPYD) | **3.65** [2.33, 7.08] | 0.67 | embryonic 1.29 | **+2.97** [1.80, 6.05] | +3.60 (0.94), palindromic | shared |
| rs2238057 (CACNA1C) | 1.53 [0.99, 3.02] | 1.05 | **heart 3.92** | +0.48 [0.03, 1.64] | +1.44 (0.81) | shared |
| rs11191580 (NT5C2) | 1.50 [1.04, 2.80] | 1.30 | liver 1.98 | +0.19 [−0.14, 0.97] | −1.50 (1.00) | shared |
| rs2710323 (ITIH1) | 2.51 [1.71, 5.06] | 5.79 | **immune 17.61** | −3.28 [−6.32, −2.36] | −1.59 (0.62) | shared |
| rs4702 (FURIN) | 0.62 | 0.69 | stromal 1.93 | −0.07 | +0.59 (0.88) | other-tissue |
| rs4129585 (TSNARE1) | 0.72 | 1.01 | kidney 1.67 | −0.30 | −0.70 (0.94) | other-tissue |
| rs12129573 | 0.51 | 1.34 | heart 2.77 | −0.83 | −0.36 (0.56) | other-tissue |
| rs13107325 (SLC39A8) | 0.48 | 1.34 | liver 2.81 | −0.86 | −0.23 (0.56) | other-tissue |
| rs58120505 | 0.29 | 0.20 | stromal 0.35 | +0.09 | +0.29 (1.00) | none |
| rs778371 (SNORC) | 0.17 | 0.24 | lung 0.52 | −0.07 | +0.12 (0.69) | none |

- **Low-denominator tracks dropped** (188 tracks): every pilot keeps its class.
  - rs1198588: brain 3.85, specificity +3.20 [1.96, 6.48]
  - rs2238057's brain advantage shrinks to +0.09 [−0.22, 0.65]
- **Random SCZ SNPs, same normalisation:**
  - brain size median 0.58 (range 0.04–1.40) vs pilots 0.67
  - 2 of 10 above 1 in brain (rs62225795 1.40, rs72730413 1.19)
  - highest specificity +0.31
  - pilots above the highest random specificity: rs1198588 (+2.97) and, narrowly, rs2238057 (+0.48)

### What the results say about the pilots
1. **rs1198588 is the one robust brain-leaning signal.**
   - Brain sits 2.97 above the body mean, with an interval far from 0, in both runs.
   - It is strongest in adult basal ganglia, pons and hippocampus, and **opposite in sign** in stem-cell-derived neural progenitors and RPE. This is the D2/D3 pattern, now normalised.
   - "Brain-specific" is too strong: the neural-lineage cells respond as strongly, just in the other direction. "Adult-brain-leaning, with a developmental/neural-lineage sign flip" is the accurate description.
   - Direction assumes the + strand (palindromic).
2. **rs2238057 (CACNA1C) is not brain-specific.** Heart (3.9) and skeletal muscle (2.3) exceed brain (1.5), all in the same direction. Its brain advantage over the body mean disappears when low-denominator tracks are dropped. This is the third downgrade (D3: brain-leaning; now shared with muscle).
3. **rs11191580 (NT5C2)** closes chromatin in essentially every tissue (brain concordance 1.00): broad, not brain-specific.
4. **rs2710323 (ITIH1)** remains the strongest non-brain signal: immune −17.6, all immune tracks agreeing. This is despite the T1D denominator being *largest* in immune tracks, so the immune signal is, if anything, understated.
5. The other six pilots are small everywhere.

### What it tells us about control-GWAS normalisation
**Assumption being tested:** T1D variants have no special brain effects, so dividing by them puts every track on a common "ordinary-variant" scale. Four checks:

1. **Brain was not inflated by T1D.** CNS denominators are 1.14× the median track. But the denominator must be neutral in **every** category we compare brain against, not just in brain, because specificity is brain minus the body mean. T1D's are not: immune 1.57×, kidney 1.34×, liver 0.69× and vascular 0.73× the median track. This tilts the body side of every comparison.
2. **The denominators are not stable properties of the tracks** (checked Nov 17 on the committed data):
   - **Split-half reliability:** dividing the 20 T1D SNPs into random halves of 10 (500 splits), the two halves rank the tracks' denominators with median Spearman 0.36. The full 20-SNP set is about 0.53 (Spearman–Brown). The 10 random SNPs give a similar 0.54.
   - **The two reference sets disagree:** T1D and random-SNP denominators correlate at **−0.18** across tracks. A denominator that measured each track's own scale would agree between any two reasonable reference sets.
   - **A few SNPs dominate:** one T1D SNP (rs10751776) is the largest contributor in 50 of 198 tracks, and in **28 of 31 immune tracks**, where it supplies a median 39% of the sum. The immune excess is largely one variant's biology.
   - **The reference changes brain scaling:** with random SNPs as the reference, CNS denominators were *small* (0.72× the median) and rs1198588's brain size was 5.84. With T1D they are slightly large (1.14×) and it is 3.65. Same SNP, same scores: a 1.6-fold change from the choice of 10–20 reference SNPs.
   - **How many SNPs would be needed:** extrapolating the split-half figure, ~75 control SNPs give reliability ~0.8 and ~170 give ~0.9.
3. **Rankings, though, are robust.** Leaving each T1D SNP out in turn (20 runs) moves no pilot's class:
   - rs1198588's specificity stays +2.84 to +3.89
   - rs2710323's −4.01 to −3.13
   - rs2238057's +0.40 to +0.89

   The SNPs that stood out on raw scores (D2) stand out here. Normalisation changes sizes and thresholds, not who is at the extremes.
4. **The threshold has no null meaning.** "Size > 1 = bigger than a typical T1D lead" doesn't separate pilots from random non-significant SNPs (medians 0.67 vs 0.58; random SNPs reach 1.40). A T1D lead is not much larger than a random variant in most tracks, so "1" is an arbitrary yardstick, not a calibrated cut-off. Classes based on it (4 "shared" pilots) shouldn't be reported as findings.

**Summary:** with 20 SNPs, the control-GWAS denominator ranks SNPs consistently, but it measures "which tissues these 20 T1D variants happen to act in" more than "how loud each track is". It also carries the control trait's own tissue biology into every body comparison. Its numbers (sizes, classes) are therefore not interpretable on their own.

### Decision (proposed, to discuss with Choo): move to matched control SNPs now
**Recommendation:** switch the main reference to matched non-associated SNPs now. Keep the T1D run as a documented comparison, and keep the control-GWAS option in the pipeline.

**Why now rather than later:**
1. **Reliability needs more control SNPs anyway.** A stable denominator needs ~75–170 control variants. The full 79-lead T1D set would just reach the low end. Once that many variants must be scored, there's no compute saving left in using a control trait, and matched SNPs cost the same.
2. **Matched SNPs remove the control-trait bias.** Non-associated variants have no shared tissue biology. A control trait always has some (T1D: immune, pancreas), and the general pipeline would need a "neutral" trait for every disease, which doesn't exist.
3. **Matching gives the threshold a meaning.** A matched null answers "is this SNP's effect in brain larger than that of comparable non-associated variants?" and allows real per-category tests at the 173-lead scale. "Bigger than a typical T1D lead" can't.
4. **It's cheap at pilot scale and needs almost no new analysis code.** 10 pilots × 5 controls = 50 SNPs, about one Mac run. `tissue_profiles.py --method ratio --reference matched` works unchanged; only a selection script is new.

**What the control-GWAS step still gave us:** the normalisation code, outputs and figures; the general lead selector (reproduces the PGC3 leads exactly); the tissue categories; the evidence that rankings are robust to the denominator; and the reliability numbers that size the matched set.

### Proposed matching logic (Choo's spec: European allele frequency + distance to nearest gene)
1. **Candidate pool:** the PGC3 European file itself, which carries the European control allele frequency (FCON) for every variant on the same build as the leads. Keep variants that are:
   - not associated: PVAL > 0.5
   - single-base, autosomal
   - MAF ≥ 0.01
   - outside the MHC
   - more than 1 Mb from every one of the 173 leads, so they're not in LD with a lead
   - well imputed: IMPINFO ≥ 0.9
2. **Distance to nearest gene:** distance from each variant to the nearest protein-coding gene TSS, computed locally from a GENCODE GRCh37 annotation (millions of variants is too many for the Ensembl REST calls `annotate_snps.py` makes). Leads and candidates use the same annotation and build.
3. **Matching per lead:**
   - folded MAF within ±0.02, and log10(distance + 1) within ±0.1 (about ±25%)
   - if fewer than k candidates qualify, widen both tolerances step by step and record the widening
   - draw **k = 5** at random (fixed seed), without replacement across leads
   - require matched controls to be >1 Mb from each other, so they're independent
   - record each lead's matched set (`lead_rsid` column) and its match quality (ΔMAF, Δlog-distance)
4. **Check match quality before any model run:** compare the lead and matched distributions of MAF and distance (summary table + QQ plot). They should overlap; if not, tighten the match.
5. **Score:** liftover → ISM (same script) → `score_tracks.py` with a `matched` set.
6. **Normalise two ways:**
   - **Main: pooled matched denominator.** All matched controls of all leads form one reference set per track (`--reference matched`). With 5 per lead that's 50 SNPs for the pilots (reliability ≈ 0.75) and 865 for the 173 leads (≈ 0.98).
   - **Sensitivity: Choo's per-SNP ratio.** Each lead's score ÷ the mean |score| of its own 5 controls. This needs a small new option. Five SNPs give reliability ≈ 0.26, so per-SNP ratios will be noisy and should be read alongside the pooled version.
7. **Tests at 173-lead scale:** leads vs their matched controls per tissue category (paired by matched set; permutation or Wilcoxon, corrected across categories).

**Pilot-scale cost now:** 50 SNPs on the Mac. **Full scale (later, GPU):** 173 + 865 SNPs.

### Status
✅ T1D-normalised pilot profiles; robustness (bootstrap, leave-one-out, low-denominator drop); reliability of the denominator measured.
✅ Conclusion: rankings robust (rs1198588 adult-brain-leaning; rs2710323 immune; rs2238057 shared with muscle); sizes and classes not interpretable with a 20-SNP control-trait denominator.
⬜ Choo's decision on moving to matched controls now.
⬜ If yes: `phase1/select_matched_controls.py` → match-quality check → ISM on 50 matched SNPs → pooled + per-SNP ratio.
⬜ Fix the clipped slope-plot title.