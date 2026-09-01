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
Options like `DIFF_SUM_LOG2` are read right to left (innermost operation first):
- `DIFF_SUM_LOG2` → apply log2 first, then sum, then compute ALT − REF
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
1. **The 6,400 ISM scores from Session 8 are valid and do not need rerunning** — they were computed at 1 bp resolution on human predictions (see note below). They do need **filtering**: drop `track_index >= 167` for ATAC and `>= 305` for DNase before any summary statistic. All max/mean values in the Session 8 table were computed over full width and are therefore diluted.
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

## Session 10 — August 25, 2026
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

