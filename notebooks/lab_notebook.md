# AlphaGenome Project — Lab Notebook
**Researcher:** Aryana Satheesh  
**PI:** Prof. Chongyuan Luo  
**Mentor:** Cuining (Choo) Liu

This notebook documents every work session in chronological order: what was run, what the output was, decisions made, and any blockers. After each completed phase, these notes get converted into thesis-ready methods prose.

## Session 3 - August 3, 2026

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