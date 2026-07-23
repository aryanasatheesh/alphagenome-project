# AlphaGenome Variant Effect Project

**Aryana Satheesh · Luo Lab (PI: Prof. Chongyuan Luo) · UCLA**
**Mentor: Cuining (Choo) Liu**

---

## Project overview

This project uses [AlphaGenome](https://github.com/genomicsxai/alphagenome-pytorch) (PyTorch port) to predict whether non-coding genetic variants associated with schizophrenia (SCZ) and autism spectrum disorder (ASD) alter chromatin accessibility in specific fetal brain cell types.

The core approach is **in silico mutagenesis (ISM)**: for each GWAS risk variant, we extract a ~131 kb window of hg38 sequence centered on the SNP, run it through AlphaGenome with the reference allele, swap in the alternate allele and run again, then compare the predicted chromatin accessibility signals between the two. The difference — the delta track — tells us whether and where the variant is predicted to open or close chromatin, and in which cell types.

### Why this matters
Most GWAS hits fall in non-coding regions, making it hard to go from a statistical association to a mechanistic explanation. Sequence-to-function models like AlphaGenome let us ask computationally which variants are likely to be regulatory, in which tissues, without requiring experimental validation for every SNP.

---

## Repository structure

```
alphagenome_project/
├── README.md                         # this file
├── notebooks/
│   └── lab_notebook.md               # chronological session log
├── docs/
│   └── methods_draft.md              # thesis-ready methods prose (built from notebook)
├── phase1/                           # GWAS processing and liftover
│   ├── get_lead_snps.py
│   ├── extract_sequences.py
│   ├── liftover_to_hg38.py
│   └── data/
│       ├── scz_lead_snps.tsv         # 141 lead SNPs, hg19 (original)
│       ├── scz_lead_snps_hg38.tsv    # 141 lead SNPs, hg38 (lifted)
│       ├── scz_pilot_5snps.tsv       # 5 pilot SNPs, hg19 (original)
│       └── scz_pilot_5snps_hg38.tsv  # 5 pilot SNPs, hg38 (lifted)
└── phase2/                           # ISM inference and visualization
    ├── run_ism.py                     # (in progress)
    └── visualize_ism.py              # (in progress)
```

---

## Pipeline overview

### Phase 1 — GWAS processing (complete)
1. Download scz2018clozuk GWAS summary statistics from PGC
2. Filter to genome-wide significant SNPs (p < 5e-8), remove indels
3. Apply 1 Mb clumping window per chromosome to get independent loci; exclude MHC (chr6)
4. Output: 141 independent lead SNPs → `scz_lead_snps_hg38.tsv`
5. Select 5 pilot SNPs for pipeline development → `scz_pilot_5snps_hg38.tsv`
6. Liftover hg19 → hg38 (confirmed: GWAS uses hg19 coordinates)

### Phase 2 — ISM with pretrained model (in progress)
1. For each SNP: extract 131,072 bp hg38 sequence window centered on variant
2. Build ref and alt copies (swap single base at center)
3. Run both through AlphaGenome pretrained model
4. Apply spatial mask (501 bp window for ATAC/DNase) and compute log2-ratio variant effect score
5. Compare brain vs non-brain tracks; visualize delta tracks

### Phase 3 — Fine-tuning (future)
- Fine-tune AlphaGenome on Ziffra et al. fetal brain ATAC-seq using LoRA
- Re-run ISM with fine-tuned model
- Compare pretrained vs fine-tuned predictions across brain cell types

---

## Environment setup

### Mac (local development)
```bash
# requires Python 3.12+ (install via homebrew if needed)
brew install python@3.12
python3.12 -m venv ~/alphagenome-env
source ~/alphagenome-env/bin/activate
pip install alphagenome-pytorch

# download model weights (~921 MB, one time)
pip install huggingface_hub
hf download gtca/alphagenome_pytorch model_fold_0.safetensors --local-dir ~/alphagenome-weights
```

Run scripts with MPS fallback enabled:
```bash
PYTORCH_ENABLE_MPS_FALLBACK=1 python3 your_script.py
```

### Hoffman2 (scale compute)
```bash
module load anaconda3
conda activate alphagenome_env   # Python 3.12, created July 2026
pip install alphagenome-pytorch  # reinstall if environment was reset
```

Project files live at: `/u/project/cluo/aryasath/alphagenome_project/`
hg38 FASTA: `/u/project/cluo/aryasath/hg38.fa`

---

## Key data notes

- **GWAS:** scz2018clozuk (CLOZUK + PGC2 SCZ meta-analysis, 2018)
  - 8,064,800 total SNPs; 18,088 genome-wide significant (p < 5e-8)
  - Coordinates are **hg19** — must liftover before use with AlphaGenome
  - A1 = effect/alt allele, A2 = other/ref allele
- **AlphaGenome input:** 131,072 bp one-hot encoded sequence, hg38
- **AlphaGenome output:** predictions at every bp for 11 modalities; ATAC shape `[1, 131072, 256]`
- **Brain tracks in pretrained model:** ATAC=1, DNase=21, CHIP_HISTONE=100 (out of 335 total brain tracks across all modalities)
- **Track metadata:** `google-deepmind/alphagenome_research` → `OutputMetadataResponse_ORGANISM_HOMO_SAPIENS.textproto`

---

## References
- Ziffra et al. 2021 — fetal brain snATAC-seq; chromatin diverges at IPC stage; SCZ enrichment in excitatory/inhibitory neuron enhancers
- AlphaGenome paper (Google DeepMind)
