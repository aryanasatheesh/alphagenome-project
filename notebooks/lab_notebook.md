# AlphaGenome Project — Lab Notebook
**Researcher:** Aryana Satheesh  
**PI:** Prof. Chongyuan Luo  
**Mentor:** Cuining (Choo) Liu

## Project Goal
Fine-tune AlphaGenome (PyTorch) on prenatal, cell type-specific fetal brain 
ATAC-seq data (ZIFFRA dataset) and apply in silico mutagenesis at ASD/SCZ 
risk variants to predict regulatory effects in specific brain cell types.

---

## Environment Setup
- **Cluster:** Hoffman2 (UCLA)
- **Conda env:** `alphagenome_env` (Python 3.12)
- **Key packages:** alphagenome-pytorch==0.3.0, torch==2.6.0+cu118
- **Weights:** `/u/project/cluo/aryasath/alphagenome/alphagenome_weights/model_all_folds.safetensors`

## Completed Steps
### [DATE] Step 0 — Installation and smoke test
- Installed alphagenome-pytorch 0.3.0 with Python 3.12
- Key fix: PyTorch cu118 required (Hoffman2 driver too old for cu121)
- Verified forward pass on chr21:14000000-14131072
- Output modalities confirmed: atac, dnase, rna_seq, cage, chip_histone, chip_tf, procap, splice sites, pair_activations

## Upcoming Steps
- [ ] Step 1: Identify brain-relevant tracks in pretrained model
- [ ] Step 2: Obtain and prepare ZIFFRA data
- [ ] Step 3: Fine-tune with LoRA on fetal brain ATAC
- [ ] Step 4: Evaluate fine-tuning (pretrained vs fine-tuned Pearson r)
- [ ] Step 5: ISM at ASD/SCZ risk variants
