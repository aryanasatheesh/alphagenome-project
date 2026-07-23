"""
Phase 2 (steps 2a/2b): Extract hg38 windows around SCZ pilot SNPs and
create reference + alt (SNP-substituted) sequences.

Output: a .npz per SNP containing ref_seq and alt_seq as strings, plus metadata.
The AlphaGenome prediction step (2c) is a SEPARATE script we'll write next.
"""

import pysam
import pandas as pd
import numpy as np
import os

# ---- config ----
FASTA = "/u/project/cluo/aryasath/hg38.fa"
PILOT = "/u/project/cluo/aryasath/alphagenome_project/phase2/scz_pilot_5snps.tsv"
OUTDIR = "/u/project/cluo/aryasath/alphagenome_project/phase2/outputs"
WINDOW = 131072  # AlphaGenome input size
HALF = WINDOW // 2  # 65536

os.makedirs(OUTDIR, exist_ok=True)

# ---- load ----
genome = pysam.FastaFile(FASTA)
snps = pd.read_csv(PILOT, sep="\t")
print(f"Loaded {len(snps)} pilot SNPs\n")

for _, row in snps.iterrows():
    rsid = row["rsid"]
    chrom = f"chr{row['chr']}"      # GWAS uses '12', FASTA uses 'chr12'
    pos_1based = int(row["pos"])    # GWAS coords are 1-based
    ref_allele = row["ref"]
    alt_allele = row["alt"]

    # --- coordinate conversion ---
    # pysam.fetch is 0-based, half-open. A 1-based position P sits at
    # 0-based index P-1. We center the window on the SNP.
    pos_0based = pos_1based - 1
    start = pos_0based - HALF
    end = start + WINDOW            # half-open, so this yields exactly WINDOW bases

    if start < 0:
        print(f"  ! {rsid}: window runs off start of {chrom}, skipping")
        continue

    # --- extract reference window ---
    ref_seq = genome.fetch(chrom, start, end).upper()
    assert len(ref_seq) == WINDOW, f"{rsid}: got {len(ref_seq)} bp, expected {WINDOW}"

    # --- verify the base at SNP position matches expected ref allele ---
    snp_idx_in_window = pos_0based - start   # should be exactly HALF
    actual_base = ref_seq[snp_idx_in_window]

    if actual_base != ref_allele.upper():
        print(f"  ! {rsid}: ref mismatch — FASTA has '{actual_base}', "
              f"GWAS ref is '{ref_allele}'. Flipping ref/alt.")
        # GWAS A1/A2 assignment can be flipped relative to the + strand;
        # if FASTA matches the alt allele, swap so we substitute correctly.
        if actual_base == alt_allele.upper():
            ref_allele, alt_allele = alt_allele, ref_allele
        else:
            print(f"  !! {rsid}: FASTA base matches NEITHER allele — skipping")
            continue

    # --- build alt sequence ---
    alt_list = list(ref_seq)
    alt_list[snp_idx_in_window] = alt_allele.upper()
    alt_seq = "".join(alt_list)

    # sanity: alt and ref differ at exactly one position
    diffs = sum(1 for a, b in zip(ref_seq, alt_seq) if a != b)
    assert diffs == 1, f"{rsid}: expected 1 diff, got {diffs}"

    # --- save ---
    out = os.path.join(OUTDIR, f"{rsid}_sequences.npz")
    np.savez(
        out,
        rsid=rsid,
        chrom=chrom,
        pos_1based=pos_1based,
        window_start_0based=start,
        window_end=end,
        snp_idx_in_window=snp_idx_in_window,
        ref_allele=ref_allele,
        alt_allele=alt_allele,
        ref_seq=ref_seq,
        alt_seq=alt_seq,
    )
    print(f"  ✓ {rsid} ({chrom}:{pos_1based}) {ref_allele}>{alt_allele} "
          f"— window {start}-{end}, saved {os.path.basename(out)}")

genome.close()
print("\nDone. Sequence files written to:", OUTDIR)
