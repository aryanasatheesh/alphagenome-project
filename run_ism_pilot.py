"""
run_ism_pilot.py
----------------
Runs in silico mutagenesis (ISM) on pilot SCZ SNPs using AlphaGenome.

For each SNP, this script:
1. Fetches the 131,072 bp hg38 sequence window centered on the variant
2. Verifies which GWAS allele matches the hg38 reference base, and assigns
   the true reference and alternate alleles accordingly
3. One-hot encodes both ref and alt sequences
4. Runs both through AlphaGenome (pretrained model)
5. Applies a 501 bp spatial mask centered on the variant
6. Computes variant effect scores: log2[(sum(alt)+1) / (sum(ref)+1)] per track
7. Saves per-track scores and summary statistics

IMPORTANT NOTE ON REF/ALT ASSIGNMENT:
    In GWAS summary statistics, A1 (effect allele) and A2 (other allele) do NOT
    correspond to "alternate" and "reference" relative to the genome assembly.
    A1 is simply the allele whose effect is measured by BETA — it can be either
    the reference genome allele or the true alternate. This script checks the
    actual hg38 base at each SNP position and assigns ref/alt correctly:
    - ref_allele = whichever of A1/A2 matches hg38
    - alt_allele = the other one
    This was discovered when all 10 PGC3 pilot SNPs showed A1 (labeled "alt"
    in our pipeline) matching the hg38 reference base.

Input:
    scz_pgc3_pilot_10snps_hg38.tsv — pilot SNPs with hg38 coordinates

Output:
    results/ism/pilot_ism_scores.tsv — variant effect scores per track per SNP
    results/ism/pilot_ism_summary.tsv — summary: max effect, brain vs non-brain

Usage:
    # on Mac (local development, small runs):
    PYTORCH_ENABLE_MPS_FALLBACK=1 python run_ism_pilot.py

    # on Hoffman2 (batch runs):
    # submit as GPU job (see job script template)

Requires:
    - alphagenome-pytorch installed
    - Model weights at ~/alphagenome-weights/model_fold_0.safetensors
    - Internet access (to fetch sequences from Ensembl REST API)
    - For Hoffman2: hg38.fa available locally (avoids API rate limits for batch)
"""

import torch
import numpy as np
import csv
import os
import json
import urllib.request
import time
import sys

from alphagenome_pytorch import AlphaGenome

# === Configuration ===
PILOT_FILE = os.path.join('phase1', 'data', 'scz_pgc3_pilot_10snps_hg38.tsv')
WEIGHTS_PATH = os.path.expanduser('~/alphagenome-weights/model_fold_0.safetensors')
OUTPUT_DIR = os.path.join('results', 'ism')
SEQ_LEN = 131_072          # AlphaGenome input length
MASK_WINDOW = 501           # bp around variant for scoring (recommended for chromatin)
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
INDEX_TO_BASE = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}


def fetch_hg38_sequence(chrom, start, end):
    """Fetch sequence from Ensembl REST API (hg38/GRCh38).
    
    Coordinates are 1-based inclusive (Ensembl convention).
    Returns uppercase DNA string.
    """
    url = (f'https://rest.ensembl.org/sequence/region/human/'
           f'{chrom}:{start}..{end}?content-type=application/json')
    req = urllib.request.Request(url, headers={'Content-Type': 'application/json'})
    
    max_retries = 3
    for attempt in range(max_retries):
        try:
            resp = urllib.request.urlopen(req, timeout=60)
            data = json.loads(resp.read())
            return data['seq'].upper()
        except Exception as e:
            if attempt < max_retries - 1:
                print(f'    Retry {attempt+1}/{max_retries}: {e}')
                time.sleep(2)
            else:
                raise


def one_hot_encode(seq_str):
    """Convert DNA string to one-hot tensor [1, seq_len, 4].
    
    Encoding: A=0, C=1, G=2, T=3
    N or unknown bases get all-zero encoding.
    """
    indices = np.array([BASE_TO_INDEX.get(b, -1) for b in seq_str])
    one_hot = np.zeros((len(seq_str), 4), dtype=np.float32)
    valid = indices >= 0
    one_hot[valid, indices[valid]] = 1.0
    return torch.tensor(one_hot).unsqueeze(0)  # [1, seq_len, 4]


def compute_variant_effect(ref_pred, alt_pred, center_pos, mask_window=MASK_WINDOW):
    """Compute log2 fold-change variant effect score per track.
    
    Applies a spatial mask of `mask_window` bp centered on the variant,
    then computes: log2[(sum(alt) + 1) / (sum(ref) + 1)]
    
    This is the recommended aggregation for chromatin accessibility
    (ATAC/DNase) from the AlphaGenome variant scoring documentation.
    The +1 pseudocount avoids log(0) for tracks with no signal.
    
    Args:
        ref_pred: tensor [1, seq_len, n_tracks] — ref predictions
        alt_pred: tensor [1, seq_len, n_tracks] — alt predictions  
        center_pos: int — position of variant in the sequence (0-indexed)
        mask_window: int — width of spatial mask in bp
    
    Returns:
        scores: numpy array [n_tracks] — log2 fold-change per track
    """
    half = mask_window // 2
    start = max(0, center_pos - half)
    end = min(ref_pred.shape[1], center_pos + half + 1)
    
    ref_masked = ref_pred[0, start:end, :].cpu().numpy()
    alt_masked = alt_pred[0, start:end, :].cpu().numpy()
    
    ref_sum = ref_masked.sum(axis=0) + 1  # pseudocount
    alt_sum = alt_masked.sum(axis=0) + 1
    
    scores = np.log2(alt_sum / ref_sum)
    return scores


def main():
    # --- Setup ---
    device = 'mps' if torch.backends.mps.is_available() else (
        'cuda' if torch.cuda.is_available() else 'cpu')
    print(f'Device: {device}')
    
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # --- Load model ---
    print(f'Loading AlphaGenome from {WEIGHTS_PATH}...')
    model = AlphaGenome.from_pretrained(WEIGHTS_PATH, device=device)
    model.eval()
    print('Model loaded.')
    
    # --- Load pilot SNPs ---
    with open(PILOT_FILE) as f:
        snps = list(csv.DictReader(f, delimiter='\t'))
    print(f'Loaded {len(snps)} pilot SNPs.')
    
    # --- Process each SNP ---
    all_scores = []
    center = SEQ_LEN // 2  # variant position in the extracted window (0-indexed)
    
    for i, snp in enumerate(snps):
        rsid = snp['rsid']
        chrom = snp['chr']
        pos = int(snp['pos'])        # 1-based hg38 position
        gwas_ref = snp['ref']        # A2 from GWAS (may not match hg38!)
        gwas_alt = snp['alt']        # A1 from GWAS (may not match hg38!)
        
        print(f'\n--- [{i+1}/{len(snps)}] {rsid} chr{chrom}:{pos} ---')
        
        # Step 1: Fetch 131,072 bp hg38 window centered on SNP
        # 1-based coordinates for Ensembl
        window_start = pos - center
        window_end = pos + center - 1  # inclusive
        print(f'  Fetching {SEQ_LEN:,} bp window (chr{chrom}:{window_start}-{window_end})...')
        seq = fetch_hg38_sequence(chrom, window_start, window_end)
        
        if len(seq) != SEQ_LEN:
            print(f'  ERROR: expected {SEQ_LEN} bp, got {len(seq)} bp. Skipping.')
            continue
        
        # Step 2: Verify ref/alt against hg38 reference
        hg38_base = seq[center]
        
        if hg38_base == gwas_ref:
            true_ref = gwas_ref
            true_alt = gwas_alt
            swapped = False
        elif hg38_base == gwas_alt:
            true_ref = gwas_alt
            true_alt = gwas_ref
            swapped = True
        else:
            print(f'  ERROR: hg38 base ({hg38_base}) matches neither '
                  f'GWAS ref ({gwas_ref}) nor alt ({gwas_alt}). Skipping.')
            continue
        
        if swapped:
            print(f'  NOTE: GWAS A1/A2 swapped — hg38 ref is {true_ref} '
                  f'(was labeled as GWAS alt/A1)')
        print(f'  Ref (hg38): {true_ref}, Alt (mutated): {true_alt}')
        
        # Step 3: Build ref and alt sequences
        ref_seq = seq  # already contains hg38 reference
        alt_seq = seq[:center] + true_alt + seq[center+1:]
        
        # Sanity check
        assert ref_seq[center] == true_ref
        assert alt_seq[center] == true_alt
        assert len(ref_seq) == len(alt_seq) == SEQ_LEN
        
        # Step 4: One-hot encode
        ref_tensor = one_hot_encode(ref_seq).to(device)
        alt_tensor = one_hot_encode(alt_seq).to(device)
        
        # Step 5: Run AlphaGenome
        print(f'  Running forward pass (ref)...')
        with torch.inference_mode():
            ref_out = model.predict(ref_tensor, organism_index=0)
            print(f'  Running forward pass (alt)...')
            alt_out = model.predict(alt_tensor, organism_index=0)
        
        # Step 6: Compute variant effect scores for ATAC and DNase
        for modality in ['atac', 'dnase']:
            ref_pred = ref_out[modality][1]  # [1] = predictions tensor
            alt_pred = alt_out[modality][1]
            n_tracks = ref_pred.shape[-1]
            
            scores = compute_variant_effect(ref_pred, alt_pred, center)
            
            max_abs = np.max(np.abs(scores))
            max_idx = np.argmax(np.abs(scores))
            mean_abs = np.mean(np.abs(scores))
            
            print(f'  {modality.upper()} ({n_tracks} tracks): '
                  f'max |score| = {max_abs:.4f} (track {max_idx}), '
                  f'mean |score| = {mean_abs:.4f}')
            
            for t in range(n_tracks):
                all_scores.append({
                    'rsid': rsid,
                    'chr': chrom,
                    'pos': pos,
                    'true_ref': true_ref,
                    'true_alt': true_alt,
                    'gwas_alleles_swapped': swapped,
                    'modality': modality,
                    'track_index': t,
                    'variant_effect_score': f'{scores[t]:.6f}',
                })
        
        # Free GPU memory between SNPs
        del ref_tensor, alt_tensor, ref_out, alt_out
        if device == 'mps':
            torch.mps.empty_cache()
        elif device == 'cuda':
            torch.cuda.empty_cache()
    
    # --- Save results ---
    scores_out = os.path.join(OUTPUT_DIR, 'pilot_ism_scores.tsv')
    fieldnames = ['rsid', 'chr', 'pos', 'true_ref', 'true_alt',
                  'gwas_alleles_swapped', 'modality', 'track_index',
                  'variant_effect_score']
    with open(scores_out, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(all_scores)
    
    print(f'\n=== DONE ===')
    print(f'Saved {len(all_scores)} scores to {scores_out}')
    print(f'({len(snps)} SNPs × tracks across ATAC + DNase)')


if __name__ == '__main__':
    main()