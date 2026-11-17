"""
run_ism_pilot.py
----------------
In silico mutagenesis (ISM) of GWAS SNPs with pretrained AlphaGenome.

For each SNP, this script:
1. Fetches the 131,072 bp hg38 sequence window centred on the variant
2. Works out which GWAS allele (A1 or A2) is the hg38 reference base
3. Builds the reference sequence and a copy with the other allele substituted
4. Runs both through AlphaGenome
5. Scores each output track: log2[(sum(alt)+1) / (sum(ref)+1)] over a 501 bp
   window centred on the variant
6. Also writes the same score oriented by disease risk (see below)

ALLELE CONVENTIONS (read this before interpreting any score)
    GWAS alleles are statistical roles, not genome roles:
        A1 = effect allele. BETA is the log odds ratio for carrying A1.
        A2 = the other allele.
    Neither is defined as the reference base of the genome. This script reads
    the hg38 base at the variant position and records:
        hg38_ref     the base present in hg38 at this position
        hg38_alt     the other GWAS allele (the base we substitute in)
        hg38_ref_is  'A1' or 'A2': which GWAS allele hg38 carries
    In PGC3, hg38_ref_is was 'A1' for all 20 pilot and control SNPs checked by
    October 2026, so substitution usually replaces the effect allele with the
    other allele. Do not assume this; use the column.

    variant_effect_score = log2(signal with hg38_alt / signal with hg38_ref)
        Positive: substituting the alternate base opens chromatin.

    risk_oriented_score  = log2(signal with risk allele / signal with non-risk allele)
        The risk allele is A1 if beta > 0 and A2 if beta < 0.
        Positive: the schizophrenia risk allele opens chromatin.
        It equals variant_effect_score when hg38_alt is the risk allele, and
        its negative when hg38_ref is the risk allele. Left blank if beta is
        0 or missing, because then there is no risk direction.

    Palindromic SNPs (A/T or C/G) are flagged. For them, matching against hg38
    cannot reveal which DNA strand the GWAS reported, so hg38_ref_is and
    risk_oriented_score assume the GWAS used the + strand. Every non-palindromic
    PGC3 SNP checked so far matched on the + strand, which supports this, but it
    is an assumption. Magnitude (|variant_effect_score|) is unaffected.

Input (TSV with header; produced by the selection scripts + liftover_to_hg38.py)
    required columns: rsid, chr, pos (1-based hg38), a1, a2, beta
    other columns (e.g. a1_freq_controls) are carried into the variant table

Output
    results/ism/<out_name>
        one row per SNP x modality x track (padding channels included;
        downstream analysis must drop them)
    results/ism/<out_name stem>.variants.tsv
        one row per input SNP: alleles, hg38 assignment, risk allele, and
        whether it was scored or skipped (with the reason)

Usage
    PYTORCH_ENABLE_MPS_FALLBACK=1 python run_ism_pilot.py [input_tsv] [out_name]
    defaults: phase1/data/scz_pgc3_pilot_10snps_hg38.tsv  pilot_ism_scores.tsv

Exit status
    0 if every SNP was scored, 1 if any were skipped (scores for the rest are
    still written, so a long run is not lost to one bad SNP).

Requires
    alphagenome-pytorch, weights at ~/alphagenome-weights/model_fold_0.safetensors,
    internet access (Ensembl REST API)
"""

import csv
import json
import math
import os
import sys
import time
import urllib.request

import numpy as np
import torch

from alphagenome_pytorch import AlphaGenome

# === Configuration ===
INPUT_FILE = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    'phase1', 'data', 'scz_pgc3_pilot_10snps_hg38.tsv')
OUT_NAME = sys.argv[2] if len(sys.argv) > 2 else 'pilot_ism_scores.tsv'
WEIGHTS_PATH = os.path.expanduser('~/alphagenome-weights/model_fold_0.safetensors')
OUTPUT_DIR = os.path.join('results', 'ism')
SEQ_LEN = 131_072          # AlphaGenome input length
MASK_WINDOW = 501           # bp around variant for scoring (recommended for chromatin)
MODALITIES = ['atac', 'dnase']
RESOLUTION = 1              # output dict key: 1 = 1 bp resolution (128 = 128 bp bins)
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
BASES = set(BASE_TO_INDEX)
PALINDROMIC_PAIRS = ({'A', 'T'}, {'C', 'G'})
REQUIRED_COLUMNS = ['rsid', 'chr', 'pos', 'a1', 'a2', 'beta']

SCORE_FIELDS = ['rsid', 'chr', 'pos', 'a1', 'a2', 'beta',
                'hg38_ref', 'hg38_alt', 'hg38_ref_is', 'palindromic',
                'modality', 'track_index',
                'variant_effect_score', 'risk_oriented_score']


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


def parse_beta(raw):
    """Return beta as a finite float, or None if missing/non-numeric/non-finite."""
    try:
        b = float(raw)
    except (TypeError, ValueError):
        return None
    return b if math.isfinite(b) else None


def assign_alleles(hg38_base, a1, a2, beta):
    """Decide which GWAS allele is the hg38 reference, and which is the risk allele.

    Pure function (no I/O), so it can be tested on its own.

    Returns a dict. If the SNP cannot be scored, 'skip_reason' is set and the
    remaining fields may be empty. Otherwise:
        hg38_ref, hg38_alt  bases for the reference and substituted sequences
        hg38_ref_is         'A1' or 'A2'
        palindromic         True for A/T and C/G SNPs
        risk_allele         A1 if beta > 0, A2 if beta < 0, None if no direction
        orientation         +1 if hg38_alt is the risk allele, -1 if hg38_ref is,
                            None if there is no risk direction
    """
    a1, a2, hg38_base = a1.strip().upper(), a2.strip().upper(), hg38_base.upper()
    out = {'hg38_ref': '', 'hg38_alt': '', 'hg38_ref_is': '', 'palindromic': '',
           'risk_allele': None, 'orientation': None, 'skip_reason': ''}

    if a1 not in BASES or a2 not in BASES:
        out['skip_reason'] = f'alleles are not single bases (A1={a1!r}, A2={a2!r})'
        return out
    if a1 == a2:
        out['skip_reason'] = f'A1 and A2 are the same base ({a1})'
        return out
    if hg38_base not in BASES:
        out['skip_reason'] = f'hg38 base at the variant is {hg38_base!r}, not A/C/G/T'
        return out

    if hg38_base == a1:
        out.update(hg38_ref=a1, hg38_alt=a2, hg38_ref_is='A1')
    elif hg38_base == a2:
        out.update(hg38_ref=a2, hg38_alt=a1, hg38_ref_is='A2')
    else:
        out['skip_reason'] = (f'hg38 base {hg38_base} matches neither A1 ({a1}) nor '
                              f'A2 ({a2}): coordinate or strand problem')
        return out

    out['palindromic'] = {a1, a2} in PALINDROMIC_PAIRS

    if beta is not None and beta != 0:
        out['risk_allele'] = a1 if beta > 0 else a2
        out['orientation'] = 1 if out['hg38_alt'] == out['risk_allele'] else -1
    return out


def fmt_score(x):
    """Six decimal places; turn -0.0 into 0.0 so zero never prints as '-0.000000'."""
    return f'{x + 0.0:.6f}'


def main():
    # --- Setup ---
    device = 'mps' if torch.backends.mps.is_available() else (
        'cuda' if torch.cuda.is_available() else 'cpu')
    print(f'Device: {device}')
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    scores_out = os.path.join(OUTPUT_DIR, OUT_NAME)
    variants_out = os.path.join(OUTPUT_DIR, os.path.splitext(OUT_NAME)[0] + '.variants.tsv')

    # --- Load SNPs and check columns BEFORE loading the model ---
    print(f'INPUT : {INPUT_FILE}')
    print(f'OUTPUT: {scores_out}')
    print(f'        {variants_out}')
    with open(INPUT_FILE, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        in_fields = list(reader.fieldnames or [])
        snps = list(reader)

    missing = [c for c in REQUIRED_COLUMNS if c not in in_fields]
    if missing:
        hint = ''
        if 'ref' in in_fields and 'a1' not in in_fields:
            hint = ('\nThis file uses the old ref/alt column names. Regenerate it with '
                    'the current selection scripts and liftover_to_hg38.py.')
        sys.exit(f'{INPUT_FILE} is missing columns {missing}. Columns are: {in_fields}{hint}')
    if not snps:
        sys.exit(f'{INPUT_FILE} has no SNP rows.')
    print(f'Loaded {len(snps)} SNPs.')

    extra_fields = [c for c in in_fields if c not in REQUIRED_COLUMNS]
    variant_fields = (['rsid', 'chr', 'pos', 'a1', 'a2', 'beta'] + extra_fields +
                      ['hg38_ref', 'hg38_alt', 'hg38_ref_is', 'palindromic',
                       'risk_allele', 'status', 'note'])

    # --- Load model ---
    print(f'Loading AlphaGenome from {WEIGHTS_PATH}...')
    model = AlphaGenome.from_pretrained(WEIGHTS_PATH, device=device)
    model.eval()
    print('Model loaded.')

    # --- Process each SNP ---
    all_scores, variant_rows = [], []
    center = SEQ_LEN // 2  # variant position in the extracted window (0-indexed)

    for i, snp in enumerate(snps):
        rsid = snp['rsid']
        chrom = snp['chr'].strip()
        chrom = chrom[3:] if chrom.lower().startswith('chr') else chrom  # Ensembl wants '7'
        pos = int(snp['pos'])        # 1-based hg38 position
        beta = parse_beta(snp['beta'])
        vrow = {k: snp.get(k, '') for k in variant_fields}
        vrow.update(hg38_ref='', hg38_alt='', hg38_ref_is='', palindromic='',
                    risk_allele='', status='', note='')

        print(f'\n--- [{i+1}/{len(snps)}] {rsid} chr{chrom}:{pos} ---')

        # Step 1: Fetch 131,072 bp hg38 window centred on the SNP (1-based, inclusive)
        window_start = pos - center
        window_end = pos + center - 1
        print(f'  Fetching {SEQ_LEN:,} bp window (chr{chrom}:{window_start}-{window_end})...')
        seq = fetch_hg38_sequence(chrom, window_start, window_end)
        if len(seq) != SEQ_LEN:
            note = f'expected {SEQ_LEN} bp from Ensembl, got {len(seq)}'
            print(f'  SKIPPED: {note}')
            variant_rows.append({**vrow, 'status': 'skipped', 'note': note})
            continue

        # Step 2: Which GWAS allele is the hg38 reference? Which is the risk allele?
        a = assign_alleles(seq[center], snp['a1'], snp['a2'], beta)
        if a['skip_reason']:
            print(f'  SKIPPED: {a["skip_reason"]}')
            variant_rows.append({**vrow, 'status': 'skipped', 'note': a['skip_reason']})
            continue
        hg38_ref, hg38_alt = a['hg38_ref'], a['hg38_alt']
        hg38_alt_is = 'A2' if a['hg38_ref_is'] == 'A1' else 'A1'
        risk_txt = ('none (beta is 0 or missing)' if a['risk_allele'] is None else
                    f"{a['risk_allele']} ({'hg38 alt' if a['orientation'] == 1 else 'hg38 ref'})")
        print(f"  GWAS : A1={snp['a1']} (effect)  A2={snp['a2']}  beta={snp['beta']}")
        print(f"  hg38 : ref={hg38_ref} ({a['hg38_ref_is']})  "
              f"alt={hg38_alt} ({hg38_alt_is})  risk allele: {risk_txt}")
        if a['palindromic']:
            print('  NOTE : palindromic SNP: strand cannot be checked against hg38; '
                  'direction assumes the GWAS reported the + strand')

        # Step 3: Build ref and alt sequences
        ref_seq = seq
        alt_seq = seq[:center] + hg38_alt + seq[center + 1:]
        assert ref_seq[center] == hg38_ref
        assert alt_seq[center] == hg38_alt
        assert len(ref_seq) == len(alt_seq) == SEQ_LEN
        assert sum(r != s for r, s in zip(ref_seq, alt_seq)) == 1  # exactly one base differs

        # Step 4: One-hot encode
        ref_tensor = one_hot_encode(ref_seq).to(device)
        alt_tensor = one_hot_encode(alt_seq).to(device)

        # Step 5: Run AlphaGenome
        print('  Running forward pass (ref)...')
        with torch.inference_mode():
            ref_out = model.predict(ref_tensor, organism_index=0)
            print('  Running forward pass (alt)...')
            alt_out = model.predict(alt_tensor, organism_index=0)

        # Step 6: Score every track
        for modality in MODALITIES:
            ref_pred = ref_out[modality][RESOLUTION]   # [1, 131072, n_tracks]
            alt_pred = alt_out[modality][RESOLUTION]
            n_tracks = ref_pred.shape[-1]
            scores = compute_variant_effect(ref_pred, alt_pred, center)

            print(f'  {modality.upper()} ({n_tracks} channels incl. padding): '
                  f'max |score| = {np.max(np.abs(scores)):.4f} '
                  f'(track {int(np.argmax(np.abs(scores)))}), '
                  f'mean |score| = {np.mean(np.abs(scores)):.4f}')

            for t in range(n_tracks):
                s = float(scores[t])
                all_scores.append({
                    'rsid': rsid, 'chr': chrom, 'pos': pos,
                    'a1': snp['a1'], 'a2': snp['a2'], 'beta': snp['beta'],
                    'hg38_ref': hg38_ref, 'hg38_alt': hg38_alt,
                    'hg38_ref_is': a['hg38_ref_is'], 'palindromic': a['palindromic'],
                    'modality': modality, 'track_index': t,
                    'variant_effect_score': fmt_score(s),
                    'risk_oriented_score': ('' if a['orientation'] is None
                                            else fmt_score(a['orientation'] * s)),
                })

        variant_rows.append({**vrow, 'hg38_ref': hg38_ref, 'hg38_alt': hg38_alt,
                             'hg38_ref_is': a['hg38_ref_is'],
                             'palindromic': a['palindromic'],
                             'risk_allele': a['risk_allele'] or '',
                             'status': 'scored'})

        # Free memory between SNPs
        del ref_tensor, alt_tensor, ref_out, alt_out
        if device == 'mps':
            torch.mps.empty_cache()
        elif device == 'cuda':
            torch.cuda.empty_cache()

    # --- Save results ---
    with open(scores_out, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=SCORE_FIELDS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(all_scores)
    with open(variants_out, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=variant_fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(variant_rows)

    scored = [v for v in variant_rows if v['status'] == 'scored']
    skipped = [v for v in variant_rows if v['status'] == 'skipped']
    print('\n=== DONE ===')
    print(f'Scored {len(scored)} of {len(snps)} SNPs; {len(all_scores)} track scores '
          f'written to {scores_out}')
    print(f'hg38 reference is A1 for {sum(v["hg38_ref_is"] == "A1" for v in scored)}, '
          f'A2 for {sum(v["hg38_ref_is"] == "A2" for v in scored)}; '
          f'{sum(v["palindromic"] is True for v in scored)} palindromic')
    print(f'Per-SNP allele table: {variants_out}')
    if skipped:
        print(f'\nWARNING: {len(skipped)} SNP(s) skipped:')
        for v in skipped:
            print(f"  {v['rsid']}: {v['note']}")
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())