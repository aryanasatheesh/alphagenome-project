"""
get_control_snps.py
-------------------
Selects internal negative-control SNPs from the PGC3 SCZ GWAS for ISM comparison.

Rationale:
    Variant effect scores from AlphaGenome have no inherent scale. A CNS mean
    |score| of 0.16 is only interpretable relative to what a variant with no
    disease association looks like. This script builds that reference set.

    An INTERNAL control (drawn from the same GWAS) is better matched than an
    external one: same ancestry, same imputation panel, same allele-frequency
    spectrum, same file conventions. The only systematic difference is disease
    association. That isolates the variable of interest.

Selection criteria:
    1. p > 0.5 — clearly non-associated. Note this is deliberately far above the
       5e-8 significance threshold: SNPs at p ~ 1e-5 may be real but underpowered
       signal, which would contaminate the null.
    2. Single-nucleotide substitutions only.
    3. Outside the extended MHC (chr6:25-34 Mb, hg19) — same exclusion as the
       lead SNP pipeline, since MHC LD structure is atypical.
    4. At least 1 Mb from any of the 173 lead SNPs. Without this, a "control"
       could sit inside the LD block of a genuine association.
    5. Spread across chromosomes (at most one per chromosome by default) to
       avoid clustering in a single genomic neighbourhood.

    Selection is otherwise RANDOM. A future version should match controls to the
    pilot SNPs on allele frequency and GC content — Li & Ernst (Genome Biology,
    2025) showed GC content alone explains much of the apparent variant signal in
    deep learning models, so unmatched controls may differ in sequence
    composition rather than in disease relevance.

Caveat on sample size:
    The default n=10 matches the pilot set but yields a weak null. Estimating a
    mean and spread from 10 points gives wide uncertainty. Treat n=10 as a
    machinery check; use 50+ before drawing conclusions.

Input:
    data/PGC3_SCZ_european.vcf.tsv.gz  — PGC3 summary statistics (hg19)
    data/scz_pgc3_lead_snps.tsv        — 173 lead SNPs, used for exclusion

Output:
    data/scz_pgc3_control_10snps.tsv       — control SNPs, hg19
    data/scz_pgc3_control_10snps_hg38.tsv  — control SNPs, hg38 (if pyliftover available)

    Same column schema as the pilot files so the ISM script can consume either.

Usage:
    python get_control_snps.py [n_controls] [random_seed]
"""

import csv
import gzip
import os
import random
import sys

# === Configuration ===
GWAS_FILE = os.path.join('data', 'PGC3_SCZ_european.vcf.tsv.gz')
LEAD_FILE = os.path.join('data', 'scz_pgc3_lead_snps.tsv')
OUTPUT_DIR = 'data'
CHAIN_FILE = os.path.join('data', 'hg19ToHg38.over.chain.gz')

P_NULL_MIN = 0.5           # SNPs must be clearly non-associated
MHC_CHR = '6'
MHC_START = 25_000_000
MHC_END = 34_000_000
LEAD_EXCLUSION = 1_000_000  # keep controls >=1 Mb from any lead SNP
MAX_PER_CHROM = 1
VALID_BASES = {'A', 'C', 'G', 'T'}

N_CONTROLS = int(sys.argv[1]) if len(sys.argv) > 1 else 10
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 42


def load_lead_positions(path):
    """Return {chrom: [positions]} for the lead SNPs, used for exclusion."""
    leads = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            leads.setdefault(row['chr'], []).append(int(row['pos']))
    return leads


def near_lead(chrom, pos, leads, window=LEAD_EXCLUSION):
    """True if this position is within `window` bp of any lead SNP."""
    return any(abs(pos - lp) <= window for lp in leads.get(chrom, ()))


def is_snv(a1, a2):
    return len(a1) == 1 and len(a2) == 1 and a1 in VALID_BASES and a2 in VALID_BASES


def in_mhc(chrom, pos):
    return chrom == MHC_CHR and MHC_START <= pos <= MHC_END


def main():
    rng = random.Random(SEED)
    leads = load_lead_positions(LEAD_FILE)
    n_leads = sum(len(v) for v in leads.values())
    print(f'Loaded {n_leads} lead SNPs for exclusion (>={LEAD_EXCLUSION:,} bp)')

    # Reservoir sample per chromosome while streaming, so we never hold
    # all 7.6M rows in memory (login nodes have limited RAM).
    reservoir = {}   # chrom -> candidate row
    seen = {}        # chrom -> count of eligible SNPs seen
    total = n_eligible = 0
    header = None

    print(f'Streaming {GWAS_FILE}...')
    with gzip.open(GWAS_FILE, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('CHROM'):
                header = line.strip().split('\t')
                continue
            total += 1
            fields = line.rstrip('\n').split('\t')
            row = dict(zip(header, fields))

            try:
                pval = float(row['PVAL'])
            except (ValueError, KeyError):
                continue
            if pval <= P_NULL_MIN:
                continue

            a1, a2 = row['A1'].upper(), row['A2'].upper()
            if not is_snv(a1, a2):
                continue

            chrom = row['CHROM']
            try:
                pos = int(row['POS'])
            except ValueError:
                continue

            if in_mhc(chrom, pos):
                continue
            if near_lead(chrom, pos, leads):
                continue

            n_eligible += 1
            # Reservoir sampling with k=1 per chromosome:
            # each eligible SNP on a chromosome has equal probability of
            # being the one retained, without storing them all.
            seen[chrom] = seen.get(chrom, 0) + 1
            if rng.random() < 1.0 / seen[chrom]:
                reservoir[chrom] = {
                    'rsid': row['ID'],
                    'chr': chrom,
                    'pos': pos,
                    'ref': a2,
                    'alt': a1,
                    'beta': float(row['BETA']),
                    'pval': pval,
                }

    print(f'Total SNPs scanned:   {total:,}')
    print(f'Eligible as controls: {n_eligible:,}')
    print(f'Chromosomes with a candidate: {len(reservoir)}')

    # Pick N chromosomes at random, one control each
    chroms = sorted(reservoir, key=lambda c: int(c) if c.isdigit() else 99)
    if len(chroms) < N_CONTROLS:
        print(f'WARNING: only {len(chroms)} chromosomes available, '
              f'requested {N_CONTROLS}')
    chosen = rng.sample(chroms, min(N_CONTROLS, len(chroms)))
    chosen.sort(key=lambda c: int(c) if c.isdigit() else 99)
    controls = [reservoir[c] for c in chosen]

    fieldnames = ['rsid', 'chr', 'pos', 'ref', 'alt', 'beta', 'pval']
    out_hg19 = os.path.join(OUTPUT_DIR, f'scz_pgc3_control_{N_CONTROLS}snps.tsv')
    with open(out_hg19, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        w.writeheader()
        w.writerows(controls)
    print(f'\nSaved {len(controls)} control SNPs (hg19) to {out_hg19}')

    print(f'\n{"rsid":15s} {"chr":>4s} {"pos":>13s} {"ref":>4s} {"alt":>4s} {"pval":>8s}')
    print('-' * 52)
    for s in controls:
        print(f'{s["rsid"]:15s} {s["chr"]:>4s} {s["pos"]:>13,d} '
              f'{s["ref"]:>4s} {s["alt"]:>4s} {s["pval"]:>8.3f}')

    # --- Liftover to hg38 ---
    try:
        from pyliftover import LiftOver
    except ImportError:
        print('\npyliftover not installed; skipping liftover.')
        print('  pip install pyliftover, then run liftover_to_hg38.py')
        return

    if not os.path.exists(CHAIN_FILE):
        print(f'\nChain file not found at {CHAIN_FILE}; skipping liftover.')
        return

    lo = LiftOver(CHAIN_FILE)
    lifted, failed = [], []
    for row in controls:
        r = dict(row)
        res = lo.convert_coordinate('chr' + r['chr'], r['pos'] - 1)  # 1-based -> 0-based
        if res:
            new_chrom, new_pos, _, _ = res[0]
            r['chr'] = new_chrom.replace('chr', '')
            r['pos'] = new_pos + 1                                    # back to 1-based
            lifted.append(r)
        else:
            failed.append(r)

    out_hg38 = os.path.join(OUTPUT_DIR, f'scz_pgc3_control_{N_CONTROLS}snps_hg38.tsv')
    with open(out_hg38, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        w.writeheader()
        w.writerows(lifted)
    print(f'\nLiftover: {len(lifted)} succeeded, {len(failed)} failed')
    print(f'Saved to {out_hg38}')
    if failed:
        for r in failed:
            print(f'  FAILED: {r["rsid"]} chr{r["chr"]}:{r["pos"]}')


if __name__ == '__main__':
    main()