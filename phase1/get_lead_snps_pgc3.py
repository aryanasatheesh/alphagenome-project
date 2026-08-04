"""
get_lead_snps_pgc3.py
---------------------
Processes the PGC3 SCZ GWAS summary statistics (Trubetskoy et al., 2022, Nature)
to extract independent lead SNPs for in silico mutagenesis with AlphaGenome.

Background:
    The PGC3 SCZ GWAS is the largest schizophrenia GWAS to date, with ~69,000 cases
    and ~237,000 controls (European ancestry subset). The summary statistics contain
    7,659,767 SNPs, of which 20,457 reach genome-wide significance (p < 5e-8).

    A naive approach of selecting the top N SNPs by p-value fails because the Major
    Histocompatibility Complex (MHC) region on chromosome 6 (~25-34 Mb, hg19)
    dominates the signal: 6,104 of 20,457 genome-wide significant SNPs (30%) are
    on chr6, and the top 1,168 SNPs by p-value are ALL in the MHC. The first
    non-chr6 SNP appears at rank #1,169. This is a known phenomenon — the MHC
    harbors extreme genetic diversity and complex LD structure, producing highly
    significant but difficult-to-interpret associations.

    To obtain a set of independent loci spread across the genome, this script:
    1. Filters to genome-wide significant SNPs (p < 5e-8)
    2. Removes indels (keeps only single-base substitutions, to standardize input for AlphaGenome)
    3. Documents the MHC dominance by reporting how many top SNPs fall in chr6
    4. Excludes the MHC region (chr6:25-34 Mb in hg19) from lead SNP selection
    5. Applies a 1 Mb clumping window per chromosome to define independent loci
    6. Selects pilot SNPs from the resulting independent loci

    The clumping procedure: sort all significant SNPs by p-value (ascending),
    take the most significant as the first lead SNP, then exclude all SNPs within
    1 Mb on the same chromosome, repeat until no SNPs remain. The 1 Mb window is
    standard practice for defining independent GWAS loci (see e.g. PGC methods).

Input:
    PGC3_SCZ_european.vcf.tsv.gz — PGC3 European-ancestry summary statistics
    Downloaded from: https://doi.org/10.6084/m9.figshare.19426775

    File format (PGCsumstatsVCFv1.0):
        CHROM   chromosome (no 'chr' prefix)
        ID      rsID
        POS     base pair position (hg19 / GRCh37)
        A1      effect allele (alt)
        A2      other allele (ref)
        FCAS    frequency of A1 in cases
        FCON    frequency of A1 in controls
        IMPINFO imputation INFO score
        BETA    log odds ratio
        SE      standard error
        PVAL    p-value
        NCAS    number of cases
        NCON    number of controls
        NEFF    effective sample size

    Coordinate system: hg19 (GRCh37) — confirmed by checking rs2007044 position
    (2,344,960 matches hg19; hg38 position is 2,235,794). Liftover to hg38 is
    required before use with AlphaGenome.

Output:
    scz_pgc3_mhc_analysis.tsv  — top 50 SNPs showing MHC dominance (for documentation)
    scz_pgc3_lead_snps.tsv     — all independent lead SNPs, hg19 coordinates
    scz_pgc3_pilot_10snps.tsv  — 10 pilot SNPs for pipeline development
    
    Output columns: rsid, chr, pos, ref, alt, beta, pval

Usage:
    python get_lead_snps_pgc3.py

    Expects PGC3_SCZ_european.vcf.tsv.gz in the data/ subdirectory.
    To change to the primary (multi-ancestry) file, edit INPUT_FILE below.
"""

import csv
import gzip
import os

# === Configuration ===
INPUT_FILE = os.path.join('data', 'PGC3_SCZ_european.vcf.tsv.gz')
OUTPUT_DIR = 'data'
P_THRESHOLD = 5e-8
CLUMP_WINDOW = 1_000_000  # 1 Mb
MHC_CHR = '6'
MHC_START = 25_000_000    # hg19 coordinates
MHC_END = 34_000_000      # hg19 coordinates
N_PILOT = 10
VALID_BASES = {'A', 'C', 'G', 'T'}


def read_gwas(filepath):
    """Read PGC3 VCF-format summary statistics, skipping ## header lines."""
    snps = []
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('CHROM'):
                header = line.strip().split('\t')
                continue
            fields = line.strip().split('\t')
            row = dict(zip(header, fields))
            snps.append(row)
    return snps


def is_snv(a1, a2):
    """Check if variant is a single-nucleotide substitution (not indel)."""
    return (len(a1) == 1 and len(a2) == 1 and
            a1.upper() in VALID_BASES and a2.upper() in VALID_BASES)


def in_mhc(chrom, pos):
    """Check if position falls in the MHC region (chr6:25-34 Mb, hg19)."""
    return chrom == MHC_CHR and MHC_START <= pos <= MHC_END


def clump_snps(snps, window=CLUMP_WINDOW):
    """
    Greedy clumping: select lead SNPs separated by at least `window` bp.
    
    Algorithm:
        1. Sort SNPs by p-value (ascending = most significant first)
        2. Take the top SNP as a lead SNP
        3. Remove all SNPs within `window` bp on the same chromosome
        4. Repeat until no SNPs remain
    
    This is a standard approach for defining independent GWAS loci.
    The 1 Mb window is conventional and accounts for typical LD extent
    in European populations.
    """
    # Sort by p-value
    snps_sorted = sorted(snps, key=lambda x: x['pval'])
    
    leads = []
    used = set()
    
    for i, snp in enumerate(snps_sorted):
        if i in used:
            continue
        leads.append(snp)
        # Mark all SNPs within window on same chromosome as used
        for j, other in enumerate(snps_sorted):
            if j in used or j == i:
                continue
            if (other['chr'] == snp['chr'] and
                abs(other['pos'] - snp['pos']) <= window):
                used.add(j)
    
    return leads


def main():
    print(f"Reading {INPUT_FILE}...")
    raw_snps = read_gwas(INPUT_FILE)
    print(f"Total SNPs in file: {len(raw_snps):,}")
    
    # === Step 1: Parse and filter to genome-wide significance ===
    parsed = []
    for row in raw_snps:
        try:
            pval = float(row['PVAL'])
        except (ValueError, KeyError):
            continue
        if pval >= P_THRESHOLD:
            continue
        
        a1 = row['A1'].upper()
        a2 = row['A2'].upper()
        chrom = row['CHROM']
        pos = int(row['POS'])
        
        parsed.append({
            'rsid': row['ID'],
            'chr': chrom,
            'pos': pos,
            'ref': a2,       # A2 = other allele = reference
            'alt': a1,       # A1 = effect allele = alternate
            'beta': float(row['BETA']),
            'pval': pval,
        })
    
    print(f"Genome-wide significant (p < {P_THRESHOLD}): {len(parsed):,}")
    
    # === Step 2: Remove indels ===
    snvs = [s for s in parsed if is_snv(s['alt'], s['ref'])]
    n_indels = len(parsed) - len(snvs)
    print(f"After removing indels: {len(snvs):,} SNVs ({n_indels:,} indels removed)")
    
    # === Step 3: Document MHC dominance ===
    snvs_sorted = sorted(snvs, key=lambda x: x['pval'])
    
    # How many of top 50 are chr6?
    top50 = snvs_sorted[:50]
    top50_chr6 = sum(1 for s in top50 if s['chr'] == MHC_CHR)
    print(f"\n--- MHC Dominance Analysis ---")
    print(f"Top 50 SNPs by p-value: {top50_chr6}/50 are on chr6")
    
    # Find rank of first non-chr6 SNP
    for i, s in enumerate(snvs_sorted):
        if s['chr'] != MHC_CHR:
            print(f"First non-chr6 SNP: rank #{i+1} ({s['rsid']}, "
                  f"chr{s['chr']}:{s['pos']:,}, p={s['pval']:.3e})")
            break
    
    # Count chr6 vs non-chr6
    n_chr6 = sum(1 for s in snvs if s['chr'] == MHC_CHR)
    n_mhc = sum(1 for s in snvs if in_mhc(s['chr'], s['pos']))
    print(f"GW-sig SNVs on chr6: {n_chr6:,} / {len(snvs):,} "
          f"({100*n_chr6/len(snvs):.1f}%)")
    print(f"GW-sig SNVs in MHC region (chr6:25-34Mb): {n_mhc:,}")
    
    # Save top 50 for documentation
    mhc_out = os.path.join(OUTPUT_DIR, 'scz_pgc3_mhc_analysis.tsv')
    with open(mhc_out, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['rsid','chr','pos','ref','alt','beta','pval'],
                                delimiter='\t')
        writer.writeheader()
        writer.writerows(top50)
    print(f"Saved top 50 SNPs to {mhc_out}")
    
    # === Step 4: Exclude MHC and clump ===
    non_mhc = [s for s in snvs if not in_mhc(s['chr'], s['pos'])]
    print(f"\n--- Lead SNP Selection ---")
    print(f"After excluding MHC: {len(non_mhc):,} SNVs")
    
    leads = clump_snps(non_mhc, window=CLUMP_WINDOW)
    print(f"After 1 Mb clumping: {len(leads)} independent lead SNPs")
    
    # === Step 5: Save outputs ===
    fieldnames = ['rsid', 'chr', 'pos', 'ref', 'alt', 'beta', 'pval']
    
    # All lead SNPs
    leads_out = os.path.join(OUTPUT_DIR, 'scz_pgc3_lead_snps.tsv')
    with open(leads_out, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(leads)
    print(f"Saved {len(leads)} lead SNPs to {leads_out}")
    
    # Pilot SNPs (top N by p-value from leads)
    pilot = leads[:N_PILOT]
    pilot_out = os.path.join(OUTPUT_DIR, 'scz_pgc3_pilot_10snps.tsv')
    with open(pilot_out, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(pilot)
    print(f"Saved {len(pilot)} pilot SNPs to {pilot_out}")
    
    # Print pilot SNPs
    print(f"\nPilot SNPs (top {N_PILOT} lead SNPs by p-value, excluding MHC):")
    print(f"{'rsid':15s} {'chr':>4s} {'pos':>12s} {'ref':>4s} {'alt':>4s} {'pval':>12s}")
    print("-" * 55)
    for s in pilot:
        print(f"{s['rsid']:15s} {s['chr']:>4s} {s['pos']:>12,d} "
              f"{s['ref']:>4s} {s['alt']:>4s} {s['pval']:>12.3e}")
    
    # Reminder about coordinates
    print(f"\nNOTE: All positions are hg19 (GRCh37). Run liftover_to_hg38.py")
    print(f"before using with AlphaGenome.")


if __name__ == '__main__':
    main()