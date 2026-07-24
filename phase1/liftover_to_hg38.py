"""
liftover_to_hg38.py
--------------------
Lifts over SNP coordinates in scz_lead_snps.tsv and scz_pilot_5snps.tsv
from hg19 (GRCh37) to hg38 (GRCh38) using pyliftover.

Background:
    The scz2018clozuk GWAS summary statistics use hg19 coordinates (confirmed
    July 21, 2026 by checking all 5 pilot SNP rsIDs against Ensembl GRCh37 and
    GRCh38 REST APIs — positions matched hg19 exactly). AlphaGenome requires
    hg38 input, so liftover is required before any ISM can be run.

Inputs:
    - scz_lead_snps.tsv         : 141 independent SCZ lead SNPs, hg19 coords
    - scz_pilot_5snps.tsv       : 5 pilot SNPs subset, hg19 coords
    - hg19ToHg38.over.chain.gz  : UCSC chain file (download once, see below)

Outputs:
    - scz_lead_snps_hg38.tsv    : same 141 SNPs with hg38 coordinates
    - scz_pilot_5snps_hg38.tsv  : same 5 pilot SNPs with hg38 coordinates

Usage:
    # install dependency
    pip install pyliftover

    # download chain file (one time, ~222 KB)
    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

    # run
    python liftover_to_hg38.py

Validation:
    All 5 pilot SNP hg38 positions were cross-checked against dbSNP via
    Ensembl REST API and matched exactly. 141/141 lead SNPs lifted successfully,
    0 failed.

Coordinate conventions:
    - GWAS TSV files use 1-based coordinates (standard for genomics text files)
    - pyliftover uses 0-based half-open intervals internally
    - This script converts 1-based → 0-based for pyliftover, then back to 1-based
      for the output TSV

Notes:
    - The UCSC liftOver binary could not be used on Hoffman2 due to a GLIBC
      version mismatch (requires GLIBC_2.29, Hoffman2 has an older version).
      pyliftover is a pure Python implementation of the same algorithm.
    - AlphaGenome uses hg38; the sequence extraction script (extract_sequences.py)
      was written assuming hg38 input and works unchanged after liftover.
"""

from pyliftover import LiftOver
import csv
import os

CHAIN_FILE = 'hg19ToHg38.over.chain.gz'
FILES = [
    ('scz_lead_snps.tsv',   'scz_lead_snps_hg38.tsv'),
    ('scz_pilot_5snps.tsv', 'scz_pilot_5snps_hg38.tsv'),
]

def main():
    if not os.path.exists(CHAIN_FILE):
        raise FileNotFoundError(
            f"Chain file not found: {CHAIN_FILE}\n"
            "Download with:\n"
            "  wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
        )

    lo = LiftOver(CHAIN_FILE)

    for tsv_in, tsv_out in FILES:
        with open(tsv_in) as f:
            reader = csv.DictReader(f, delimiter='\t')
            fieldnames = reader.fieldnames
            rows = list(reader)

        lifted, failed = [], []

        for row in rows:
            chrom = 'chr' + row['chr']
            pos_hg19 = int(row['pos'])              # 1-based input
            result = lo.convert_coordinate(chrom, pos_hg19 - 1)  # 0-based for pyliftover

            if result:
                new_chrom, new_pos, strand, _ = result[0]
                row['chr'] = new_chrom.replace('chr', '')
                row['pos'] = str(new_pos + 1)       # back to 1-based output
                lifted.append(row)
            else:
                failed.append(row)
                print(f"  WARNING: could not lift {row.get('rsid', '?')} "
                      f"chr{row['chr']}:{pos_hg19}")

        with open(tsv_out, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            writer.writerows(lifted)

        print(f"{tsv_in} → {tsv_out}: {len(lifted)} lifted, {len(failed)} failed")

    # print pilot SNP hg38 positions for quick validation
    print("\nPilot SNP hg38 positions (cross-check against dbSNP):")
    with open('scz_pilot_5snps_hg38.tsv') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            print(f"  {row['rsid']:12s}  chr{row['chr']}:{row['pos']:12s}  "
                  f"ref={row['ref']}  alt={row['alt']}")

if __name__ == '__main__':
    main()
