#!/usr/bin/env python3
"""Draw a random sample of non-associated GWAS variants as matched-pool candidates.

Step 6 of the matched-SNP design (Choo, Nov 18): the pool for each bin holds
10-20 SNPs with no association to the trait (p > 0.5). This script makes the
candidate list. It streams the whole summary-statistics file once, keeps the
variants that pass every rule below, and keeps a random sample of them
(reservoir sampling: every eligible variant has the same chance, fixed seed).
The sample then goes through the same steps as the leads (liftover, distance,
MAF, bins), so leads and pool are described identically.

Why a random sample and not all eligible variants: a few million variants
would need a genome-wide MAF extraction from 1000 Genomes (hours), while
20,000 give hundreds of candidates per bin. The sample has a second use: as
random common variants, its MAF and distance quantiles are the "genome-wide"
bin edges that would let one pool serve every trait.

Rules (each counted in the log):
  - p-value > --min-p (default 0.5): not associated
  - A1 and A2 single bases A/C/G/T (SNVs)
  - autosomes 1-22
  - MAF >= --min-maf (default 0.01), from the GWAS control frequency
  - imputation INFO >= --min-info (default 0.9), if the file has the column
  - outside the MHC (--mhc, file's build; default hg19 chr6:25-34 Mb)
  - more than --lead-window bp (default 1 Mb) from every lead in --leads
    (lead positions in the file's build; for the SCZ hg38 lead table these are
    the chr_hg19 / pos_hg19 columns), so candidates are not in LD with a lead

Output: the lead-table format (rsid chr pos a1 a2 a1_freq_controls beta pval,
plus impinfo), in the file's build, ready for liftover_to_hg38.py.

Usage (PGC3 defaults; from the repo root)
  python3 phase1/sample_pool_candidates.py --gwas phase1/data/PGC3_SCZ_european.vcf.tsv.gz \\
      --leads phase1/data/scz_pgc3_lead_snps_hg38.tsv \\
      --out phase1/data/pool_candidates_scz_20k_hg19.tsv --n 20000 --seed 2026
  Column options (defaults are PGC3's): --chr-col CHROM --pos-col POS --rsid-col ID
      --a1-col A1 --a2-col A2 --freq-col FCON --info-col IMPINFO --beta-col BETA --p-col PVAL
  Lead columns: --lead-chr-col chr_hg19 --lead-pos-col pos_hg19

Requires: Python 3 standard library.
"""

import argparse
import bisect
import csv
import gzip
import hashlib
import random
import sys
from datetime import datetime
from pathlib import Path

BASES = set("ACGT")


def norm_chr(c):
    c = c.strip()
    return c[3:] if c.lower().startswith("chr") else c


def md5(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if str(path).endswith(".gz") \
        else open(path, encoding="utf-8", newline="")


def main():
    ap = argparse.ArgumentParser(description="Random sample of non-associated GWAS variants (pool candidates).")
    ap.add_argument("--gwas", required=True)
    ap.add_argument("--leads", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--min-p", type=float, default=0.5)
    ap.add_argument("--min-maf", type=float, default=0.01)
    ap.add_argument("--min-info", type=float, default=0.9)
    ap.add_argument("--mhc", default="chr6:25000000-34000000")
    ap.add_argument("--lead-window", type=int, default=1_000_000)
    for k, v in [("chr", "CHROM"), ("pos", "POS"), ("rsid", "ID"), ("a1", "A1"), ("a2", "A2"),
                 ("freq", "FCON"), ("info", "IMPINFO"), ("beta", "BETA"), ("p", "PVAL")]:
        ap.add_argument(f"--{k}-col", default=v)
    ap.add_argument("--lead-chr-col", default="chr_hg19")
    ap.add_argument("--lead-pos-col", default="pos_hg19")
    args = ap.parse_args()

    gwas, leads, out = Path(args.gwas), Path(args.leads), Path(args.out)
    for p in (gwas, leads):
        if not p.exists():
            sys.exit(f"Not found: {p}")
    log_path = out.with_name((out.name[:-4] if out.name.endswith(".tsv") else out.name) + ".log.txt")
    lines = []

    def log(msg=""):
        print(msg)
        lines.append(msg)

    log(f"sample_pool_candidates.py  {datetime.now().isoformat(timespec='seconds')}")
    log(f"GWAS: {gwas}  ({gwas.stat().st_size:,} bytes, md5 {md5(gwas)})")

    # leads, in the GWAS file's build
    with open(leads, newline="", encoding="utf-8") as f:
        rd = csv.DictReader(f, delimiter="\t")
        for c in (args.lead_chr_col, args.lead_pos_col):
            if c not in rd.fieldnames:
                sys.exit(f"{leads} has no '{c}' column; set --lead-chr-col/--lead-pos-col "
                         f"to lead positions in the GWAS file's build")
        lead_pos = {}
        n_leads = 0
        for r in rd:
            lead_pos.setdefault(norm_chr(r[args.lead_chr_col]), []).append(int(r[args.lead_pos_col]))
            n_leads += 1
    for v in lead_pos.values():
        v.sort()
    log(f"Leads: {leads} ({n_leads}; columns {args.lead_chr_col}/{args.lead_pos_col}); "
        f"window +/- {args.lead_window:,} bp")
    mhc_c, rng = args.mhc.split(":")
    mhc_c = norm_chr(mhc_c)
    mhc_lo, mhc_hi = (int(x) for x in rng.split("-"))
    log(f"MHC excluded: chr{mhc_c}:{mhc_lo:,}-{mhc_hi:,}")

    def near_lead(c, p):
        v = lead_pos.get(c)
        if not v:
            return False
        j = bisect.bisect_left(v, p)
        return (j < len(v) and v[j] - p <= args.lead_window) or (j > 0 and p - v[j - 1] <= args.lead_window)

    rnd = random.Random(args.seed)
    keep = []
    n_elig = 0
    cnt = {"rows": 0, "p <= min-p or missing": 0, "not SNV": 0, "not autosomal": 0,
           "MAF < min or missing": 0, "INFO < min or missing": 0, "in MHC": 0, "near a lead": 0}
    with open_text(gwas) as f:
        header = None
        for line in f:
            if line.startswith("##"):
                continue
            header = line.rstrip("\n").lstrip("#").split("\t")
            break
        if header is None:
            sys.exit("No header line found")
        H = {c: i for i, c in enumerate(header)}
        need = [args.chr_col, args.pos_col, args.rsid_col, args.a1_col, args.a2_col,
                args.freq_col, args.beta_col, args.p_col]
        miss = [c for c in need if c not in H]
        if miss:
            sys.exit(f"Missing columns {miss}; header is: {header}")
        has_info = args.info_col in H
        log(f"Columns: {', '.join(header)}")
        log(f"INFO filter: {'>= ' + str(args.min_info) + ' on ' + args.info_col if has_info else 'column absent, not applied'}")
        ci, pi, ri, a1i, a2i, fi, bi, pvi = (H[c] for c in need)
        ii = H.get(args.info_col)
        autos = {str(i) for i in range(1, 23)}
        for line in f:
            r = line.rstrip("\n").split("\t")
            cnt["rows"] += 1
            try:
                pv = float(r[pvi])
            except (ValueError, IndexError):
                cnt["p <= min-p or missing"] += 1
                continue
            if not pv > args.min_p:
                cnt["p <= min-p or missing"] += 1
                continue
            a1, a2 = r[a1i].upper(), r[a2i].upper()
            if len(a1) != 1 or len(a2) != 1 or a1 not in BASES or a2 not in BASES or a1 == a2:
                cnt["not SNV"] += 1
                continue
            c = norm_chr(r[ci])
            if c not in autos:
                cnt["not autosomal"] += 1
                continue
            try:
                fr = float(r[fi])
            except ValueError:
                cnt["MAF < min or missing"] += 1
                continue
            if min(fr, 1 - fr) < args.min_maf:
                cnt["MAF < min or missing"] += 1
                continue
            if has_info:
                try:
                    ok = float(r[ii]) >= args.min_info
                except ValueError:
                    ok = False
                if not ok:
                    cnt["INFO < min or missing"] += 1
                    continue
            p = int(r[pi])
            if c == mhc_c and mhc_lo <= p <= mhc_hi:
                cnt["in MHC"] += 1
                continue
            if near_lead(c, p):
                cnt["near a lead"] += 1
                continue
            n_elig += 1
            rec = [r[ri], c, str(p), a1, a2, r[fi], r[bi], r[pvi], r[ii] if has_info else ""]
            if len(keep) < args.n:
                keep.append(rec)
            else:
                j = rnd.randrange(n_elig)
                if j < args.n:
                    keep[j] = rec
    log(f"Rows read: {cnt.pop('rows'):,}")
    for k, v in cnt.items():
        log(f"  removed, {k}: {v:,}")
    log(f"Eligible: {n_elig:,}; sampled {len(keep):,} (seed {args.seed})")

    keep.sort(key=lambda x: (int(x[1]), int(x[2])))
    with open(out, "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f, delimiter="\t", lineterminator="\n")
        wr.writerow(["rsid", "chr", "pos", "a1", "a2", "a1_freq_controls", "beta", "pval", "impinfo"])
        wr.writerows(keep)
    log(f"Written: {out}")
    log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Log: {log_path}")


if __name__ == "__main__":
    main()
