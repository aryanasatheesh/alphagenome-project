#!/usr/bin/env python3
"""
survey_liftover_edge_cases.py
-----------------------------
Count how often each liftover edge case occurs across an entire GWAS summary
statistics file (hg19 -> hg38), using the SAME lift_one() function as
liftover_to_hg38.py, so the survey tests exactly the code the pipeline runs.

Why a separate script
    liftover_to_hg38.py keeps every row in memory, which is right for the few
    hundred SNPs the pipeline lifts but would run out of memory on millions of
    rows. This script streams the file one line at a time and keeps only counts,
    so memory stays flat no matter how big the GWAS is.

What it reports, for three subsets of SNPs
    all            every SNP in the file
    gw_sig         p < 5e-8 (the pool lead SNPs are drawn from)
    null_p_gt_0.5  p > 0.5 (the pool internal control SNPs are drawn from)
  - number of SNPs in each liftover status (lifted, unmapped, multi, ...)
  - number of palindromic (A/T, C/G) SNPs, which matters for allele orientation

Outputs (prefix set by --out-prefix)
    <prefix>.summary.tsv          status x subset counts and percentages (small; commit it)
    <prefix>.gw_sig_failures.tsv  every genome-wide significant SNP that did not lift (small; commit it)
    <prefix>.all_failures.tsv.gz  every SNP that did not lift (large; gitignored)
    <prefix>.log.txt              provenance: command, chain md5, versions, runtime

Usage (Hoffman2, from the repo root; takes a few minutes)
    python phase1/survey_liftover_edge_cases.py \\
        --gwas data/PGC3_SCZ_european.vcf.tsv.gz \\
        --out-prefix results/liftover_survey/pgc3_european

Column names default to the PGC3 VCF-style format (CHROM, ID, POS, A1, A2, PVAL)
and can be changed with --chr-col, --id-col, --pos-col, --a1-col, --a2-col,
--p-col. Lines starting with '##' are skipped; the header line may start with
'#' or not.
"""

import argparse
import datetime
import gzip
import shlex
import sys
import time
from collections import Counter
from pathlib import Path

# Import from the strict liftover script that sits in the same folder.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from liftover_to_hg38 import (DEFAULT_CHAIN, STATUS_REASONS, LiftOver, lift_one,  # noqa: E402
                              md5sum, pyliftover)

SUBSETS = ["all", "gw_sig", "null_p_gt_0.5"]
GW_SIG = 5e-8
NULL_P = 0.5
PALINDROMIC = ({"A", "T"}, {"C", "G"})
FAIL_FIELDS = ["id", "chr_hg19", "pos_hg19", "a1", "a2", "p", "status",
               "n_mappings", "chr_hg38", "pos_hg38", "strand", "note"]


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def main():
    ap = argparse.ArgumentParser(description="Count liftover edge cases across a whole GWAS file.")
    ap.add_argument("--gwas", required=True, help="GWAS summary statistics (.tsv or .tsv.gz)")
    ap.add_argument("--out-prefix", required=True, help="output path prefix")
    ap.add_argument("--chain", default=str(DEFAULT_CHAIN))
    ap.add_argument("--chr-col", default="CHROM")
    ap.add_argument("--id-col", default="ID")
    ap.add_argument("--pos-col", default="POS")
    ap.add_argument("--a1-col", default="A1")
    ap.add_argument("--a2-col", default="A2")
    ap.add_argument("--p-col", default="PVAL")
    ap.add_argument("--progress-every", type=int, default=1_000_000)
    args = ap.parse_args()

    gwas, chain = Path(args.gwas), Path(args.chain)
    if not gwas.exists():
        sys.exit(f"GWAS file not found: {gwas}")
    if not chain.exists():
        sys.exit(f"Chain file not found: {chain}")
    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    print(f"GWAS  : {gwas}")
    print(f"CHAIN : {chain}")
    lo = LiftOver(str(chain))

    status_counts = Counter()   # (subset, status) -> n
    subset_totals = Counter()   # subset -> n
    palin_counts = Counter()    # subset -> n palindromic
    bad_p = 0
    header, cols = None, None

    all_fail_path = Path(f"{prefix}.all_failures.tsv.gz")
    gw_fail_path = Path(f"{prefix}.gw_sig_failures.tsv")
    with open_text(gwas) as fin, \
            gzip.open(all_fail_path, "wt") as f_all, \
            open(gw_fail_path, "w") as f_gw:
        f_all.write("\t".join(FAIL_FIELDS) + "\n")
        f_gw.write("\t".join(FAIL_FIELDS) + "\n")

        n = 0
        for line in fin:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = [h.lstrip("#") for h in fields]
                wanted = [args.chr_col, args.id_col, args.pos_col,
                          args.a1_col, args.a2_col, args.p_col]
                missing = [c for c in wanted if c not in header]
                if missing:
                    sys.exit(f"Columns {missing} not in header: {header}")
                cols = {c: header.index(c) for c in wanted}
                continue
            if len(fields) < len(header):
                continue  # truncated line; counted nowhere, reported below via totals

            n += 1
            chrom = fields[cols[args.chr_col]]
            pos = fields[cols[args.pos_col]]
            a1 = fields[cols[args.a1_col]].upper()
            a2 = fields[cols[args.a2_col]].upper()
            try:
                p = float(fields[cols[args.p_col]])
            except ValueError:
                p = None
                bad_p += 1

            subsets = ["all"]
            if p is not None and p < GW_SIG:
                subsets.append("gw_sig")
            if p is not None and p > NULL_P:
                subsets.append("null_p_gt_0.5")

            res = lift_one(lo, chrom, pos)
            palin = {a1, a2} in PALINDROMIC
            for s in subsets:
                subset_totals[s] += 1
                status_counts[(s, res["status"])] += 1
                if palin:
                    palin_counts[s] += 1

            if res["status"] != "lifted":
                row = [fields[cols[args.id_col]], chrom, pos, a1, a2,
                       fields[cols[args.p_col]], res["status"], str(res["n_mappings"]),
                       res["chr_hg38"], "" if res["pos_hg38"] is None else str(res["pos_hg38"]),
                       res["strand"], res["note"]]
                f_all.write("\t".join(row) + "\n")
                if "gw_sig" in subsets:
                    f_gw.write("\t".join(row) + "\n")

            if args.progress_every and n % args.progress_every == 0:
                print(f"  {n:>11,} SNPs  ({time.time() - t0:.0f} s)", flush=True)

    if header is None or n == 0:
        sys.exit(f"No data rows found in {gwas}")

    # ---- summary table ----
    summary_path = Path(f"{prefix}.summary.tsv")
    statuses = [s for s in STATUS_REASONS if any(status_counts[(sub, s)] for sub in SUBSETS)]
    with open(summary_path, "w") as f:
        f.write("status\t" + "\t".join(f"{s}_n\t{s}_pct" for s in SUBSETS) + "\tmeaning\n")
        for st in statuses + ["TOTAL", "palindromic"]:
            cells = []
            for sub in SUBSETS:
                tot = subset_totals[sub]
                k = (tot if st == "TOTAL" else palin_counts[sub] if st == "palindromic"
                     else status_counts[(sub, st)])
                cells += [str(k), f"{100 * k / tot:.4f}" if tot else ""]
            meaning = STATUS_REASONS.get(st, "A/T or C/G SNP; strand not checkable"
                                         if st == "palindromic" else "")
            f.write(f"{st}\t" + "\t".join(cells) + f"\t{meaning}\n")

    elapsed = time.time() - t0
    with open(f"{prefix}.log.txt", "w") as f:
        f.write(f"date              : {datetime.datetime.now().isoformat(timespec='seconds')}\n")
        f.write(f"command           : {' '.join(shlex.quote(a) for a in sys.argv)}\n")
        f.write(f"gwas              : {gwas}\n")
        f.write(f"chain file        : {chain.name}\n")
        f.write(f"chain md5         : {md5sum(chain)}\n")
        f.write(f"pyliftover version: {getattr(pyliftover, '__version__', 'unknown')}\n")
        f.write(f"python version    : {sys.version.split()[0]}\n")
        f.write(f"snps read         : {n}\n")
        f.write(f"non-numeric p     : {bad_p}\n")
        f.write(f"runtime_seconds   : {elapsed:.0f}\n")

    # ---- console report ----
    print(f"\nRead {n:,} SNPs in {elapsed:.0f} s  (non-numeric p-values: {bad_p})\n")
    print(f"{'status':<15}" + "".join(f"{s:>24}" for s in SUBSETS))
    for st in statuses + ["TOTAL", "palindromic"]:
        line = f"{st:<15}"
        for sub in SUBSETS:
            tot = subset_totals[sub]
            k = (tot if st == "TOTAL" else palin_counts[sub] if st == "palindromic"
                 else status_counts[(sub, st)])
            line += f"{k:>13,} ({100 * k / tot:6.3f}%)" if tot else f"{'':>24}"
        print(line)
    print(f"\nSUMMARY : {summary_path}")
    print(f"GW-SIG  : {gw_fail_path}")
    print(f"ALL     : {all_fail_path}  (gitignored)")
    return 0


if __name__ == "__main__":
    sys.exit(main())