#!/usr/bin/env python3
"""Select independent lead SNPs from any GWAS summary-statistics file.

This is the general version of the lead-selection step: the same rules that
produced the 173 PGC3 schizophrenia leads, applied to any trait through a
column map. It was written for the type 1 diabetes control GWAS (Chiou et al.
2021, GCST90014023) but takes any tab- or space-separated file, gzipped or not.

Rules, in order (each step's count goes in the log):
  1. p < --p-thresh (default 5e-8)
  2. single-base A/C/G/T alleles on both sides (ISM scores SNVs only)
  3. effect-allele frequency within [--min-freq, 1 - --min-freq], if a
     frequency column is given (default 0.01: common variants, like the SCZ leads)
  4. autosomes only (chr1-22), like the PGC3 autosome file (--keep-sex-chr to turn off)
  5. outside the MHC (--mhc, in the file's own build; default hg19 chr6:25-34 Mb)
  6. greedy distance clumping: walk SNPs from smallest p; keep one unless a
     kept lead on the same chromosome lies within --window bp (default 1 Mb)

Output columns (what the strict liftover step and the ISM scripts read):
  rsid  chr  pos  a1  a2  [a1_freq]  beta  pval  trait
  (same names and order as phase1/data/scz_pgc3_lead_snps.tsv; a1_freq is the
  effect-allele frequency as the file reports it)
  a1 = effect allele, a2 = other allele, beta = effect of a1 (log OR if the
  file gives an odds ratio). So beta > 0 means a1 is the risk allele, the same
  convention as PGC3.

Usage
  # 1. Look at the file first: prints the header, two rows and a guessed column map
  python3 phase1/select_leads.py --gwas t1d.tsv.gz --inspect

  # 2. Select leads (column names from step 1)
  python3 phase1/select_leads.py --gwas t1d.tsv.gz --trait T1D \\
      --chr-col chromosome --pos-col base_pair_location \\
      --ea-col effect_allele --oa-col other_allele \\
      --beta-col beta --p-col p_value --rsid-col variant_id \\
      --freq-col effect_allele_frequency \\
      --out phase1/data/t1d_chiou2021_lead_snps.tsv \\
      --sample 20 --seed 2026 --sample-out phase1/data/t1d_chiou2021_control20.tsv

  # 3. Check the rules reproduce an existing lead list (e.g. the 173 SCZ leads)
  python3 phase1/select_leads.py --gwas <PGC3 file> ... --compare phase1/data/scz_pgc3_lead_snps.tsv
"""
import argparse
import csv
import gzip
import hashlib
import math
import random
import re
import sys
from datetime import datetime
from pathlib import Path

BASES = set("ACGT")

# Common header spellings, used only by --inspect to suggest a column map.
GUESSES = {
    "chr": ["chromosome", "chrom", "chr", "#chrom", "hm_chrom", "chr_name"],
    "pos": ["base_pair_location", "pos", "bp", "position", "hm_pos", "chr_position"],
    "ea": ["effect_allele", "a1", "alt", "hm_effect_allele", "allele1", "ea"],
    "oa": ["other_allele", "a2", "ref", "hm_other_allele", "allele2", "nea"],
    "beta": ["beta", "hm_beta", "effect", "b", "logor"],
    "or": ["odds_ratio", "or", "hm_odds_ratio"],
    "p": ["p_value", "p", "pval", "pvalue", "p-value", "p.value"],
    "neglog10p": ["neg_log_10_p_value", "mlog10p", "log10p"],
    "rsid": ["variant_id", "rsid", "snp", "id", "hm_rsid", "markername", "rs_id"],
    "freq": ["effect_allele_frequency", "eaf", "freq", "af", "hm_effect_allele_frequency", "fcon"],
}


def open_text(path):
    """Open plain or gzipped text (checked by magic bytes, not the file name)."""
    with open(path, "rb") as f:
        gz = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt", newline="") if gz else open(path, newline="")


def read_header(path):
    """Return (header fields, separator, number of '##' lines skipped)."""
    skipped = 0
    with open_text(path) as f:
        for line in f:
            if line.startswith("##"):
                skipped += 1
                continue
            sep = "\t" if "\t" in line else None  # None = any whitespace
            return line.rstrip("\r\n").split(sep), sep, skipped
    sys.exit(f"{path}: no header line found")


def iter_rows(path, sep, skipped):
    with open_text(path) as f:
        for _ in range(skipped):
            next(f)
        next(f)  # header
        for line in f:
            if line.strip():
                yield line.rstrip("\r\n").split(sep)


def md5(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def norm_chr(c):
    c = c.strip()
    return c[3:] if c.lower().startswith("chr") else c


def parse_region(s):
    m = re.fullmatch(r"(?:chr)?(\w+):(\d+)-(\d+)", s.replace(",", ""))
    if not m:
        sys.exit(f"--mhc must look like chr6:25000000-34000000, got {s!r}")
    return m.group(1), int(m.group(2)), int(m.group(3))


def inspect(path):
    header, sep, skipped = read_header(path)
    print(f"FILE   : {path}")
    print(f"FORMAT : {'tab' if sep == chr(9) else 'whitespace'}-separated; {skipped} '##' lines skipped")
    print(f"COLUMNS ({len(header)}): {header}")
    rows = iter_rows(path, sep, skipped)
    for i in range(2):
        try:
            r = next(rows)
        except StopIteration:
            break
        print(f"ROW {i + 1} : " + " | ".join(f"{h}={v}" for h, v in zip(header, r)))
    lower = {h.lower(): h for h in header}
    print("\nSuggested column map (check it against the rows above):")
    flags = {"chr": "--chr-col", "pos": "--pos-col", "ea": "--ea-col", "oa": "--oa-col",
             "beta": "--beta-col", "or": "--or-col", "p": "--p-col",
             "neglog10p": "--p-col (with --p-is-neglog10)", "rsid": "--rsid-col",
             "freq": "--freq-col"}
    for key, names in GUESSES.items():
        hit = next((lower[n] for n in names if n in lower), None)
        print(f"  {flags[key]:<30} {hit if hit else '(not found)'}")
    print("\nNotes:\n  - Check which allele the effect refers to (the effect allele, usually"
          " 'effect_allele' or 'A1').\n  - Check the genome build in the study's README/yaml:"
          " --mhc must use the same build.")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gwas", required=True, help="summary statistics (.tsv / .txt, gzipped or not)")
    ap.add_argument("--inspect", action="store_true", help="print header + 2 rows + suggested columns, then exit")
    ap.add_argument("--trait", help="short trait label written to every output row, e.g. T1D")
    ap.add_argument("--chr-col"); ap.add_argument("--pos-col")
    ap.add_argument("--ea-col", help="effect allele"); ap.add_argument("--oa-col", help="other allele")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--beta-col", help="effect size on the log-odds (or linear) scale")
    g.add_argument("--or-col", help="odds ratio; converted to beta = ln(OR)")
    ap.add_argument("--p-col")
    ap.add_argument("--p-is-neglog10", action="store_true", help="--p-col holds -log10(p)")
    ap.add_argument("--rsid-col", help="variant ID; non-rs IDs are replaced by chr:pos:a1:a2")
    ap.add_argument("--freq-col", help="effect-allele frequency (optional)")
    ap.add_argument("--p-thresh", type=float, default=5e-8)
    ap.add_argument("--min-freq", type=float, default=0.01,
                    help="drop SNPs with effect-allele frequency outside [f, 1-f]; needs --freq-col")
    ap.add_argument("--mhc", default="chr6:25000000-34000000",
                    help="region to exclude, same build as the file (default hg19 extended MHC)")
    ap.add_argument("--keep-sex-chr", action="store_true", help="keep chrX/Y (default: autosomes only)")
    ap.add_argument("--window", type=int, default=1_000_000, help="clumping distance in bp")
    ap.add_argument("--out", help="TSV of all leads")
    ap.add_argument("--sample", type=int, help="also write a random sample of N leads")
    ap.add_argument("--seed", type=int, help="random seed for --sample (required with it)")
    ap.add_argument("--sample-out", help="TSV for the sample")
    ap.add_argument("--compare", help="existing lead TSV with an 'rsid' column: report overlap")
    args = ap.parse_args()

    gwas = Path(args.gwas)
    if not gwas.exists():
        sys.exit(f"File not found: {gwas}")
    if args.inspect:
        inspect(gwas)
        return

    need = {"--trait": args.trait, "--chr-col": args.chr_col, "--pos-col": args.pos_col,
            "--ea-col": args.ea_col, "--oa-col": args.oa_col, "--p-col": args.p_col, "--out": args.out}
    missing = [k for k, v in need.items() if not v]
    if not (args.beta_col or args.or_col):
        missing.append("--beta-col or --or-col")
    if missing:
        sys.exit(f"Missing: {', '.join(missing)}. Run with --inspect first to see the columns.")
    if args.sample and (args.seed is None or not args.sample_out):
        sys.exit("--sample needs --seed and --sample-out, so the sample can be reproduced.")
    out = Path(args.out)
    if out.resolve() == gwas.resolve():
        sys.exit("Refusing to overwrite the input file.")

    header, sep, skipped = read_header(gwas)
    idx = {h: i for i, h in enumerate(header)}
    cols = {"chr": args.chr_col, "pos": args.pos_col, "ea": args.ea_col, "oa": args.oa_col,
            "eff": args.beta_col or args.or_col, "p": args.p_col,
            "rsid": args.rsid_col, "freq": args.freq_col}
    for k, name in cols.items():
        if name and name not in idx:
            sys.exit(f"Column {name!r} not in the file. Columns are: {header}")
    ci = {k: idx[v] for k, v in cols.items() if v}
    mhc_chr, mhc_lo, mhc_hi = parse_region(args.mhc)
    if args.freq_col is None:
        print("NOTE: no --freq-col, so no frequency filter is applied.")

    counts = {"rows read": 0, "malformed (skipped)": 0, "p < threshold": 0,
              "not single-base A/C/G/T": 0, "frequency outside range": 0, "not autosomal": 0, "in MHC": 0}
    cands = []
    for r in iter_rows(gwas, sep, skipped):
        counts["rows read"] += 1
        try:
            raw_p = float(r[ci["p"]])
            p = 10 ** (-raw_p) if args.p_is_neglog10 else raw_p
            if not p < args.p_thresh:
                continue
            chrom, pos = norm_chr(r[ci["chr"]]), int(float(r[ci["pos"]]))
            ea, oa = r[ci["ea"]].upper(), r[ci["oa"]].upper()
            eff = float(r[ci["eff"]])
            beta = math.log(eff) if args.or_col else eff
            freq = float(r[ci["freq"]]) if "freq" in ci else None
        except (ValueError, IndexError):
            counts["malformed (skipped)"] += 1
            continue
        counts["p < threshold"] += 1
        if len(ea) != 1 or len(oa) != 1 or not {ea, oa} <= BASES or ea == oa:
            counts["not single-base A/C/G/T"] += 1
            continue
        if freq is not None and not (args.min_freq <= freq <= 1 - args.min_freq):
            counts["frequency outside range"] += 1
            continue
        if not args.keep_sex_chr and not (chrom.isdigit() and 1 <= int(chrom) <= 22):
            counts["not autosomal"] += 1
            continue
        if chrom == mhc_chr and mhc_lo <= pos <= mhc_hi:
            counts["in MHC"] += 1
            continue
        rid = r[ci["rsid"]] if "rsid" in ci else ""
        if not re.fullmatch(r"rs\d+", rid or ""):
            rid = f"{chrom}:{pos}:{ea}:{oa}"
        cands.append({"rsid": rid, "chr": chrom, "pos": pos, "a1": ea, "a2": oa,
                      "beta": beta, "p": p, "freq": freq})

    # Greedy distance clumping: smallest p first (ties: chromosome, position).
    cands.sort(key=lambda d: (d["p"], d["chr"], d["pos"]))
    kept, by_chr = [], {}
    for d in cands:
        near = by_chr.setdefault(d["chr"], [])
        if any(abs(d["pos"] - q) <= args.window for q in near):
            continue
        near.append(d["pos"])
        kept.append(d)
    kept.sort(key=lambda d: (int(d["chr"]) if d["chr"].isdigit() else 99, d["chr"], d["pos"]))
    if not kept:
        sys.exit("No leads passed the filters. Check the column map and --p-thresh.")

    fields = ["rsid", "chr", "pos", "a1", "a2"] + (["a1_freq"] if args.freq_col else []) + ["beta", "pval", "trait"]

    def write(path, rows):
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
            w.writeheader()
            for d in rows:
                w.writerow({**d, "beta": f"{d['beta']:.6g}", "pval": f"{d['p']:.4g}",
                            "a1_freq": "" if d["freq"] is None else f"{d['freq']:.4g}", "trait": args.trait})

    write(out, kept)
    sample = None
    if args.sample:
        if args.sample > len(kept):
            sys.exit(f"--sample {args.sample} is larger than the {len(kept)} leads.")
        sample = random.Random(args.seed).sample(kept, args.sample)
        sample.sort(key=lambda d: (int(d["chr"]) if d["chr"].isdigit() else 99, d["chr"], d["pos"]))
        write(Path(args.sample_out), sample)

    n_pos = sum(d["beta"] > 0 for d in kept)
    lines = [
        f"select_leads.py  {datetime.now().isoformat(timespec='seconds')}",
        f"input   : {gwas}  (md5 {md5(gwas)})",
        f"trait   : {args.trait}",
        f"columns : {cols}",
        f"rules   : p < {args.p_thresh:g}; SNVs only; "
        f"freq in [{args.min_freq}, {1 - args.min_freq}] {'(applied)' if args.freq_col else '(no freq column: not applied)'}; "
        f"{'all chromosomes' if args.keep_sex_chr else 'autosomes only'}; MHC {args.mhc} excluded; greedy clumping, {args.window:,} bp",
        *[f"  {k:<26} {v:,}" for k, v in counts.items()],
        f"  {'candidates after filters':<26} {len(cands):,}",
        f"  {'independent leads':<26} {len(kept):,}  (a1 is risk-raising for {n_pos}; protective for {len(kept) - n_pos})",
        f"output  : {out}",
    ]
    if sample is not None:
        lines.append(f"sample  : {args.sample} leads, seed {args.seed} -> {args.sample_out}")
        lines.append("          " + ", ".join(d["rsid"] for d in sample))
    if args.compare:
        with open(args.compare, newline="") as f:
            ref = {row["rsid"] for row in csv.DictReader(f, delimiter="\t")}
        mine = {d["rsid"] for d in kept}
        lines += [f"compare : {args.compare}",
                  f"  in both {len(mine & ref)} | only here {len(mine - ref)} | only in reference {len(ref - mine)}",
                  f"  only here: {sorted(mine - ref)[:20]}",
                  f"  only in reference: {sorted(ref - mine)[:20]}"]
    log = out.with_suffix(".log.txt")
    log.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"log     : {log}")


if __name__ == "__main__":
    main()