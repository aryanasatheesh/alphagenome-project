#!/usr/bin/env python3
"""Add European minor allele frequency (MAF) from 1000 Genomes to a SNP table.

Why: the matched-SNP null (Choo, Nov 18) bins SNPs by European MAF and by
distance to the nearest protein-coding gene. add_matching_features.py adds the
distance; this adds the MAF, from one reference panel for every trait, so the
same bins can be reused across GWAS inputs (each GWAS reports frequency in its
own cohorts, or not at all).

Source: the 1000 Genomes high-coverage (30x) phased panel, GRCh38, one VCF per
chromosome (Choo's copy on Hoffman2). Its INFO field already holds European
counts from the 2,504 unrelated samples (503 Europeans, AN_EUR_unrel = 1006),
so no sample list is needed:
    AF_EUR_unrel, AC_EUR_unrel, AN_EUR_unrel, HWE_EUR, ExcHet_EUR
Checked Nov 24: contigs are named chr1..chr22; rs2007044 is at chr12:2,235,794
(hg38) as 12:2235794:A:G, and nothing is at its hg19 position.

Matching: by chromosome, position AND allele set. A panel record matches when
{REF, ALT} equals {a1, a2} (either order). Indels and structural variants at
the same position never match. The a1 frequency is AC/AN if a1 is ALT, else
1 - AC/AN. MAF = min(f, 1 - f), which does not depend on strand, so the MAF is
right even for palindromic SNPs (A/T, C/G); only eur_af_a1 assumes the GWAS
alleles are on the + strand, and the comparison below flags palindromic SNPs
whose frequency fits the other strand better.

Columns added (input columns kept in order, in front):
  eur_maf          minor allele frequency in 1000G EUR (unrelated). Matching feature.
  eur_af_a1        frequency of a1 (the GWAS effect allele) in 1000G EUR
  eur_an           EUR allele number at the site (1006 when every call is present)
  eur_hwe, eur_exchet  panel's HWE and excess-heterozygosity flags (1 = good)
  kg_id kg_ref kg_alt  the matched panel record
  eur_maf_status   matched | not_in_panel | allele_mismatch | multiple_matches
Unmatched SNPs get blank frequencies; the log lists what the panel has there.

If the table has a GWAS frequency column (a1_freq_controls, a1_freq, or
--freq-col), the log compares it with eur_af_a1 (correlation, median and
largest differences, palindromic flags).

Runs on Hoffman2 (where the panel is). Needs bcftools on the PATH:
    module load bcftools/1.21
About one second per chromosome for a few hundred SNPs (indexed lookups).
For a genome-wide pool use a sites-only extraction instead (Step 6).

Usage (from the repo root)
  python3 phase1/add_eur_maf.py phase1/data/scz_pgc3_lead_snps_hg38.genes.tsv \\
      --out phase1/data/scz_pgc3_lead_snps_hg38.features.tsv
  Options: --kg-dir DIR  --kg-pattern '1kGP_high_coverage_Illumina.{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz'
           --chr-col chr --pos-col pos --a1-col a1 --a2-col a2 --freq-col NAME
           --bcftools PATH

Requires: Python 3 standard library and bcftools.
"""

import argparse
import csv
import math
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path

KG_DIR = "/u/project/cluo/chliu/Reference/BB1000Genomes"
KG_PATTERN = "1kGP_high_coverage_Illumina.{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
EXPECTED_AN = 1006          # 503 unrelated Europeans x 2
FIELDS = ["AF_EUR_unrel", "AC_EUR_unrel", "AN_EUR_unrel", "HWE_EUR", "ExcHet_EUR"]
QUERY_FMT = "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t" + "\\t".join(f"%INFO/{f}" for f in FIELDS) + "\\n"
ADDED = ["eur_maf", "eur_af_a1", "eur_an", "eur_hwe", "eur_exchet",
         "kg_id", "kg_ref", "kg_alt", "eur_maf_status"]
PALINDROMIC = {frozenset("AT"), frozenset("CG")}


def norm_chr(c):
    c = c.strip()
    return c[3:] if c.lower().startswith("chr") else c


def fmt(x, nd=6):
    return "" if x is None else f"{x:.{nd}f}"


def run(cmd):
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"Command failed ({p.returncode}): {' '.join(cmd)}\n{p.stderr.strip()}")
    return p.stdout


def pearson(x, y):
    n = len(x)
    if n < 3:
        return float("nan")
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    return sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else float("nan")


def median(v):
    s = sorted(v)
    n = len(s)
    return float("nan") if n == 0 else (s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2)


def contig_name(bcftools, vcf, chrom):
    """The contig name this file uses for chrom ('chr12' or '12')."""
    head = run([bcftools, "view", "-h", str(vcf)])
    names = {line.split("ID=", 1)[1].split(",", 1)[0].rstrip(">")
             for line in head.splitlines() if line.startswith("##contig=<ID=")}
    for cand in (f"chr{chrom}", chrom):
        if cand in names:
            return cand
    sys.exit(f"{vcf}: no contig for chromosome {chrom} in the header")


def query_chr(bcftools, vcf, contig, positions, tmpdir):
    reg = Path(tmpdir) / f"regions_{contig}.tsv"
    with open(reg, "w") as f:
        for p in sorted(set(positions)):
            f.write(f"{contig}\t{p}\t{p}\n")
    out = run([bcftools, "query", "-R", str(reg), "-f", QUERY_FMT, str(vcf)])
    recs = {}
    seen = set()
    for line in out.splitlines():
        if not line:
            continue
        c, pos, vid, ref, alt, af, ac, an, hwe, exc = line.split("\t")
        key = (int(pos), ref, alt)
        if key in seen:                        # -R can return a record once per overlapping region
            continue
        seen.add(key)
        alts = alt.split(",")
        acs = ac.split(",") if ac != "." else ["."] * len(alts)
        hwes = hwe.split(",") if hwe != "." else ["."] * len(alts)
        excs = exc.split(",") if exc != "." else ["."] * len(alts)
        for i, a in enumerate(alts):
            recs.setdefault(int(pos), []).append({
                "id": vid, "ref": ref.upper(), "alt": a.upper(),
                "ac": None if acs[i] == "." else int(acs[i]),
                "an": None if an == "." else int(an),
                "hwe": hwes[i], "exchet": excs[i]})
    return recs


def main():
    ap = argparse.ArgumentParser(description="European MAF from the 1000 Genomes high-coverage panel (GRCh38).")
    ap.add_argument("snp_file")
    ap.add_argument("--out")
    ap.add_argument("--kg-dir", default=KG_DIR)
    ap.add_argument("--kg-pattern", default=KG_PATTERN)
    ap.add_argument("--chr-col", default="chr")
    ap.add_argument("--pos-col", default="pos")
    ap.add_argument("--a1-col", default="a1")
    ap.add_argument("--a2-col", default="a2")
    ap.add_argument("--freq-col", help="GWAS a1 frequency to compare with (default: a1_freq_controls or a1_freq if present)")
    ap.add_argument("--bcftools", default="bcftools")
    args = ap.parse_args()

    snp = Path(args.snp_file)
    if not snp.exists():
        sys.exit(f"SNP file not found: {snp}")
    if args.out:
        out = Path(args.out)
    else:
        stem = snp.name[:-4] if snp.name.endswith(".tsv") else snp.name
        out = snp.with_name(stem + ".maf.tsv")
    if out.resolve() == snp.resolve():
        sys.exit("--out would overwrite the input; choose another name")
    log_path = out.with_name((out.name[:-4] if out.name.endswith(".tsv") else out.name) + ".log.txt")
    kg_dir = Path(args.kg_dir)
    if shutil.which(args.bcftools) is None and not Path(args.bcftools).exists():
        sys.exit(f"bcftools not found ('{args.bcftools}'). On Hoffman2: module load bcftools/1.21")

    lines = []

    def log(msg=""):
        print(msg)
        lines.append(msg)

    with open(snp, newline="", encoding="utf-8") as f:
        rd = csv.reader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
        head = next(rd)
        rows = [r for r in rd if r]
    for c in (args.chr_col, args.pos_col, args.a1_col, args.a2_col):
        if c not in head:
            sys.exit(f"{snp} has no '{c}' column (columns: {', '.join(head)})")
    clash = [c for c in ADDED if c in head]
    if clash:
        sys.exit(f"{snp} already has {', '.join(clash)}; was it already annotated?")
    freq_col = args.freq_col or next((c for c in ("a1_freq_controls", "a1_freq") if c in head), None)
    if args.freq_col and args.freq_col not in head:
        sys.exit(f"--freq-col {args.freq_col} is not a column of {snp}")
    ci, pi, a1i, a2i = (head.index(c) for c in (args.chr_col, args.pos_col, args.a1_col, args.a2_col))
    rs_i = head.index("rsid") if "rsid" in head else None

    bver = run([args.bcftools, "--version"]).splitlines()[0]
    log(f"add_eur_maf.py  {datetime.now().isoformat(timespec='seconds')}")
    log(f"SNP file: {snp}  ({len(rows)} SNPs)")
    log(f"Panel: {kg_dir}/{args.kg_pattern}")
    log(f"Tool: {bver}")

    by_chr = {}
    for i, r in enumerate(rows):
        by_chr.setdefault(norm_chr(r[ci]), []).append(i)

    result = [None] * len(rows)
    with tempfile.TemporaryDirectory() as tmp:
        for chrom in sorted(by_chr, key=lambda c: (not c.isdigit(), int(c) if c.isdigit() else c)):
            idx = by_chr[chrom]
            vcf = kg_dir / args.kg_pattern.format(chr=f"chr{chrom}")
            if not vcf.exists():
                vcf = kg_dir / args.kg_pattern.format(chr=chrom)
            if not vcf.exists():
                log(f"  chr{chrom}: no panel file ({len(idx)} SNPs left blank)")
                for i in idx:
                    result[i] = ("not_in_panel", None, [])
                continue
            st = vcf.stat()
            contig = contig_name(args.bcftools, vcf, chrom)
            recs = query_chr(args.bcftools, vcf, contig, [int(rows[i][pi]) for i in idx], tmp)
            log(f"  {vcf.name}  ({st.st_size:,} bytes, modified "
                f"{datetime.fromtimestamp(st.st_mtime).date()}): {len(idx)} SNPs")
            for i in idx:
                pos = int(rows[i][pi])
                alle = {rows[i][a1i].upper(), rows[i][a2i].upper()}
                here = recs.get(pos, [])
                hits = [h for h in here if {h["ref"], h["alt"]} == alle]
                if len(hits) == 1:
                    result[i] = ("matched", hits[0], here)
                elif len(hits) > 1:
                    result[i] = ("multiple_matches", None, here)
                elif here:
                    result[i] = ("allele_mismatch", None, here)
                else:
                    result[i] = ("not_in_panel", None, here)

    counts = {}
    cmp_rows = []
    unmatched = []
    low_an = []
    out_rows = []
    for i, r in enumerate(rows):
        status, h, here = result[i]
        counts[status] = counts.get(status, 0) + 1
        rid = r[rs_i] if rs_i is not None else f"{r[ci]}:{r[pi]}"
        if status != "matched" or h["ac"] is None or not h["an"]:
            if status == "matched":
                status = "allele_mismatch"      # matched alleles but no EUR counts: treat as unusable
                counts["matched"] -= 1
                counts[status] = counts.get(status, 0) + 1
            unmatched.append(f"  {rid} chr{norm_chr(r[ci])}:{r[pi]} {r[a1i]}/{r[a2i]}: {status}"
                             + (f"; panel has {', '.join(sorted({h2['ref'] + '>' + h2['alt'] for h2 in here}))}"
                                if here else ""))
            out_rows.append(r + ["", "", "", "", "", "", "", "", status])
            continue
        f_alt = h["ac"] / h["an"]
        f_a1 = f_alt if r[a1i].upper() == h["alt"] else 1.0 - f_alt
        maf = min(f_a1, 1.0 - f_a1)
        if h["an"] < 0.9 * EXPECTED_AN:
            low_an.append(f"  {rid}: AN_EUR_unrel {h['an']}")
        out_rows.append(r + [fmt(maf), fmt(f_a1), str(h["an"]), h["hwe"], h["exchet"],
                             h["id"], h["ref"], h["alt"], "matched"])
        if freq_col:
            try:
                g = float(r[head.index(freq_col)])
                pal = frozenset((r[a1i] + r[a2i]).upper()) in PALINDROMIC
                cmp_rows.append((rid, f_a1, g, pal))
            except ValueError:
                pass

    with open(out, "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE, escapechar="\\")
        wr.writerow(head + ADDED)
        wr.writerows(out_rows)

    log("Match status: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items(), key=lambda kv: -kv[1])))
    if unmatched:
        log("Unmatched (blank MAF):")
        for u in unmatched:
            log(u)
    if low_an:
        log(f"Sites with < 90% of EUR calls present ({len(low_an)}):")
        for u in low_an:
            log(u)
    mafs = [float(o[-9]) for o in out_rows if o[-1] == "matched"]
    if mafs:
        ms = sorted(mafs)
        q = [ms[int(p * (len(ms) - 1))] for p in (0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0)]
        log(f"eur_maf: min {q[0]:.3f}; 20/40/50/60/80th pct {q[1]:.3f} / {q[2]:.3f} / {q[3]:.3f} / "
            f"{q[4]:.3f} / {q[5]:.3f}; max {q[6]:.3f}; below 0.01: {sum(m < 0.01 for m in mafs)}")

    if freq_col and cmp_rows:
        d = [abs(a - b) for _, a, b, _ in cmp_rows]
        log(f"Compared with the GWAS column '{freq_col}' ({len(cmp_rows)} SNPs):")
        log(f"  Pearson r (a1 frequency) {pearson([a for _, a, _, _ in cmp_rows], [b for _, _, b, _ in cmp_rows]):.4f}; "
            f"r (MAF) {pearson([min(a, 1 - a) for _, a, _, _ in cmp_rows], [min(b, 1 - b) for _, _, b, _ in cmp_rows]):.4f}")
        log(f"  |difference| median {median(d):.4f}, max {max(d):.4f}; > 0.05: {sum(x > 0.05 for x in d)}; "
            f"> 0.10: {sum(x > 0.10 for x in d)}")
        flagged = [(rid, a, b, pal) for rid, a, b, pal in cmp_rows if abs(a - b) > 0.10]
        for rid, a, b, pal in flagged:
            note = ""
            if pal and abs((1 - a) - b) < abs(a - b):
                note = "  palindromic; fits the other strand better (MAF unaffected)"
            elif pal:
                note = "  palindromic"
            log(f"    {rid}: 1000G EUR {a:.3f} vs GWAS {b:.3f}{note}")
    log(f"Written: {out}")
    log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Log: {log_path}")


if __name__ == "__main__":
    main()
