#!/usr/bin/env python3
"""
annotate_snps.py
----------------
Genomic context for each SNP (Choo, Nov 5: "column of what genes the SNP is next to").

For every SNP in an hg38 SNP file (rsid, chr, pos, ...) this asks Ensembl (GRCh38):
  1. Genes within +/- WINDOW bp      GET /overlap/region/human/{chr}:{start}-{end}?feature=gene
     -> genes the SNP lies inside, the nearest gene (0 bp if inside), the nearest
        protein-coding gene by distance to its transcription start site (TSS), and
        how many protein-coding genes are within the window
  2. Variant consequence by rsID      GET /vep/human/id/{rsid}
     -> most severe consequence (intron_variant, 3_prime_UTR_variant, intergenic_variant, ...)
        and the genes it affects. Why it matters: a coding or UTR variant may act on
        the protein or the RNA, which a chromatin-accessibility score cannot see.

Output: <input stem>.annotation.tsv next to the input, one row per SNP:
    rsid chr pos most_severe_consequence vep_genes overlapping_genes nearest_gene
    nearest_gene_distance_bp nearest_pc_gene nearest_pc_tss_distance_bp n_pc_genes_window
Distances are SNP minus feature (negative = SNP upstream on the + strand).

Network: needs rest.ensembl.org (works from the Mac, like run_ism_pilot.py). About
2 requests per SNP; ~1 minute for 200 SNPs.

Usage (from the repo root)
    python phase1/annotate_snps.py phase1/data/scz_pgc3_lead_snps_hg38.tsv
    python phase1/annotate_snps.py phase1/data/scz_pgc3_control_10snps_hg38.tsv --window 500000
"""

import argparse
import csv
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

SERVER = "https://rest.ensembl.org"


def get_json(path, retries=4):
    req = urllib.request.Request(SERVER + path, headers={"Content-Type": "application/json",
                                                         "User-Agent": "alphagenome-capstone/annotate_snps"})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.loads(resp.read())
        except urllib.error.HTTPError as e:
            if e.code == 429:                        # rate limited: wait as told
                time.sleep(float(e.headers.get("Retry-After", 2)))
                continue
            if e.code in (400, 404):                 # e.g. rsID unknown to VEP
                return None
            if attempt == retries - 1:
                raise
            time.sleep(2)
        except Exception:
            if attempt == retries - 1:
                raise
            time.sleep(2)
    return None


def tss(g):
    return g["start"] if g.get("strand", 1) == 1 else g["end"]


def gene_name(g):
    return g.get("external_name") or g.get("gene_id") or g.get("id")


def annotate(rsid, chrom, pos, window):
    chrom = chrom[3:] if chrom.lower().startswith("chr") else chrom
    genes = get_json(f"/overlap/region/human/{chrom}:{max(1, pos - window)}-{pos + window}"
                     f"?feature=gene;content-type=application/json") or []
    inside = [g for g in genes if g["start"] <= pos <= g["end"]]
    dist = lambda g: 0 if g["start"] <= pos <= g["end"] else (pos - g["end"] if pos > g["end"] else pos - g["start"])
    nearest = min(genes, key=lambda g: (abs(dist(g)), g.get("biotype") != "protein_coding"), default=None)
    pc = [g for g in genes if g.get("biotype") == "protein_coding"]
    nearest_pc = min(pc, key=lambda g: abs(pos - tss(g)), default=None)

    vep = get_json(f"/vep/human/id/{rsid}?content-type=application/json") if rsid.startswith("rs") else None
    cons, vep_genes = "", ""
    if vep:
        v = vep[0]
        cons = v.get("most_severe_consequence", "")
        names = []
        for tc in v.get("transcript_consequences", []):
            if cons in tc.get("consequence_terms", []) and tc.get("gene_symbol") and tc["gene_symbol"] not in names:
                names.append(tc["gene_symbol"])
        vep_genes = ",".join(names)
    return {"most_severe_consequence": cons, "vep_genes": vep_genes,
            "overlapping_genes": ",".join(sorted({gene_name(g) for g in inside})),
            "nearest_gene": gene_name(nearest) if nearest else "",
            "nearest_gene_distance_bp": dist(nearest) if nearest else "",
            "nearest_pc_gene": gene_name(nearest_pc) if nearest_pc else "",
            "nearest_pc_tss_distance_bp": pos - tss(nearest_pc) if nearest_pc else "",
            "n_pc_genes_window": len(pc)}


def main():
    ap = argparse.ArgumentParser(description="Nearest genes and VEP consequence for hg38 SNPs (Ensembl REST).")
    ap.add_argument("snp_file")
    ap.add_argument("--window", type=int, default=500_000, help="bp each side to search for genes")
    ap.add_argument("--out")
    args = ap.parse_args()
    rows = list(csv.DictReader(open(args.snp_file, newline="", encoding="utf-8"), delimiter="\t"))
    for c in ("rsid", "chr", "pos"):
        if rows and c not in rows[0]:
            sys.exit(f"{args.snp_file} has no '{c}' column")
    out = Path(args.out) if args.out else Path(args.snp_file).with_suffix(".annotation.tsv")
    fields = ["rsid", "chr", "pos", "most_severe_consequence", "vep_genes", "overlapping_genes", "nearest_gene",
              "nearest_gene_distance_bp", "nearest_pc_gene", "nearest_pc_tss_distance_bp", "n_pc_genes_window"]
    res = []
    for i, r in enumerate(rows, 1):
        a = annotate(r["rsid"], r["chr"], int(r["pos"]), args.window)
        res.append({"rsid": r["rsid"], "chr": r["chr"], "pos": r["pos"], **a})
        print(f"[{i}/{len(rows)}] {r['rsid']:<12} {a['most_severe_consequence']:<28} "
              f"in: {a['overlapping_genes'] or '-':<20} nearest protein-coding: {a['nearest_pc_gene']} "
              f"({a['nearest_pc_tss_distance_bp']} bp from TSS)")
        time.sleep(0.15)                             # stay well under Ensembl's 15 requests/s
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(res)
    from collections import Counter
    print(f"\nWritten {out} ({len(res)} SNPs)")
    for k, v in Counter(r["most_severe_consequence"] or "unknown" for r in res).most_common():
        print(f"  {v:>4}  {k}")
    return 0


if __name__ == "__main__":
    sys.exit(main())