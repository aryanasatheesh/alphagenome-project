#!/usr/bin/env python3
"""Add distance-to-gene features to a SNP table, from a local GENCODE GTF.

Why: the matched-SNP null (Choo, Nov 18) bins SNPs by distance to the nearest
protein-coding gene and by European MAF, then normalises each risk SNP by
non-associated SNPs from the same bin. This script supplies the distance half.
It reads a downloaded GTF instead of a web service, so it gives the same
answer every time and scales to millions of pool candidates (annotate_snps.py
asks Ensembl about each SNP, which is fine for 173 leads and impossible for a
genome-wide pool).

Columns added (input columns are kept unchanged, in order, in front):
  pc_body_dist_bp       0 if the SNP lies inside a protein-coding gene body
                        (gene start..end, any transcript), else bp to the
                        nearest gene edge. Always >= 0. Main matching feature.
  nearest_pc_body_gene  that gene. If several genes contain the SNP, the one
                        whose TSS is nearest.
  n_pc_genes_overlapping  how many protein-coding gene bodies contain the SNP
  pc_tss_dist_bp        SNP position minus the TSS of the protein-coding gene
                        whose TSS is nearest (gene-level TSS: gene start on +,
                        gene end on -). Negative = SNP to the left on the
                        chromosome. Same definition and sign as
                        nearest_pc_tss_distance_bp in annotate_snps.py.
  nearest_pc_tss_gene   that gene
SNPs on a chromosome with no protein-coding genes in the GTF get blanks.

Decisions (also written to the log):
  - Genes: gene lines with gene_type "protein_coding" (--gene-type).
  - Readthrough genes are excluded by default. They are protein_coding loci
    that span two neighbouring genes, so including them would put intergenic
    SNPs between those two genes "inside a gene". GENCODE's own protein-coding
    gene count excludes them. A gene counts as readthrough if its gene line is
    tagged readthrough_gene, or every one of its transcripts is tagged
    readthrough_transcript. --keep-readthrough keeps them (Ensembl REST, which
    annotate_snps.py uses, keeps them).
  - Genes with IDs ending _PAR_Y (chrY copies of pseudoautosomal genes) are
    dropped.
  - Chromosome names are compared without a "chr" prefix, so 7 and chr7 match.
  - Positions are 1-based, like the GTF and the lead tables.
  - Ties (two genes equally near) go to the gene on the left.

Build: the GTF must be in the same build as the positions. --build (default
GRCh38) is checked against the GTF's ##description line and the run stops if
it does not appear there. For hg19 positions use a GENCODE "lift37" GTF with
--build GRCh37, or (preferred, so leads and pool share one build) lift the
positions to hg38 first.

Input is read twice (positions first, then rows streamed to the output), so
memory holds only positions, not the whole table. Gzipped input and output are
fine (.gz).

Usage (from the repo root)
  python3 phase1/add_matching_features.py phase1/data/scz_pgc3_lead_snps_hg38.tsv \\
      --gtf ~/reference/gencode/gencode.v50.annotation.gtf.gz \\
      --compare-annotation phase1/data/scz_pgc3_lead_snps_hg38.annotation.tsv
  -> phase1/data/scz_pgc3_lead_snps_hg38.genes.tsv  (+ .genes.log.txt)

  Options: --out PATH  --chr-col chr  --pos-col pos  --build GRCh38
           --gene-type protein_coding  --keep-readthrough

Requires: numpy.
"""

import argparse
import csv
import gzip
import hashlib
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

ADDED = ["pc_body_dist_bp", "nearest_pc_body_gene", "n_pc_genes_overlapping",
         "pc_tss_dist_bp", "nearest_pc_tss_gene"]
ATTR = re.compile(r'(\S+) "([^"]*)"')


def open_text(path, mode="rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return open(path, mode, encoding="utf-8", newline="")


def norm_chr(c):
    c = c.strip()
    return c[3:] if c.lower().startswith("chr") else c


def file_md5(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


# ---------------------------------------------------------------- GTF
def read_gtf(path, gene_type, keep_readthrough):
    """Return (genes by chromosome, header description, counts)."""
    header = []
    genes = {}                 # gene_id -> [chr, start, end, strand, name, tagged_readthrough]
    tx_total, tx_readthrough = {}, {}
    with open_text(path) as f:
        for line in f:
            if line.startswith("#"):
                header.append(line.rstrip("\n"))
                continue
            if "\tgene\t" not in line and "\ttranscript\t" not in line:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in ("gene", "transcript"):
                continue
            attrs = fields[8]
            gid = ATTR.search(attrs)
            if not gid or gid.group(1) != "gene_id":
                kv = dict(ATTR.findall(attrs))
                gid_val = kv.get("gene_id")
            else:
                gid_val = gid.group(2)
            if gid_val is None:
                continue
            if fields[2] == "gene":
                kv = {}
                tags = []
                for k, v in ATTR.findall(attrs):
                    if k == "tag":
                        tags.append(v)
                    else:
                        kv[k] = v
                if kv.get("gene_type", kv.get("gene_biotype")) != gene_type:
                    continue
                genes[gid_val] = [norm_chr(fields[0]), int(fields[3]), int(fields[4]), fields[6],
                                  kv.get("gene_name", gid_val), "readthrough_gene" in tags]
            elif gid_val in genes:       # transcripts follow their gene line in GENCODE GTFs
                tx_total[gid_val] = tx_total.get(gid_val, 0) + 1
                if 'tag "readthrough_transcript"' in attrs:
                    tx_readthrough[gid_val] = tx_readthrough.get(gid_val, 0) + 1

    counts = {"genes_of_type": len(genes), "readthrough": 0, "par_y": 0, "kept": 0}
    by_chr = {}
    for gid, (c, s, e, strand, name, tagged) in genes.items():
        if gid.endswith("_PAR_Y"):
            counts["par_y"] += 1
            continue
        rt = tagged or (tx_total.get(gid, 0) > 0 and tx_readthrough.get(gid, 0) == tx_total[gid])
        if rt:
            counts["readthrough"] += 1
            if not keep_readthrough:
                continue
        tss = s if strand == "+" else e
        by_chr.setdefault(c, []).append((s, e, tss, name))
        counts["kept"] += 1
    out = {}
    for c, lst in by_chr.items():
        lst.sort(key=lambda g: (g[0], g[1], g[3]))
        out[c] = {"start": np.array([g[0] for g in lst], dtype=np.int64),
                  "end": np.array([g[1] for g in lst], dtype=np.int64),
                  "tss": np.array([g[2] for g in lst], dtype=np.int64),
                  "name": np.array([g[3] for g in lst], dtype=object)}
    desc = next((h for h in header if h.startswith("##description")), "")
    return out, desc, header, counts


# ---------------------------------------------------------------- distances
def features_one_chr(pos, g):
    """pos: int64 array of 1-based positions. g: gene arrays for this chromosome."""
    n = len(pos)
    start, end, tss, name = g["start"], g["end"], g["tss"], g["name"]
    big = np.iinfo(np.int64).max

    # inside a gene body: loop over genes, SNPs sorted (each gene touches only its SNPs)
    order = np.argsort(pos, kind="mergesort")
    ps = pos[order]
    n_over = np.zeros(n, dtype=np.int64)
    best_d = np.full(n, big, dtype=np.int64)
    best_g = np.full(n, -1, dtype=np.int64)
    lo_all = np.searchsorted(ps, start, side="left")
    hi_all = np.searchsorted(ps, end, side="right")
    for gi in np.nonzero(hi_all > lo_all)[0]:
        lo, hi = lo_all[gi], hi_all[gi]
        n_over[lo:hi] += 1
        d = np.abs(ps[lo:hi] - tss[gi])
        better = d < best_d[lo:hi]                     # strict: ties keep the earlier (left) gene
        best_d[lo:hi][better] = d[better]
        best_g[lo:hi][better] = gi
    inv = np.empty(n, dtype=np.int64)
    inv[order] = np.arange(n)
    n_over, best_g = n_over[inv], best_g[inv]

    # outside: nearest edge on the left (largest end < pos) and right (smallest start > pos)
    by_end = np.argsort(end, kind="mergesort")
    end_sorted = end[by_end]
    il = np.searchsorted(end_sorted, pos, side="left") - 1
    has_l = il >= 0
    gl = np.where(has_l, by_end[np.clip(il, 0, None)], -1)
    dl = np.where(has_l, pos - end[np.clip(gl, 0, None)], big)
    ir = np.searchsorted(start, pos, side="right")      # start is sorted
    has_r = ir < len(start)
    gr = np.where(has_r, np.clip(ir, None, len(start) - 1), -1)
    dr = np.where(has_r, start[np.clip(gr, 0, None)] - pos, big)
    left_wins = dl <= dr
    body_d = np.where(n_over > 0, 0, np.where(left_wins, dl, dr))
    body_g = np.where(n_over > 0, best_g, np.where(left_wins, gl, gr))

    # nearest TSS
    t_order = np.argsort(tss, kind="mergesort")
    t_sorted = tss[t_order]
    j = np.searchsorted(t_sorted, pos, side="left")
    jl = np.clip(j - 1, 0, len(t_sorted) - 1)
    jr = np.clip(j, 0, len(t_sorted) - 1)
    d_l = np.where(j - 1 >= 0, pos - t_sorted[jl], big)          # >= 0
    d_r = np.where(j < len(t_sorted), t_sorted[jr] - pos, big)   # >= 0
    use_l = d_l <= d_r
    tg = np.where(use_l, t_order[jl], t_order[jr])
    tss_signed = pos - tss[tg]
    return body_d, name[body_g], n_over, tss_signed, name[tg]


# ---------------------------------------------------------------- annotation cross-check
def compare_annotation(path, rows_feat, log):
    ann = {}
    with open_text(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            ann[r["rsid"]] = r
    n = same_gene = same_dist = body_ok = body_n = 0
    mism = []
    for rsid, feat in rows_feat.items():
        a = ann.get(rsid)
        if not a or not a.get("nearest_pc_gene"):
            continue
        n += 1
        g_ok = a["nearest_pc_gene"] == feat["nearest_pc_tss_gene"]
        d_ok = str(a["nearest_pc_tss_distance_bp"]) == str(feat["pc_tss_dist_bp"])
        same_gene += g_ok
        same_dist += d_ok
        if not (g_ok and d_ok):
            mism.append(f"  {rsid}: Ensembl {a['nearest_pc_gene']} {a['nearest_pc_tss_distance_bp']}"
                        f"  | GTF {feat['nearest_pc_tss_gene']} {feat['pc_tss_dist_bp']}")
        if feat["pc_body_dist_bp"] == "0":
            body_n += 1
            body_ok += feat["nearest_pc_body_gene"] in a.get("overlapping_genes", "").split(",")
    log(f"Cross-check against {path} (Ensembl REST, +/-500 kb window):")
    log(f"  SNPs with an Ensembl nearest protein-coding gene: {n}")
    log(f"  same nearest-TSS gene: {same_gene}/{n}   same signed TSS distance: {same_dist}/{n}")
    log(f"  inside a gene here, and that gene is in Ensembl's overlapping_genes: {body_ok}/{body_n}")
    if mism:
        log(f"  differences ({len(mism)}):")
        for m in mism:
            log(m)


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description="Distance to the nearest protein-coding gene, from a GENCODE GTF.")
    ap.add_argument("snp_file")
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--out")
    ap.add_argument("--chr-col", default="chr")
    ap.add_argument("--pos-col", default="pos")
    ap.add_argument("--build", default="GRCh38", help="must appear in the GTF ##description line (default GRCh38)")
    ap.add_argument("--gene-type", default="protein_coding")
    ap.add_argument("--keep-readthrough", action="store_true")
    ap.add_argument("--compare-annotation", help="annotate_snps.py output for the same SNPs (needs an rsid column)")
    args = ap.parse_args()

    snp = Path(args.snp_file)
    gtf = Path(args.gtf).expanduser()
    if not snp.exists():
        sys.exit(f"SNP file not found: {snp}")
    if not gtf.exists():
        sys.exit(f"GTF not found: {gtf}")
    if args.out:
        out = Path(args.out)
    else:
        base = snp.name[:-3] if snp.name.endswith(".gz") else snp.name
        stem = base[:-4] if base.endswith((".tsv", ".txt")) else base
        out = snp.with_name(stem + ".genes.tsv" + (".gz" if snp.name.endswith(".gz") else ""))
    if out.resolve() == snp.resolve():
        sys.exit("--out would overwrite the input; choose another name")
    log_path = out.with_name(out.name.replace(".gz", "").rsplit(".tsv", 1)[0] + ".log.txt")

    lines = []

    def log(msg=""):
        print(msg)
        lines.append(msg)

    log(f"add_matching_features.py  {datetime.now().isoformat(timespec='seconds')}")
    log(f"SNP file: {snp}")
    log(f"GTF: {gtf}  ({gtf.stat().st_size:,} bytes, md5 {file_md5(gtf)})")

    genes, desc, header, counts = read_gtf(gtf, args.gene_type, args.keep_readthrough)
    log(f"GTF description: {desc or '(none found)'}")
    for h in header:
        if h.startswith(("##date", "##provider", "##format")):
            log(f"  {h}")
    if args.build not in desc:
        sys.exit(f"Build check failed: '{args.build}' is not in the GTF description line.\n"
                 f"  {desc or '(no ##description line)'}\n"
                 f"Use a GTF in the same build as the positions, or set --build if this is intended.")
    log(f"Build check: '{args.build}' found in the description")
    log(f"Genes with gene_type {args.gene_type}: {counts['genes_of_type']:,}; "
        f"readthrough {counts['readthrough']:,} ({'kept' if args.keep_readthrough else 'excluded'}); "
        f"_PAR_Y copies dropped {counts['par_y']:,}; used {counts['kept']:,} on {len(genes)} chromosomes")
    if counts["kept"] == 0:
        sys.exit("No genes left after filtering; check --gene-type and the GTF")

    # pass 1: positions
    chrs, poss = [], []
    with open_text(snp) as f:
        rd = csv.reader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
        head = next(rd)
        for c in (args.chr_col, args.pos_col):
            if c not in head:
                sys.exit(f"{snp} has no '{c}' column (columns: {', '.join(head)})")
        clash = [c for c in ADDED if c in head]
        if clash:
            sys.exit(f"{snp} already has {', '.join(clash)}; was it already annotated?")
        ci, pi = head.index(args.chr_col), head.index(args.pos_col)
        for r in rd:
            if not r:
                continue
            chrs.append(norm_chr(r[ci]))
            poss.append(int(r[pi]))
    n = len(poss)
    pos = np.array(poss, dtype=np.int64)
    chr_arr = np.array(chrs, dtype=object)
    del chrs, poss
    log(f"SNPs: {n:,}")

    body_d = np.full(n, -1, dtype=np.int64)
    body_g = np.full(n, "", dtype=object)
    n_over = np.full(n, -1, dtype=np.int64)
    tss_d = np.zeros(n, dtype=np.int64)
    tss_g = np.full(n, "", dtype=object)
    missing = {}
    for c in sorted(set(chr_arr.tolist())):
        idx = np.nonzero(chr_arr == c)[0]
        if c not in genes:
            missing[c] = len(idx)
            continue
        bd, bg, no, td, tg = features_one_chr(pos[idx], genes[c])
        body_d[idx], body_g[idx], n_over[idx], tss_d[idx], tss_g[idx] = bd, bg, no, td, tg
    if missing:
        log(f"No {args.gene_type} genes on: " + ", ".join(f"{c} ({k} SNPs)" for c, k in missing.items())
            + " -> blank features")
    ok = body_d >= 0
    if ok.any():
        q = np.percentile(body_d[ok], [25, 50, 75, 90])
        log(f"pc_body_dist_bp: {np.mean(body_d[ok] == 0):.1%} inside a gene; quartiles "
            f"{int(q[0]):,} / {int(q[1]):,} / {int(q[2]):,}; 90th pct {int(q[3]):,}")
        log(f"|pc_tss_dist_bp| median {int(np.median(np.abs(tss_d[ok]))):,}")

    # pass 2: stream rows out
    feats_by_rsid = {}
    want_cmp = bool(args.compare_annotation)
    with open_text(snp) as f, open_text(out, "wt") as g:
        rd = csv.reader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
        wr = csv.writer(g, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE, escapechar="\\")
        head = next(rd)
        wr.writerow(head + ADDED)
        rs_i = head.index("rsid") if "rsid" in head else None
        i = -1
        for r in rd:
            if not r:
                continue
            i += 1
            if ok[i]:
                extra = [str(body_d[i]), body_g[i], str(n_over[i]), str(tss_d[i]), tss_g[i]]
            else:
                extra = [""] * 5
            wr.writerow(r + extra)
            if want_cmp and rs_i is not None:
                feats_by_rsid[r[rs_i]] = dict(zip(ADDED, extra))
    if i + 1 != n:
        sys.exit(f"Row count changed between passes ({n} then {i + 1}); was the input file modified?")
    log(f"Written: {out}")

    if want_cmp:
        if not feats_by_rsid:
            log("Cross-check skipped: the SNP file has no rsid column")
        else:
            compare_annotation(args.compare_annotation, feats_by_rsid, log)

    log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Log: {log_path}")


if __name__ == "__main__":
    main()
