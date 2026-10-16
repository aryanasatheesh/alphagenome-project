#!/usr/bin/env python3
"""
liftover_to_hg38.py
-------------------
Lift SNP coordinates in a TSV file from hg19 (GRCh37) to hg38 (GRCh38).

Why this exists
    PGC3 schizophrenia summary statistics report hg19 coordinates, but
    AlphaGenome was trained on hg38. Every variant must be lifted before
    sequence extraction or in silico mutagenesis.

Design: strict by default
    A SNP is written to the output only if it maps to exactly one place,
    on a standard chromosome, on the same chromosome it came from, on the
    + strand. Anything else is recorded in the report with a reason, and the
    script exits with status 1 WITHOUT writing the output file, so a partial
    result can never be picked up by accident. Pass --allow-failures to write
    the successfully lifted rows anyway (failures stay in the report).

Usage
    python phase1/liftover_to_hg38.py \\
        --input  phase1/data/scz_pgc3_lead_snps.tsv \\
        --output phase1/data/scz_pgc3_lead_snps_hg38.tsv

    Optional: --chain PATH  --chr-col NAME  --pos-col NAME
              --report PATH  --allow-failures

Inputs
    A tab-separated file with a header row. Chromosome may be written as
    "6" or "chr6". Position must be a 1-based integer (the convention used by
    GWAS summary statistics, VCF, dbSNP and Ensembl).

Outputs
    <output>                        input rows with chr/pos replaced by hg38
                                    values, plus chr_hg19 and pos_hg19 columns
    <output stem>.liftover_report.tsv   one row per input SNP with its status
    <output stem>.liftover_log.txt      provenance: command, chain file and its
                                        md5, pyliftover version, counts

Coordinate conventions (the most likely place to silently corrupt data)
    pyliftover works in 0-based coordinates. This script subtracts 1 before
    calling it and adds 1 to the result, so input and output are both 1-based.

Requirements
    pip install pyliftover          (validated with 0.4.1)
    UCSC chain file, downloaded once (~222 KB, gitignored):
    curl -o phase1/data/hg19ToHg38.over.chain.gz \\
      https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

Out of scope
    This script moves coordinates only. It does not touch allele columns.
    Allele orientation (which base is reference in hg38, and which is the GWAS
    effect allele) is checked against the hg38 sequence in run_ism_pilot.py and
    documented in docs/allele_conventions.md.
"""

import argparse
import csv
import datetime
import hashlib
import shlex
import sys
from pathlib import Path

try:
    import pyliftover
    from pyliftover import LiftOver
except ImportError:
    sys.exit("pyliftover is not installed. Run: pip install pyliftover")

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_CHAIN = SCRIPT_DIR / "data" / "hg19ToHg38.over.chain.gz"
CHAIN_URL = ("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/"
             "hg19ToHg38.over.chain.gz")

# Chromosomes we accept on either side. Anything else (chrM, chr6_..._alt,
# chrUn_..., chr6_..._random) is rejected rather than passed downstream.
STANDARD_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y"}

# Every status a SNP can receive. Only "lifted" rows reach the output.
STATUS_REASONS = {
    "lifted":        "mapped uniquely to the same standard chromosome on the + strand",
    "invalid_input": "chromosome or position in the input could not be parsed",
    "unknown_chrom": "chromosome name not present in the chain file",
    "unmapped":      "position falls in a region with no hg38 equivalent",
    "multi":         "position maps to more than one hg38 location",
    "alt_contig":    "maps to a non-standard contig (alt, random, chrUn, chrM)",
    "chrom_change":  "maps to a different chromosome than it came from",
    "strand_flip":   "maps to the - strand; the hg38 base is the complement of the hg19 base",
}

REPORT_FIELDS = ["row", "id", "chr_hg19", "pos_hg19", "status", "n_mappings",
                 "chr_hg38", "pos_hg38", "strand", "note"]


def strip_chr(name):
    """'chr6' -> '6', '6' -> '6'. Case-insensitive on the prefix only."""
    return name[3:] if name.lower().startswith("chr") else name


def lift_one(lo, chrom_raw, pos_raw):
    """
    Lift a single 1-based position.

    Returns a dict with keys: status, n_mappings, chr_hg38, pos_hg38 (1-based,
    int or None), strand, note.
    """
    out = {"n_mappings": "", "chr_hg38": "", "pos_hg38": None,
           "strand": "", "note": ""}

    chrom = strip_chr(chrom_raw.strip())
    if chrom not in STANDARD_CHROMS:
        return {**out, "status": "invalid_input",
                "note": f"chromosome {chrom_raw!r} is not 1-22, X or Y"}
    try:
        pos1 = int(pos_raw.strip())
    except (ValueError, AttributeError):
        return {**out, "status": "invalid_input",
                "note": f"position {pos_raw!r} is not an integer"}
    if pos1 < 1:
        return {**out, "status": "invalid_input",
                "note": f"position {pos1} is not a valid 1-based coordinate"}

    # 1-based -> 0-based for pyliftover
    results = lo.convert_coordinate("chr" + chrom, pos1 - 1)

    # pyliftover returns None for an unknown chromosome, [] for an unmapped
    # position, and otherwise a list sorted by chain score (best first).
    if results is None:
        return {**out, "status": "unknown_chrom"}
    if len(results) == 0:
        return {**out, "status": "unmapped", "n_mappings": 0}
    if len(results) > 1:
        targets = "; ".join(f"{c}:{p + 1}({s})" for c, p, s, _ in results)
        return {**out, "status": "multi", "n_mappings": len(results),
                "note": targets}

    new_chrom_full, new_pos0, strand, _score = results[0]
    new_chrom = strip_chr(new_chrom_full)
    lifted = {**out, "n_mappings": 1, "chr_hg38": new_chrom,
              "pos_hg38": new_pos0 + 1,  # 0-based -> 1-based
              "strand": strand}

    # Order matters: a mapping onto an alt contig is reported as such even if
    # it is also on the - strand, because the contig is the bigger problem.
    if new_chrom not in STANDARD_CHROMS:
        return {**lifted, "status": "alt_contig", "note": new_chrom_full}
    if new_chrom != chrom:
        return {**lifted, "status": "chrom_change",
                "note": f"chr{chrom} -> chr{new_chrom}"}
    if strand != "+":
        return {**lifted, "status": "strand_flip"}
    return {**lifted, "status": "lifted"}


def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(
        description="Lift SNP coordinates in a TSV from hg19 to hg38 (strict by default).")
    ap.add_argument("--input", required=True, help="hg19 TSV with a header row")
    ap.add_argument("--output", required=True, help="hg38 TSV to write")
    ap.add_argument("--chain", default=str(DEFAULT_CHAIN),
                    help=f"UCSC chain file (default: {DEFAULT_CHAIN})")
    ap.add_argument("--chr-col", default="chr", help="chromosome column (default: chr)")
    ap.add_argument("--pos-col", default="pos", help="1-based position column (default: pos)")
    ap.add_argument("--id-col", default="rsid",
                    help="identifier column used in messages (default: rsid; optional)")
    ap.add_argument("--report", default=None,
                    help="per-SNP report TSV (default: <output stem>.liftover_report.tsv)")
    ap.add_argument("--allow-failures", action="store_true",
                    help="write the lifted rows even if some SNPs failed (exit 0)")
    args = ap.parse_args()

    in_path, out_path, chain_path = Path(args.input), Path(args.output), Path(args.chain)
    stem = out_path.with_suffix("")
    report_path = Path(args.report) if args.report else Path(f"{stem}.liftover_report.tsv")
    log_path = Path(f"{stem}.liftover_log.txt")

    # ---- Preconditions -------------------------------------------------------
    if not in_path.exists():
        sys.exit(f"Input file not found: {in_path}")
    if not chain_path.exists():
        sys.exit(f"Chain file not found: {chain_path}\nDownload it with:\n"
                 f"  curl -o {chain_path} {CHAIN_URL}")
    if in_path.resolve() == out_path.resolve():
        sys.exit("Refusing to overwrite the input file. Choose a different --output.")

    with open(in_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)

    for col in (args.chr_col, args.pos_col):
        if col not in fieldnames:
            sys.exit(f"Column {col!r} not found in {in_path}. Columns are: {fieldnames}")
    # Guard against lifting a file twice, which would silently shift every SNP.
    if "pos_hg19" in fieldnames or "chr_hg19" in fieldnames:
        sys.exit(f"{in_path} already has chr_hg19/pos_hg19 columns: it looks like an "
                 "already-lifted file. Lift the original hg19 file instead.")
    if not rows:
        sys.exit(f"{in_path} has a header but no data rows.")
    has_id = args.id_col in fieldnames

    # ---- Lift ----------------------------------------------------------------
    print(f"INPUT : {in_path}  ({len(rows)} SNPs)")
    print(f"CHAIN : {chain_path}")
    lo = LiftOver(str(chain_path))

    # Preserve the input's chromosome style ("6" vs "chr6") in the output.
    uses_prefix = rows[0][args.chr_col].strip().lower().startswith("chr")

    lifted_rows, report_rows = [], []
    for i, row in enumerate(rows, start=1):
        res = lift_one(lo, row[args.chr_col], row[args.pos_col])
        snp_id = row[args.id_col] if has_id else ""
        report_rows.append({
            "row": i, "id": snp_id,
            "chr_hg19": row[args.chr_col], "pos_hg19": row[args.pos_col],
            "status": res["status"], "n_mappings": res["n_mappings"],
            "chr_hg38": res["chr_hg38"],
            "pos_hg38": "" if res["pos_hg38"] is None else res["pos_hg38"],
            "strand": res["strand"], "note": res["note"],
        })
        if res["status"] == "lifted":
            new = dict(row)
            new["chr_hg19"] = row[args.chr_col]
            new["pos_hg19"] = row[args.pos_col]
            new[args.chr_col] = ("chr" if uses_prefix else "") + res["chr_hg38"]
            new[args.pos_col] = str(res["pos_hg38"])
            lifted_rows.append(new)

    # ---- Summarise -----------------------------------------------------------
    counts = {}
    for r in report_rows:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    n_failed = len(rows) - counts.get("lifted", 0)

    print(f"\n{'status':<14} {'n':>5}  meaning")
    for status, reason in STATUS_REASONS.items():
        if counts.get(status):
            print(f"{status:<14} {counts[status]:>5}  {reason}")
    if n_failed:
        print("\nSNPs that did not lift:")
        for r in report_rows:
            if r["status"] != "lifted":
                label = r["id"] or f"row {r['row']}"
                print(f"  {label:<14} chr{strip_chr(r['chr_hg19'])}:{r['pos_hg19']:<12} "
                      f"{r['status']:<13} {r['note']}")

    # The report is always written: it is how failures get inspected.
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=REPORT_FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(report_rows)

    wrote_output = n_failed == 0 or args.allow_failures
    if wrote_output:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames + ["chr_hg19", "pos_hg19"],
                               delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(lifted_rows)

    # Provenance lives in a sidecar file, not in the TSV, because comment lines
    # at the top of a TSV would be read as the header by csv.DictReader.
    with open(log_path, "w") as f:
        f.write(f"date              : {datetime.datetime.now().isoformat(timespec='seconds')}\n")
        f.write(f"command           : {' '.join(shlex.quote(a) for a in sys.argv)}\n")
        f.write(f"input             : {in_path}\n")
        f.write(f"output            : {out_path if wrote_output else '(not written)'}\n")
        f.write(f"chain file        : {chain_path.name}\n")
        f.write(f"chain md5         : {md5sum(chain_path)}\n")
        f.write(f"pyliftover version: {getattr(pyliftover, '__version__', 'unknown')}\n")
        f.write(f"python version    : {sys.version.split()[0]}\n")
        f.write(f"snps in           : {len(rows)}\n")
        for status in STATUS_REASONS:
            if counts.get(status):
                f.write(f"  {status:<16}: {counts[status]}\n")

    print(f"\nREPORT: {report_path}")
    print(f"LOG   : {log_path}")
    if wrote_output:
        print(f"OUTPUT: {out_path}  ({len(lifted_rows)} SNPs)")
        if n_failed:
            print(f"WARNING: {n_failed} SNP(s) excluded (--allow-failures). See the report.")
        return 0
    print(f"\nFAILED: {n_failed} of {len(rows)} SNPs did not lift. Output NOT written.")
    print("Inspect the report, then rerun with --allow-failures to proceed without them.")
    return 1


if __name__ == "__main__":
    sys.exit(main())