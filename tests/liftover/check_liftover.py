#!/usr/bin/env python3
"""
check_liftover.py
-----------------
Self-test for phase1/liftover_to_hg38.py using the REAL UCSC chain file and
real PGC3 SNPs whose outcome is known: two that lift (with known hg38
positions), and one each that are unmapped, land on an alt contig, or lie in a
region inverted in hg38, plus one invalid input.

Every expected value comes from earlier runs: the PGC3 pilot liftover
(August 2026, checked against Ensembl) and the whole-genome edge-case survey
(October 2026, results/liftover_survey/).

Usage (from the repo root, after downloading the chain file):
    python tests/liftover/check_liftover.py

Options:
    --chain PATH   chain file (default: phase1/data/hg19ToHg38.over.chain.gz)
    --cases PATH   test SNPs (default: tests/liftover/edge_cases_hg19.tsv)

Prints PASS/FAIL for each check and exits 0 only if every check passes.
Writes nothing inside the repo (all outputs go to a temporary folder).
"""

import argparse
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "phase1" / "liftover_to_hg38.py"

results = []


def check(name, ok, detail=""):
    results.append(ok)
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  ({detail})" if detail and not ok else ""))


def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def run(args):
    return subprocess.run([sys.executable, str(SCRIPT)] + args,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          universal_newlines=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--chain", default=str(REPO / "phase1" / "data" / "hg19ToHg38.over.chain.gz"))
    ap.add_argument("--cases", default=str(REPO / "tests" / "liftover" / "edge_cases_hg19.tsv"))
    args = ap.parse_args()

    # The liftover script runs under this same Python, so check it can import
    # pyliftover first; otherwise every check fails with a confusing message.
    try:
        import pyliftover  # noqa: F401
    except ImportError:
        sys.exit(f"pyliftover is not installed for this Python ({sys.executable}).\n"
                 "Activate the environment that has it first, e.g. on the Mac:\n"
                 "    source ~/alphagenome-env/bin/activate\n"
                 "or install it:  pip install pyliftover==0.4.1")

    for p, what in ((SCRIPT, "liftover script"), (Path(args.chain), "chain file"),
                    (Path(args.cases), "test cases")):
        if not p.exists():
            sys.exit(f"Missing {what}: {p}")

    cases = read_tsv(args.cases)
    expected = {c["rsid"]: c for c in cases}
    n_lift = sum(c["expected_status"] == "lifted" for c in cases)

    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "out_hg38.tsv"
        report = Path(tmp) / "out_hg38.liftover_report.tsv"

        print("1. Strict mode (default)")
        r = run(["--input", args.cases, "--output", str(out), "--chain", args.chain])
        # Require the script's own failure summary, so a crash or import error
        # (which also exits 1) cannot pass this check.
        check("exits with status 1 because some SNPs cannot be lifted",
              r.returncode == 1 and "Output NOT written" in r.stdout,
              f"exit {r.returncode}; output:\n{r.stdout}")
        check("does NOT write the output file", report.exists() and not out.exists(),
              "no report either: the script did not run" if not report.exists() else "")
        check("writes the per-SNP report", report.exists())
        if report.exists():
            rep = {row["id"]: row for row in read_tsv(report)}
            check("report has one row per input SNP", len(rep) == len(cases),
                  f"{len(rep)} rows for {len(cases)} SNPs")
            for rsid, c in expected.items():
                got = rep.get(rsid, {})
                check(f"{rsid}: status {c['expected_status']}",
                      got.get("status") == c["expected_status"],
                      f"got {got.get('status')!r}")
                if c["expected_pos_hg38"]:
                    check(f"{rsid}: hg38 position {c['expected_pos_hg38']}",
                          got.get("pos_hg38") == c["expected_pos_hg38"],
                          f"got {got.get('pos_hg38')!r}")

        print("2. --allow-failures")
        r = run(["--input", args.cases, "--output", str(out), "--chain", args.chain,
                 "--allow-failures"])
        check("exits with status 0", r.returncode == 0, f"exit {r.returncode}")
        if out.exists():
            rows = read_tsv(out)
            check(f"output holds only the {n_lift} lifted SNPs", len(rows) == n_lift,
                  f"{len(rows)} rows")
            for row in rows:
                c = expected[row["rsid"]]
                check(f"{row['rsid']}: pos is hg38, pos_hg19 keeps the original",
                      row["pos"] == c["expected_pos_hg38"] and row["pos_hg19"] == c["pos"],
                      f"pos {row['pos']}, pos_hg19 {row['pos_hg19']}")
                check(f"{row['rsid']}: extra input columns carried through",
                      row.get("note") == c["note"])
        else:
            check("writes the output file", False)

        print("3. Guards")
        if out.exists():
            r = run(["--input", str(out), "--output", str(Path(tmp) / "twice.tsv"),
                     "--chain", args.chain])
            check("refuses to lift an already-lifted file",
                  r.returncode == 1 and "already-lifted" in r.stdout, r.stdout.strip())
        r = run(["--input", args.cases, "--output", args.cases, "--chain", args.chain])
        check("refuses to overwrite its input", r.returncode == 1 and "overwrite" in r.stdout,
              r.stdout.strip())
        r = run(["--input", args.cases, "--output", str(out), "--chain", args.chain,
                 "--pos-col", "position"])
        check("names the missing column", r.returncode == 1 and "'position'" in r.stdout,
              r.stdout.strip())
        r = run(["--input", args.cases, "--output", str(out),
                 "--chain", str(Path(tmp) / "missing.chain.gz")])
        check("explains how to get a missing chain file",
              r.returncode == 1 and "curl" in r.stdout, r.stdout.strip())

    n_fail = results.count(False)
    print(f"\n{len(results) - n_fail}/{len(results)} checks passed"
          + ("" if not n_fail else f"  ({n_fail} FAILED)"))
    print("ALL CHECKS PASSED" if not n_fail else "SOME CHECKS FAILED: see FAIL lines above")
    return 0 if not n_fail else 1


if __name__ == "__main__":
    sys.exit(main())