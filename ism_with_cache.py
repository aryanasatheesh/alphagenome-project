#!/usr/bin/env python3
"""
ism_with_cache.py
-----------------
Run run_ism_pilot.py with the hg38 sequence windows read from a local cache
instead of the Ensembl website.

Why
    GPU compute nodes on Hoffman2 may not be able to reach the internet, and
    hundreds of Ensembl requests in one job are slow and can hit rate limits.
    So the 131,072 bp windows are downloaded ONCE (on the login node, which has
    internet) and saved; the GPU job then reads them from disk. Everything else
    (allele checks, model, scoring, output files) is run_ism_pilot.main()
    unchanged; only the source of the sequence differs. If a window is missing
    from the cache, it falls back to Ensembl and saves it.

Cache
    data/seq_cache/chr<chrom>_<start>_<end>.txt.gz   (1-based inclusive hg38
    coordinates, exactly the window run_ism_pilot requests; ~35 KB each; *.gz is
    gitignored). Saving the sequences also documents exactly what the model saw.

Usage (repo root, alphagenome_env active)
    python ism_with_cache.py prefetch phase1/data/scz_pgc3_lead_snps_hg38.tsv [more files]
    python ism_with_cache.py run phase1/data/lead_chunks/leads_01.tsv leads_01_ism_scores.tsv
"""

import csv
import gzip
import os
import sys
import time

CACHE = os.path.join("data", "seq_cache")
SEQ_LEN = 131_072
CENTER = SEQ_LEN // 2


def cache_path(chrom, start, end):
    return os.path.join(CACHE, f"chr{chrom}_{start}_{end}.txt.gz")


def load_ism(argv):
    sys.argv = argv                       # run_ism_pilot reads sys.argv when imported
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import run_ism_pilot
    return run_ism_pilot


def main():
    if len(sys.argv) < 3 or sys.argv[1] not in ("prefetch", "run"):
        sys.exit(__doc__)
    mode, args = sys.argv[1], sys.argv[2:]          # keep: load_ism() overwrites sys.argv
    os.makedirs(CACHE, exist_ok=True)

    if mode == "prefetch":
        ism = load_ism(["run_ism_pilot.py"])
        n_new = n_have = 0
        for path in args:
            for r in csv.DictReader(open(path, newline=""), delimiter="\t"):
                chrom = r["chr"][3:] if r["chr"].lower().startswith("chr") else r["chr"]
                pos = int(r["pos"])
                start, end = pos - CENTER, pos + CENTER - 1      # same window as run_ism_pilot
                p = cache_path(chrom, start, end)
                if os.path.exists(p):
                    n_have += 1
                    continue
                seq = ism.fetch_hg38_sequence(chrom, start, end)
                if len(seq) != SEQ_LEN or set(seq) - set("ACGTN"):
                    print(f"  WARNING {r['rsid']}: got {len(seq)} bp / unexpected letters; not cached")
                    continue
                with gzip.open(p + ".tmp", "wt") as f:
                    f.write(seq)
                os.replace(p + ".tmp", p)
                n_new += 1
                print(f"  cached {r['rsid']:<12} chr{chrom}:{start}-{end}")
                time.sleep(0.1)
        print(f"\nPrefetch done: {n_new} downloaded, {n_have} already cached, in {CACHE}/")
        return 0

    # mode == "run"
    if len(args) != 2:
        sys.exit("usage: python ism_with_cache.py run SNP_FILE OUT_NAME")
    ism = load_ism(["run_ism_pilot.py", args[0], args[1]])
    from_network = ism.fetch_hg38_sequence
    stats = {"cache": 0, "network": 0}

    def cached_fetch(chrom, start, end):
        p = cache_path(chrom, start, end)
        if os.path.exists(p):
            stats["cache"] += 1
            with gzip.open(p, "rt") as f:
                return f.read()
        stats["network"] += 1
        print("    (not in cache: fetching from Ensembl)")
        seq = from_network(chrom, start, end)
        if len(seq) == SEQ_LEN:
            with gzip.open(p, "wt") as f:
                f.write(seq)
        return seq

    ism.fetch_hg38_sequence = cached_fetch        # main() looks this name up at call time
    rc = ism.main()
    print(f"Sequences: {stats['cache']} from cache, {stats['network']} from Ensembl")
    return rc


if __name__ == "__main__":
    sys.exit(main())