import gzip
import pandas as pd

# load the file
print("Loading GWAS file...")
df = pd.read_csv(
    "CLOZUK_PGC2noclo.METAL.assoc.dosage.fix.gz",
    sep="\s+",
    compression="gzip"
)

print(f"Total SNPs: {len(df)}")

# filter to genome-wide significant, skip indels (keep single-base alleles only)
sig = df[
    (df["P"] < 5e-8) &
    (df["A1"].str.len() == 1) &
    (df["A2"].str.len() == 1)
].copy()

print(f"Significant SNPs (no indels): {len(sig)}")

# sort by p-value
sig = sig.sort_values("P")

# pick lead SNP per locus: most significant SNP, then exclude anything
# within 1Mb on the same chromosome
lead_snps = []
used_regions = []  # list of (chr, start, end)

for _, row in sig.iterrows():
    chrom = row["CHR"]
    pos = row["BP"]
    # check if this SNP is within 1Mb of an already-selected lead SNP
    too_close = any(
        c == chrom and abs(pos - p) < 1_000_000
        for c, p in used_regions
    )
    if not too_close:
        lead_snps.append(row)
        used_regions.append((chrom, pos))

lead_df = pd.DataFrame(lead_snps)[["SNP", "CHR", "BP", "A2", "A1", "P"]]
lead_df.columns = ["rsid", "chr", "pos", "ref", "alt", "p_value"]

print(f"\nTotal independent loci found: {len(lead_df)}")
print("\nTop 20 lead SNPs:")
print(lead_df.head(20).to_string(index=False))

# save full lead SNP list
lead_df.to_csv("scz_lead_snps.tsv", sep="\t", index=False)
print("\nSaved: scz_lead_snps.tsv")

# save a 5-SNP pilot set (skip chr6 MHC as first pick, take diverse chromosomes)
pilot = lead_df[lead_df["chr"] != 6].head(5)
pilot.to_csv("scz_pilot_5snps.tsv", sep="\t", index=False)
print("\nPilot 5 SNPs (non-MHC):")
print(pilot.to_string(index=False))
