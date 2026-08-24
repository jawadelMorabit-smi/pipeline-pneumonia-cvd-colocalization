#!/usr/bin/env python3
"""
Step 2: scan the ENTIRE FinnGen file (not just the first 50k rows)
to find the real pneumonia signals.

Goal: count genome-wide significant variants (p < 5e-8) and list the
strongest loci with their nearest genes. If signals match the validation
papers (MUC5AC, TNFRSF1A, HYKK, PBX3, CHRNA3/5), the file is confirmed good.

Reads in chunks so the large file fits comfortably in memory.

Usage:
    python3 02_scan_full.py finngen_R13_J10_PNEUMONIA.gz
"""
import sys
import pandas as pd

GW_SIG = 5e-8          # genome-wide significance threshold
CHUNK = 500_000        # rows per chunk

def main(path):
    print(f"== Scanning full file: {path} ==\n")

    total_rows = 0
    sig_rows = []      # collect only the significant hits (small set)

    reader = pd.read_csv(path, sep="\t", chunksize=CHUNK)
    for i, chunk in enumerate(reader):
        total_rows += len(chunk)
        hits = chunk[chunk["pval"] < GW_SIG]
        if len(hits):
            sig_rows.append(hits)
        print(f"  chunk {i+1}: read {len(chunk):>7} rows "
              f"(running total {total_rows:,}), "
              f"sig hits so far: {sum(len(h) for h in sig_rows)}")

    print(f"\nTotal variants in file: {total_rows:,}")

    if not sig_rows:
        print("No genome-wide significant hits found. That would be a red flag — "
              "re-check the file.")
        return

    sig = pd.concat(sig_rows, ignore_index=True)
    print(f"Genome-wide significant variants (p < 5e-8): {len(sig):,}\n")

    # How many distinct genomic regions? Count significant hits per chromosome.
    print("== Significant hits per chromosome ==")
    print(sig["#chrom"].value_counts().sort_index().to_string())
    print()

    # Top 20 strongest signals with nearest gene
    cols = ["#chrom", "pos", "ref", "alt", "rsids",
            "nearest_genes", "pval", "beta"]
    top = sig.sort_values("pval").head(20)[cols]
    print("== Top 20 strongest signals ==")
    with pd.option_context("display.max_rows", None,
                           "display.width", 200):
        print(top.to_string(index=False))

    # Flag any validation genes that appear
    targets = ["MUC5AC", "MUC5B", "TNFRSF1A", "HYKK", "PBX3",
               "CHRNA3", "CHRNA5"]
    genes_seen = " ".join(sig["nearest_genes"].dropna().astype(str).tolist())
    print("\n== Validation-gene check ==")
    for g in targets:
        present = g in genes_seen
        mark = "FOUND" if present else "not in nearest_genes"
        print(f"  {g:>9}: {mark}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1])
