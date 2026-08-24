#!/usr/bin/env python3
"""
Extract the FinnGen chr15 slice (locus window) into a small reusable file.
Reads the big FinnGen file in chunks, keeps only chr15 within the window.

Usage:
    python3 05_extract_gwas_chr15.py finngen_R13_J10_PNEUMONIA.gz
"""
import sys
import os
import pandas as pd

os.makedirs("results", exist_ok=True)

CHROM = 15
START = 78000000
END   = 79500000
CHUNK = 500_000
OUT   = "results/gwas_chr15.tsv"


def main(path):
    kept = []
    total = 0
    for chunk in pd.read_csv(path, sep="\t", chunksize=CHUNK):
        total += len(chunk)
        # FinnGen column names: #chrom, pos, ref, alt, ...
        sub = chunk[(chunk["#chrom"] == CHROM) &
                    (chunk["pos"] >= START) & (chunk["pos"] <= END)]
        if len(sub):
            kept.append(sub)
    if not kept:
        print("No chr15 rows found in window — check column names/window.")
        return
    gwas = pd.concat(kept, ignore_index=True)
    gwas.to_csv(OUT, sep="\t", index=False)
    print(f"Scanned {total:,} rows.")
    print(f"chr15 window {START}-{END}: {len(gwas):,} variants kept")
    print(f"position range: {gwas['pos'].min()}-{gwas['pos'].max()}")
    print(f"genome-wide significant here: {(gwas['pval'] < 5e-8).sum()}")
    print(f"saved: {OUT}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    main(sys.argv[1])
