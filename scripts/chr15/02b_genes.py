#!/usr/bin/env python3
"""
Step 2b: re-scan ONLY for the significant hits and report the genes,
without re-reading 21M rows for display. Still streams the file in chunks
but only keeps p < 5e-8 rows, then summarizes.

Usage:
    python3 02b_genes.py finngen_R13_J10_PNEUMONIA.gz
"""
import sys
import pandas as pd

GW_SIG = 5e-8
CHUNK = 1_000_000

def main(path):
    keep = []
    for chunk in pd.read_csv(path, sep="\t", chunksize=CHUNK):
        keep.append(chunk[chunk["pval"] < GW_SIG])
    sig = pd.concat(keep, ignore_index=True)
    print(f"Significant variants: {len(sig)}\n")

    # Genes present at each chromosome's locus
    print("== Genes named at each significant locus ==")
    for chrom, grp in sig.groupby("#chrom"):
        genes = sorted(set(grp["nearest_genes"].dropna().astype(str)))
        # split combined entries like "GENE1,GENE2"
        flat = sorted({g for entry in genes for g in entry.split(",")})
        print(f"  chr{chrom}: {', '.join(flat) if flat else '(none annotated)'}")
    print()

    # Validation-gene check (blanks dropped)
    targets = ["MUC5AC", "MUC5B", "TNFRSF1A", "HYKK", "PBX3",
               "CHRNA3", "CHRNA5"]
    genes_seen = " ".join(sig["nearest_genes"].dropna().astype(str).tolist())
    print("== Validation-gene check ==")
    for g in targets:
        mark = "FOUND" if g in genes_seen else "not in nearest_genes"
        print(f"  {g:>9}: {mark}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    main(sys.argv[1])
