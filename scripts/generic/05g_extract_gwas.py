#!/usr/bin/env python3
"""
Extract a FinnGen locus slice into a small reusable file. Generic version.

Lancer depuis la racine Project/. Sortie : results/<locus>/gwas_<locus>.tsv

Usage:
    python3 scripts/generic/05g_extract_gwas.py <gwas_file> <locus_name> <chrom> <start> <end>
Example (CRP, chr1):
    python3 scripts/generic/05g_extract_gwas.py data/finngen_R13_J10_PNEUMONIA.gz CRP 1 159200000 160200000
"""
import sys
import os
import pandas as pd

CHUNK = 500_000


def main(path, locus, chrom, start, end):
    chrom = int(chrom); start = int(start); end = int(end)
    os.makedirs(f"results/{locus}", exist_ok=True)
    out = f"results/{locus}/gwas_{locus}.tsv"
    kept = []
    total = 0
    for chunk in pd.read_csv(path, sep="\t", chunksize=CHUNK):
        total += len(chunk)
        cnum = pd.to_numeric(chunk["#chrom"], errors="coerce")
        sub = chunk[(cnum == chrom) &
                    (chunk["pos"] >= start) & (chunk["pos"] <= end)]
        if len(sub):
            kept.append(sub)
        # arret anticipe : fichier trie par chrom croissant.
        # des qu'un chunk est entierement au-dela du chrom cible, on stoppe.
        if kept and (cnum > chrom).all():
            break
    if not kept:
        print("No rows found in window — check chrom/window."); return
    g = pd.concat(kept, ignore_index=True)
    g.to_csv(out, sep="\t", index=False)
    print(f"Scanned {total:,} rows.")
    print(f"{locus} window chr{chrom}:{start}-{end}: {len(g):,} variants kept")
    print(f"position range: {g['pos'].min()}-{g['pos'].max()}")
    print(f"genome-wide significant here: {(g['pval'] < 5e-8).sum()}")
    print(f"saved: {out}")


if __name__ == "__main__":
    if len(sys.argv) < 6:
        print(__doc__); sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
