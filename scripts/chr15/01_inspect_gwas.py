#!/usr/bin/env python3
"""
Step 1 QC: inspect the FinnGen DF13 Pneumonia summary statistics.
Run this the moment the file lands. Confirms build, columns, case proportion,
computes MAF from af_alt, and flags ambiguous strand-flip SNPs.

Usage:
    python3 01_inspect_gwas.py finngen_R13_J10_PNEUMONIA.gz
"""
import sys
import gzip
import pandas as pd

# FinnGen total sample size for DF13 — UPDATE from the manifest/release page
# DF13 cases/controls for J10_PNEUMONAS go here once you read the manifest row.
# The manifest TSV has num_cases and num_controls columns per endpoint.

def sniff_header(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
    return header

def main(path):
    print(f"== Inspecting {path} ==\n")
    header = sniff_header(path)
    print("Columns found:")
    for c in header:
        print(f"  - {c}")
    print()

    # Read a sample for fast inspection, then the full file for stats.
    df = pd.read_csv(path, sep="\t", nrows=50000)
    cols = {c.lower().lstrip("#"): c for c in df.columns}

    # --- Build check ---
    chrom_col = cols.get("chrom") or cols.get("chr")
    pos_col   = cols.get("pos")
    ref_col   = cols.get("ref")
    alt_col   = cols.get("alt")
    print("== BUILD CHECK ==")
    print(df[[chrom_col, pos_col, ref_col, alt_col]].head())
    print("If positions match hg38 coordinates (e.g. chr-style or 1-25),")
    print("and you sourced FinnGen DF13, this is GRCh38 — no liftOver needed.\n")

    # --- MAF from af_alt ---
    af_col = cols.get("af_alt") or cols.get("af")
    if af_col:
        af = df[af_col].astype(float)
        maf = af.where(af <= 0.5, 1 - af)
        print("== MAF (derived from af_alt) ==")
        print(maf.describe())
        print(f"SNPs with MAF < 0.01 in sample: {(maf < 0.01).sum()}\n")

    # --- Ambiguous strand SNPs (A/T, C/G) ---
    amb = {("A","T"),("T","A"),("C","G"),("G","C")}
    is_amb = df.apply(lambda r: (str(r[ref_col]).upper(), str(r[alt_col]).upper()) in amb, axis=1)
    print("== AMBIGUOUS STRAND SNPs (A/T, C/G) ==")
    print(f"In 50k sample: {is_amb.sum()} ({100*is_amb.mean():.1f}%)")
    print("These get extra scrutiny at harmonization when MAF is near 0.5.\n")

    # --- p-value sanity ---
    p_col = cols.get("pval") or cols.get("p")
    if p_col:
        print("== P-VALUE RANGE ==")
        print(df[p_col].describe())
        print(f"Genome-wide significant (p<5e-8) in sample: {(df[p_col] < 5e-8).sum()}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1])
