#!/usr/bin/env python3
"""
Harmonize the FinnGen chr15 GWAS slice against a GTEx eQTL table (one tissue),
per gene, producing coloc-ready inputs.

Three checks (the diagram):
  1. Match variants by position + sorted allele pair (order-independent).
  2. Align effect allele: if ref/alt are swapped between GWAS and eQTL,
     flip the eQTL beta sign and use 1-MAF. (eQTL effect allele = ALT.)
  3. Remove ambiguous strand SNPs (A/T, C/G) with MAF near 0.5.

Outputs one harmonized TSV per gene into results/harmonized_<tissue>/.

Usage:
    python3 06_harmonize_chr15.py results/gwas_chr15.tsv results/eqtl_chr15_blood.tsv blood
    python3 06_harmonize_chr15.py results/gwas_chr15.tsv results/eqtl_chr15_lung.tsv  lung
"""
import sys
import os
import pandas as pd

AMBIG = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
MAF_AMBIG_CUTOFF = 0.40   # treat as "near 0.5" if MAF > this


def norm_key(chrom, pos, a1, a2):
    """Order-independent variant key: chr_pos_<sorted alleles>."""
    x, y = sorted([str(a1).upper(), str(a2).upper()])
    return f"{chrom}_{pos}_{x}_{y}"


def maf_from_af(af):
    af = float(af)
    return af if af <= 0.5 else 1 - af


def main(gwas_path, eqtl_path, tissue):
    outdir = f"results/harmonized_{tissue}"
    os.makedirs(outdir, exist_ok=True)

    # --- load GWAS slice (FinnGen columns) ---
    g = pd.read_csv(gwas_path, sep="\t")
    g["key"] = [norm_key(c, p, r, a) for c, p, r, a in
                zip(g["#chrom"], g["pos"], g["ref"], g["alt"])]
    g["gwas_maf"] = g["af_alt"].apply(maf_from_af)
    # FinnGen effect allele = alt; keep ref/alt to resolve orientation
    g = g.rename(columns={"ref": "gwas_ref", "alt": "gwas_alt",
                          "pval": "gwas_p", "beta": "gwas_beta",
                          "sebeta": "gwas_se"})
    gsub = g[["key", "pos", "gwas_ref", "gwas_alt", "gwas_p",
              "gwas_beta", "gwas_se", "gwas_maf", "rsids"]]

    # --- load eQTL (Catalogue columns); effect allele = alt ---
    e = pd.read_csv(eqtl_path, sep="\t")
    e["key"] = [norm_key(c, p, r, a) for c, p, r, a in
                zip(e["chromosome"], e["position"], e["ref"], e["alt"])]

    genes = sorted(e["gene_id"].unique())
    print(f"[{tissue}] {len(genes)} genes to harmonize\n")

    summary = []
    for gene in genes:
        eg = e[e["gene_id"] == gene].copy()
        eg = eg.rename(columns={"ref": "eqtl_ref", "alt": "eqtl_alt",
                                "pvalue": "eqtl_p", "beta": "eqtl_beta",
                                "se": "eqtl_se", "maf": "eqtl_maf"})
        merged = gsub.merge(
            eg[["key", "eqtl_ref", "eqtl_alt", "eqtl_p",
                "eqtl_beta", "eqtl_se", "eqtl_maf", "an"]],
            on="key", how="inner")
        if merged.empty:
            continue

        # --- check 3: drop ambiguous strand SNPs near MAF 0.5 ---
        def is_amb(r):
            pair = (str(r["gwas_ref"]).upper(), str(r["gwas_alt"]).upper())
            return pair in AMBIG and r["gwas_maf"] > MAF_AMBIG_CUTOFF
        n_before = len(merged)
        merged = merged[~merged.apply(is_amb, axis=1)]
        n_amb = n_before - len(merged)

        # --- check 2: align effect allele ---
        # eQTL effect allele = eqtl_alt; GWAS effect allele = gwas_alt.
        # If they differ (ref/alt swapped), flip eQTL beta and MAF.
        def aligned_eqtl_beta(r):
            if str(r["eqtl_alt"]).upper() == str(r["gwas_alt"]).upper():
                return r["eqtl_beta"]
            else:
                return -r["eqtl_beta"]
        merged["eqtl_beta_aligned"] = merged.apply(aligned_eqtl_beta, axis=1)

        merged["eqtl_N"] = (merged["an"] / 2).round().astype("Int64")

        out = f"{outdir}/{gene}.tsv"
        cols = ["key", "pos", "rsids",
                "gwas_ref", "gwas_alt", "gwas_p", "gwas_beta", "gwas_se", "gwas_maf",
                "eqtl_ref", "eqtl_alt", "eqtl_p", "eqtl_beta_aligned",
                "eqtl_se", "eqtl_maf", "eqtl_N"]
        merged[cols].to_csv(out, sep="\t", index=False)

        min_gwas_p = merged["gwas_p"].min()
        min_eqtl_p = merged["eqtl_p"].min()
        summary.append((gene, len(merged), n_amb, min_gwas_p, min_eqtl_p))

    print("gene\t\tshared_SNPs\tambig_dropped\tmin_gwas_p\tmin_eqtl_p")
    for gene, n, na, gp, ep in summary:
        print(f"{gene}\t{n}\t{na}\t\t{gp:.2e}\t{ep:.2e}")
    print(f"\nHarmonized files in {outdir}/")
    print("A gene is a promising coloc candidate if BOTH min_gwas_p and "
          "min_eqtl_p are small (both signals present in the shared SNPs).")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(__doc__); sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
