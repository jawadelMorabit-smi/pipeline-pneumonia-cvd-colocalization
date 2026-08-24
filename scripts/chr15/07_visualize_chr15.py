#!/usr/bin/env python3
"""
Visualize the harmonized chr15 results for one tissue.

Produces:
  1. gene_symbols.tsv         : Ensembl ID -> gene symbol (via Ensembl REST)
  2. scatter_<tissue>.png     : per-gene min -log10(GWAS p) vs min -log10(eQTL p)
                                (top-right corner = strong double signal)
  3. locus_<tissue>_<gene>.png: for each promising gene, GWAS and eQTL
                                -log10(p) along genomic position (do the peaks
                                line up? = visual hint of colocalization)

Requires: pandas, matplotlib, requests
    pip install matplotlib requests

Usage:
    python3 07_visualize_chr15.py results/harmonized_blood blood
    python3 07_visualize_chr15.py results/harmonized_lung  lung
"""
import sys
import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")           # no display needed; saves to file
import matplotlib.pyplot as plt
import requests

ENSEMBL = "https://rest.ensembl.org"


def gene_symbol(ensg):
    """Look up a gene symbol from an Ensembl gene ID (best-effort)."""
    try:
        r = requests.get(f"{ENSEMBL}/lookup/id/{ensg}",
                         headers={"Content-Type": "application/json"},
                         timeout=20)
        if r.status_code == 200:
            return r.json().get("display_name", ensg)
    except Exception:
        pass
    return ensg


def neglog10(p):
    p = pd.to_numeric(p, errors="coerce").clip(lower=1e-300)
    return -np.log10(p)


def main(harm_dir, tissue):
    files = sorted(glob.glob(f"{harm_dir}/*.tsv"))
    if not files:
        print(f"No harmonized files in {harm_dir}")
        return
    outdir = f"results/figures_{tissue}"
    os.makedirs(outdir, exist_ok=True)

    # --- map Ensembl IDs to symbols ---
    ids = [os.path.splitext(os.path.basename(f))[0] for f in files]
    print("Looking up gene symbols...")
    sym = {i: gene_symbol(i) for i in ids}
    pd.DataFrame({"gene_id": list(sym), "symbol": list(sym.values())}
                 ).to_csv(f"{outdir}/gene_symbols.tsv", sep="\t", index=False)
    for i in ids:
        print(f"  {i} -> {sym[i]}")

    # --- collect per-gene summary for the scatter ---
    rows = []
    for f, i in zip(files, ids):
        d = pd.read_csv(f, sep="\t")
        if d.empty:
            continue
        rows.append({
            "gene": sym[i],
            "gwas": neglog10(d["gwas_p"]).max(),
            "eqtl": neglog10(d["eqtl_p"]).max(),
            "n": len(d),
        })
    summ = pd.DataFrame(rows)

    # --- (2) scatter: double-signal overview ---
    plt.figure(figsize=(7, 6))
    plt.scatter(summ["gwas"], summ["eqtl"], s=40)
    for _, r in summ.iterrows():
        plt.annotate(r["gene"], (r["gwas"], r["eqtl"]),
                     fontsize=8, xytext=(4, 4), textcoords="offset points")
    plt.axvline(-np.log10(5e-8), ls="--", lw=0.8, color="grey")
    plt.axhline(-np.log10(5e-8), ls="--", lw=0.8, color="grey")
    plt.xlabel("max -log10(GWAS p)  [pneumonia signal]")
    plt.ylabel("max -log10(eQTL p)  [expression signal]")
    plt.title(f"Double-signal overview — {tissue}\n(top-right = coloc candidates)")
    plt.tight_layout()
    plt.savefig(f"{outdir}/scatter_{tissue}.png", dpi=130)
    plt.close()
    print(f"\nsaved {outdir}/scatter_{tissue}.png")

    # --- (3) locus plots for promising genes (both signals genome-wide) ---
    GW = -np.log10(5e-8)
    promising = summ[(summ["gwas"] >= GW) & (summ["eqtl"] >= GW)]["gene"].tolist()
    print(f"Promising genes (both signals >= 5e-8): {promising}")

    for f, i in zip(files, ids):
        if sym[i] not in promising:
            continue
        d = pd.read_csv(f, sep="\t").sort_values("pos")
        fig, ax1 = plt.subplots(figsize=(9, 4.5))
        ax1.scatter(d["pos"] / 1e6, neglog10(d["gwas_p"]),
                    s=10, color="tab:blue", label="GWAS (pneumonia)")
        ax1.set_xlabel("chr15 position (Mb)")
        ax1.set_ylabel("-log10(GWAS p)", color="tab:blue")
        ax1.tick_params(axis="y", labelcolor="tab:blue")
        ax2 = ax1.twinx()
        ax2.scatter(d["pos"] / 1e6, neglog10(d["eqtl_p"]),
                    s=10, color="tab:red", alpha=0.6, label="eQTL (expression)")
        ax2.set_ylabel("-log10(eQTL p)", color="tab:red")
        ax2.tick_params(axis="y", labelcolor="tab:red")
        plt.title(f"{sym[i]} ({i}) — {tissue}\nDo the GWAS and eQTL peaks coincide?")
        fig.tight_layout()
        out = f"{outdir}/locus_{tissue}_{sym[i]}.png"
        plt.savefig(out, dpi=130)
        plt.close()
        print(f"saved {out}")

    print(f"\nAll figures in {outdir}/")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__); sys.exit(1)
    main(sys.argv[1], sys.argv[2])
