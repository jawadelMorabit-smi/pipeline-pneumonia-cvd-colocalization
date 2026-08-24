#!/usr/bin/env python3
# Phase 5: overlap of pneumonia GWAS genes with CVD gene list + Venn diagram.
# Run from Project/. Outputs: results/phase5_cvd_overlap.tsv, figures/venn_pneumonia_cvd.png
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

PNEU = {
    "chr1 (CRP)": ["CRP", "CRPP1"],
    "chr2 (IFIH1)": ["IFIH1"],
    "chr5 (CDH6)": ["CDH6"],
    "chr6 (MHC)": ["HLA-DQA1", "HLA-DQB1", "HCP5"],
    "chr9 (TNFSF15)": ["TNFSF15", "DEC1"],
    "chr12 (HNF1A)": ["HNF1A", "HNF1A-AS1"],
    "chr15": ["CHRNA3", "CHRNA5", "HYKK", "IREB2", "PSMA4"],
}
cvd = set(l.strip() for l in open("data/cvd_genes.txt") if l.strip())

rows, pneu_all = [], set()
for loc, genes in PNEU.items():
    for g in genes:
        pneu_all.add(g)
        rows.append({"locus": loc, "gene": g, "in_CVD": g in cvd})
os.makedirs("results", exist_ok=True)
pd.DataFrame(rows).sort_values(["locus", "gene"]).to_csv(
    "results/phase5_cvd_overlap.tsv", sep="\t", index=False)

inter = sorted(pneu_all & cvd)
a, b, n = len(pneu_all) - len(inter), len(cvd) - len(inter), len(inter)
print("Pneumonia: %d | CVD: %d | Overlap: %d -> %s" % (
    len(pneu_all), len(cvd), n, ", ".join(inter)))

fig, ax = plt.subplots(figsize=(8, 5.5))
ax.add_patch(Circle((-0.55, 0), 1.0, color="#d1495b", alpha=0.50))
ax.add_patch(Circle((0.55, 0), 1.0, color="#66a182", alpha=0.50))
ax.text(-1.05, 0, str(a), ha="center", va="center", fontsize=20, fontweight="bold")
ax.text(1.05, 0, "{:,}".format(b), ha="center", va="center", fontsize=18, fontweight="bold")
ax.text(0, 0, str(n), ha="center", va="center", fontsize=22, fontweight="bold")
ax.text(-0.9, 1.15, "Pneumonia locus genes\n(GWAS FinnGen)", ha="center", fontsize=11)
ax.text(0.9, 1.15, "CVD genes\n(cvd_genes.txt)", ha="center", fontsize=11)
ax.text(0, -1.55, "Shared (7): " + ", ".join(inter), ha="center",
        fontsize=11, fontweight="bold", color="#7a1f2b")
ax.set_title("Pneumonia GWAS genes vs CVD genes  -  7 / 16 shared", fontsize=13)
ax.set_xlim(-2.2, 2.2)
ax.set_ylim(-1.9, 1.6)
ax.set_aspect("equal")
ax.axis("off")
os.makedirs("figures", exist_ok=True)
plt.savefig("figures/venn_pneumonia_cvd.png", dpi=150, bbox_inches="tight")
print("OK")
