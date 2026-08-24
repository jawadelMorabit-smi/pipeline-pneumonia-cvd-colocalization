# Pneumonia–CVD GWAS–eQTL Colocalization — Results Synthesis

**Team:** Jaouad El Morabit (Data & QC), Cheikh Taha (Colocalization), Loua Celistin (Visualization), Hamza Arbaoui (Interpretation & Writing) · **Supervisor:** Prof. Cherng Yang
*Updated 16 July 2026 — includes the liver follow-up (D3). This is the up-to-date base for the slide deck.*

## Objective
Test whether pneumonia-susceptibility loci and cardiovascular disease (CVD) share **causal variants** that act through tissue-specific gene expression (cis-eQTL), using statistical **colocalization**.

Logic chain: pneumonia GWAS signal → colocalizes with a tissue eQTL → implicates a candidate gene → check if that gene is in the CVD gene list → biological interpretation.

## Supervisor decisions (validated 16/07)
- **D1** — Primary GWAS = FinnGen DF13 (alternative GCST90624053 is not independent: contains FinnGen R10).
- **D2** — No coloc.susie; single-causal-variant assumption presented as a limitation, susie as future work.
- **D3** — Liver added for CRP and HNF1A (see below).

## Data & methods
- **GWAS:** FinnGen Data Freeze 13, endpoint `J10_PNEUMONIA` (European, GRCh38/hg38). 21.3M variants; **501 genome-wide significant** (p < 5e-8). Coloc parameters: N = 500,186; cases = 84,364; s = 0.1687; type = "cc".
- **eQTL:** GTEx v8 via EBI eQTL Catalogue API v2 — whole blood (QTD000356, n = 670), lung (QTD000271, n = 510), and **liver (QTD000266, n = 208)** for the CRP/HNF1A loci.
- **Harmonization (3 checks):** SNP match by position + sorted allele pair; effect-allele alignment (beta sign flipped when ref/alt swapped); removal of ambiguous A/T & C/G SNPs near MAF 0.5.
- **Colocalization:** `coloc.abf` (beta/varbeta parameterization), per gene per tissue. Significance threshold **PP.H4 > 0.5**. Build-matched (both hg38, no liftOver).

## Loci tested (colocalization)
| Locus | GWAS hits | Key candidate genes (strong cis-eQTL) | Best PP.H4 | Colocalization |
|---|---|---|---|---|
| chr15 (15q25) | 140 | CTSH, IREB2, CIB2, TBC1D2B (blood); RASGRF1 (lung) | ~0.02 | **No** (PP.H3 ≈ 1) |
| chr1 (CRP) | 82 | SLAMF8, IGSF8, COPA, AIM2, MNDA | 0.10 | **No** (PP.H3 ≈ 1) |
| chr12 (HNF1A) | 12 | P2RX4, COQ5, C12orf43 | 0.29 | **No** (PP.H3 ≈ 1) |
| chr9 (TNFSF15) | 11 | TNFSF15, ORM1 | <0.01 | **No** (PP.H3 ≈ 1) |

**Master result:** across all four loci and **74 gene–tissue tests (blood, lung, liver)**, **no colocalization** (max PP.H4 = 0.29, LINC01089 in lung). Genes with strong GWAS **and** strong eQTL signals consistently give **PP.H3 ≈ 1** — both signals are real but driven by **distinct causal variants**. Full table: `results/MASTER_coloc_results.tsv`.

## Liver follow-up (D3)
Liver eQTL (GTEx v8) tested for the CRP and HNF1A loci, on the supervisor's suggestion (both genes are primarily liver-expressed).
- **CRP and HNF1A still have no cis-eQTL in liver** → the lead genes remain untestable even in their primary tissue.
- Where a neighbouring gene was testable, the result is again **PP.H3 ≈ 1** (no colocalization).
- **The result holds across three tissues (blood, lung, liver)** — the limitation "lead genes lack cis-eQTL" is now demonstrated empirically, not just asserted.

## Loci not testable
- **chr2 (IFIH1):** a single genome-wide significant variant, and GTEx has **no cis-eQTL for IFIH1** in blood, lung or liver → colocalization not possible.
- **chr5 (CDH6):** too few hits (underpowered).
- **chr6 (MHC):** excluded by design (extreme, extended LD makes coloc uninterpretable).

## Gene-level overlap with CVD (Phase 5)
The CVD list (`cvd_genes.txt`) has **4,311 lines but 2 duplicated symbols (CDKN2A, TNS1) → 4,309 unique genes** (~20% of protein-coding genes; a hit is suggestive, not proof).
- **View A — GWAS-locus genes ∩ CVD: 7/16** → CHRNA3, CHRNA5, CRP, HCP5, HLA-DQB1, HNF1A, PSMA4. (`figures/venn_pneumonia_cvd.png`)
- **View B — eQTL candidate genes (double signal) ∩ CVD: 7/26** → CASQ1, COQ5, CTSH, IGSF8, ORM1, P2RX4, PIGM. (`figures/venn_candidates_cvd.png`)

## Conclusion
Pneumonia and CVD **share genomic loci and genes** (notably CRP and HNF1A, direct inflammation↔CVD links), but colocalization finds **no evidence of shared causal variants** acting through cis-regulation of gene expression in whole blood, lung or liver, under the single-causal-variant assumption of `coloc.abf`. The shared genetic architecture appears to be at the level of overlapping regions/genes rather than identical causal variants. This is a consistent, reproducible result across four loci and three tissues — and the pipeline is validated (it recovers the known 15q25 genes HYKK and CHRNA3/5).

## Limitations & next steps
- **Single-causal-variant assumption:** `coloc.abf` assumes one causal variant per region. Several loci (chr15, CRP, HNF1A) are multi-signal, where H3 can mask a hidden colocalization → **`coloc.susie`** (requires an LD reference panel, e.g. 1000 Genomes EUR) is the rigorous follow-up.
- **Three tissues** (blood, lung, liver); further disease-relevant tissues not tested (e.g. artery/endothelial for CVD, stomach).
- GTEx lacks cis-eQTLs for the lead genes themselves (CRP, IFIH1, HNF1A), so their direct expression-mediated effect cannot be assessed here.
- Report PP.H4 at both 0.5 and 0.8 (sensitivity); cross-check whole-blood hits against eQTLGen; extend the validated pipeline to a second infection trait.
