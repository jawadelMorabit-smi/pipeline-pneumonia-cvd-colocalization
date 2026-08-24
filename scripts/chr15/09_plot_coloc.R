#!/usr/bin/env Rscript
# ============================================================
# Visualisation des resultats coloc : barres empilees PP.H0-PP.H4 par gene.
# Montre visuellement la dominance de H3 (signaux distincts) vs H4 (coloc).
#
# Usage (console RStudio) :
#   commandArgs <- function(trailingOnly=TRUE) c("results/coloc_blood.tsv","blood")
#   source("09_plot_coloc.R")
#   commandArgs <- function(trailingOnly=TRUE) c("results/coloc_lung.tsv","lung")
#   source("09_plot_coloc.R")
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
coloc_file <- args[1]
tissue     <- args[2]

# table de correspondance ID Ensembl -> symbole (depuis les figures)
symfile <- file.path(paste0("results/figures_", tissue), "gene_symbols.tsv")
sym <- NULL
if (file.exists(symfile)) {
  sym <- read.delim(symfile, stringsAsFactors = FALSE)
}

d <- read.delim(coloc_file, stringsAsFactors = FALSE)

# ajouter le symbole si dispo
label <- d$gene
if (!is.null(sym)) {
  m <- match(d$gene, sym$gene_id)
  label <- ifelse(!is.na(m) & sym$symbol[m] != d$gene, sym$symbol[m], d$gene)
}
d$label <- label

# matrice des 5 probabilites (lignes = genes, colonnes = H0..H4)
M <- t(as.matrix(d[, c("PP.H0","PP.H1","PP.H2","PP.H3","PP.H4")]))
colnames(M) <- d$label

cols <- c("grey80", "#8dd3c7", "#ffffb3", "#fb8072", "#80b1d3")  # H0..H4
# H3 = rouge (fb8072), H4 = bleu (80b1d3) : les 2 a surveiller

outpng <- file.path("results", paste0("coloc_barplot_", tissue, ".png"))
png(outpng, width = 1100, height = 650, res = 130)
par(mar = c(8, 4, 4, 8), xpd = TRUE)
barplot(M, col = cols, las = 2, ylab = "Probabilite posterieure",
        main = paste0("Colocalisation par gene — ", tissue,
                      "\n(H3 rouge = signaux distincts ; H4 bleu = colocalisation)"),
        cex.names = 0.8)
abline(h = 0.5, lty = 2, lwd = 1)   # seuil PPH4 = 0.5
legend("topright", inset = c(-0.18, 0),
       legend = c("H0 rien","H1 GWAS seul","H2 eQTL seul",
                  "H3 variants distincts","H4 meme variant"),
       fill = cols, bty = "n", cex = 0.8)
dev.off()

cat(sprintf("Figure sauvegardee : %s\n", outpng))
cat("\nGenes tries par PP.H4 :\n")
print(d[order(-d$PP.H4), c("label","PP.H3","PP.H4")], row.names = FALSE)
cat(sprintf("\nMax PP.H4 observe (%s) : %.4f  -> aucun gene ne colocalise (seuil 0.5)\n",
            tissue, max(d$PP.H4)))
