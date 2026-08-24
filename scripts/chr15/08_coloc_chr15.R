#!/usr/bin/env Rscript
# ============================================================
# Phase 4 : colocalisation GWAS (pneumonie) <-> eQTL (GTEx)
# Lance coloc.abf par gene pour un tissu, et rassemble PPH0-PPH4.
#
# Entrees : les fichiers harmonises produits par 06_harmonize_chr15.py
#   results/harmonized_blood/<gene>.tsv
#   results/harmonized_lung/<gene>.tsv
#
# Usage :
#   Rscript 08_coloc_chr15.R results/harmonized_blood blood
#   Rscript 08_coloc_chr15.R results/harmonized_lung  lung
# ============================================================

# --- Parametres GWAS FinnGen (verrouilles) ---
GWAS_N <- 500186      # taille totale
GWAS_S <- 0.1687      # proportion de cas
PPH4_THRESHOLD <- 0.5 # seuil fixe par Prof. Yang

# --- lecture des arguments ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 08_coloc_chr15.R <harmonized_dir> <tissue>")
}
harm_dir <- args[1]
tissue   <- args[2]

# --- charger coloc (installer si besoin) ---
if (!requireNamespace("coloc", quietly = TRUE)) {
  install.packages("coloc", repos = "https://cloud.r-project.org")
}
suppressMessages(library(coloc))

# --- lister les fichiers harmonises (un par gene) ---
files <- list.files(harm_dir, pattern = "\\.tsv$", full.names = TRUE)
if (length(files) == 0) stop(paste("Aucun fichier .tsv dans", harm_dir))

results <- list()

for (f in files) {
  gene <- sub("\\.tsv$", "", basename(f))
  d <- tryCatch(read.delim(f, stringsAsFactors = FALSE),
                error = function(e) NULL)
  if (is.null(d) || nrow(d) < 5) {
    cat(sprintf("  [skip] %s : trop peu de SNP\n", gene)); next
  }

  # nettoyer : garder les lignes completes et varbeta > 0
  d$gwas_se <- as.numeric(d$gwas_se)
  d$eqtl_se <- as.numeric(d$eqtl_se)
  d$gwas_beta <- as.numeric(d$gwas_beta)
  d$eqtl_beta_aligned <- as.numeric(d$eqtl_beta_aligned)
  d$gwas_maf <- as.numeric(d$gwas_maf)
  d$eqtl_maf <- as.numeric(d$eqtl_maf)
  d$eqtl_N   <- as.numeric(d$eqtl_N)

  d <- d[complete.cases(d$gwas_beta, d$gwas_se, d$eqtl_beta_aligned,
                        d$eqtl_se, d$gwas_maf, d$eqtl_maf) &
         d$gwas_se > 0 & d$eqtl_se > 0 &
         d$gwas_maf > 0 & d$gwas_maf < 1 &
         d$eqtl_maf > 0 & d$eqtl_maf < 1, ]
  if (nrow(d) < 5) {
    cat(sprintf("  [skip] %s : trop peu de SNP apres nettoyage\n", gene)); next
  }

  # --- jeu GWAS (cas-temoins) ---
  ds_gwas <- list(
    snp     = d$key,
    beta    = d$gwas_beta,
    varbeta = d$gwas_se^2,
    type    = "cc",
    s       = GWAS_S,
    N       = GWAS_N
  )

  # --- jeu eQTL (quantitatif) ---
  # sdY inconnu -> coloc l'estime depuis MAF + varbeta + N (warning normal)
  ds_eqtl <- list(
    snp     = d$key,
    beta    = d$eqtl_beta_aligned,
    varbeta = d$eqtl_se^2,
    type    = "quant",
    N       = round(mean(d$eqtl_N, na.rm = TRUE)),
    MAF     = d$eqtl_maf
  )

  res <- tryCatch(
    suppressWarnings(coloc.abf(dataset1 = ds_gwas, dataset2 = ds_eqtl)),
    error = function(e) { cat(sprintf("  [err] %s : %s\n", gene, e$message)); NULL }
  )
  if (is.null(res)) next

  s <- res$summary
  results[[gene]] <- data.frame(
    gene    = gene,
    nsnps   = as.integer(s["nsnps"]),
    PP.H0   = s["PP.H0.abf"],
    PP.H1   = s["PP.H1.abf"],
    PP.H2   = s["PP.H2.abf"],
    PP.H3   = s["PP.H3.abf"],
    PP.H4   = s["PP.H4.abf"],
    row.names = NULL
  )
  cat(sprintf("  [ok] %s : nsnps=%d  PP.H4=%.3f\n",
              gene, as.integer(s["nsnps"]), s["PP.H4.abf"]))
}

if (length(results) == 0) stop("Aucun resultat coloc produit.")

out <- do.call(rbind, results)
out <- out[order(-out$PP.H4), ]

# --- sauvegarde ---
outfile <- file.path("results", paste0("coloc_", tissue, ".tsv"))
write.table(out, outfile, sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n================ RESULTATS COLOC (", tissue, ") ================\n")
print(out, row.names = FALSE, digits = 3)
cat(sprintf("\nGenes colocalises (PP.H4 > %.2f) : %s\n",
            PPH4_THRESHOLD,
            paste(out$gene[out$PP.H4 > PPH4_THRESHOLD], collapse = ", ")))
cat(sprintf("Sauvegarde : %s\n", outfile))
