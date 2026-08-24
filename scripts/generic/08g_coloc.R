#!/usr/bin/env Rscript
# ============================================================
# Coloc generique : coloc.abf par gene pour un locus + tissu.
# Lancer depuis la racine Project/. Sortie : results/<locus>/coloc_<tissu>.tsv
# Usage (dans RStudio, via source apres avoir redefini commandArgs) :
#   setwd("D:/000my_Carrier/BIAM_MASTER/S2/05_Bioinformatics_1/Project")
#   commandArgs <- function(trailingOnly=TRUE) c("results/CRP/harmonized_blood","CRP","blood")
#   source("scripts/generic/08g_coloc.R")
# ============================================================

GWAS_N <- 500186
GWAS_S <- 0.1687
PPH4_THRESHOLD <- 0.5

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: <harmonized_dir> <locus> <tissue>")
harm_dir <- args[1]; locus <- args[2]; tissue <- args[3]

if (!requireNamespace("coloc", quietly = TRUE))
  install.packages("coloc", repos = "https://cloud.r-project.org")
suppressMessages(library(coloc))

files <- list.files(harm_dir, pattern = "\\.tsv$", full.names = TRUE)
if (length(files) == 0) stop(paste("Aucun .tsv dans", harm_dir))

results <- list()
for (f in files) {
  gene <- sub("\\.tsv$", "", basename(f))
  d <- tryCatch(read.delim(f, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(d) || nrow(d) < 5) next
  for (col in c("gwas_beta","gwas_se","eqtl_beta_aligned","eqtl_se",
                "gwas_maf","eqtl_maf","eqtl_N"))
    d[[col]] <- as.numeric(d[[col]])
  d <- d[complete.cases(d$gwas_beta,d$gwas_se,d$eqtl_beta_aligned,
                        d$eqtl_se,d$gwas_maf,d$eqtl_maf) &
         d$gwas_se>0 & d$eqtl_se>0 &
         d$gwas_maf>0 & d$gwas_maf<1 & d$eqtl_maf>0 & d$eqtl_maf<1, ]
  if (nrow(d) < 5) next

  ds_gwas <- list(snp=d$key, beta=d$gwas_beta, varbeta=d$gwas_se^2,
                  type="cc", s=GWAS_S, N=GWAS_N)
  ds_eqtl <- list(snp=d$key, beta=d$eqtl_beta_aligned, varbeta=d$eqtl_se^2,
                  type="quant", N=round(mean(d$eqtl_N, na.rm=TRUE)), MAF=d$eqtl_maf)

  res <- tryCatch(suppressWarnings(coloc.abf(ds_gwas, ds_eqtl)),
                  error=function(e){cat(sprintf("  [err] %s: %s\n",gene,e$message));NULL})
  if (is.null(res)) next
  s <- res$summary
  results[[gene]] <- data.frame(gene=gene, nsnps=as.integer(s["nsnps"]),
    PP.H0=s["PP.H0.abf"], PP.H1=s["PP.H1.abf"], PP.H2=s["PP.H2.abf"],
    PP.H3=s["PP.H3.abf"], PP.H4=s["PP.H4.abf"], row.names=NULL)
  cat(sprintf("  [ok] %s : nsnps=%d  PP.H4=%.3f\n",
              gene, as.integer(s["nsnps"]), s["PP.H4.abf"]))
}

if (length(results)==0) stop("Aucun resultat coloc.")
out <- do.call(rbind, results); out <- out[order(-out$PP.H4), ]
dir.create(file.path("results", locus), showWarnings = FALSE, recursive = TRUE)
outfile <- file.path("results", locus, paste0("coloc_", tissue, ".tsv"))
write.table(out, outfile, sep="\t", quote=FALSE, row.names=FALSE)

cat(sprintf("\n======= COLOC %s (%s) =======\n", locus, tissue))
print(out, row.names=FALSE, digits=3)
cat(sprintf("\nColoc (PP.H4 > %.2f) : %s\n", PPH4_THRESHOLD,
    paste(out$gene[out$PP.H4 > PPH4_THRESHOLD], collapse=", ")))
cat(sprintf("Sauvegarde : %s\n", outfile))
