############################################################
## Korbous (Aïn Kanassira) - common objects sourced by every script
## Input: exported, taxonomy,
## validated chemistry (data/), PD values from picante (Table S1)
## Output: figures/ and tables/ used in the manuscript
## R 4.3, vegan, zCompositions, compositions, ComplexHeatmap, ggplot2
############################################################
suppressPackageStartupMessages({
  library(vegan); library(zCompositions); library(compositions)
  library(ggplot2); library(patchwork); library(ggrepel); library(dplyr); library(tidyr)
  library(ComplexHeatmap); library(circlize); library(grid); library(scales)
})
set.seed(123)
dir.create("figures", showWarnings = FALSE); dir.create("tables", showWarnings = FALSE)

## ---------------- data ----------------
## MWC is the chemical seawater reference only: its library (failed direct amplification) is dominated by
## reagent-associated genera and contains no SAR11, so it is excluded from all community analyses.
## raw_all / rar_all / meta_all keep the 15 samples (MWC chemistry, MWC contaminant fraction).
sample_all   <- c("MWC","MSed","TWS","MAWS","MBWS","MASed","MA10","MA40","MBSed","MB10","MB40","TSed1","TSed2","T10","T40")
sample_order <- setdiff(sample_all, "MWC")
zones <- c("Marine reference","Hydrothermal water","Submarine hydrothermal sediment","Terrestrial hydrothermal sediment")
zone_cols <- setNames(c("#66c2a5","#fc8d62","#8da0cb","#e78ac3"), zones)
zone_shapes <- setNames(c(16,17,15,18), zones)
raw_all <- t(as.matrix(read.csv("data/asv_counts_raw_corrected.csv", row.names = 1, check.names = FALSE)))[sample_all, ]
rar_all <- t(as.matrix(read.csv("data/asv_counts_rarefied_corrected.csv", row.names = 1, check.names = FALSE)))[sample_all, ]
raw <- raw_all[sample_order, ]; raw <- raw[, colSums(raw) > 0]
rar <- rar_all[sample_order, ]; rar <- rar[, colSums(rar) > 0]
tax <- read.csv("data/taxonomy.csv", row.names = 1)
meta_all <- read.csv("data/metadata.csv", row.names = 1)[sample_all, ]
meta <- meta_all[sample_order, ]
meta$Zone <- factor(meta$Zone, levels = zones)
rel <- 100 * rar / rowSums(rar)
theme_k <- function(b = 8) theme_classic(base_size = b) + theme(axis.text = element_text(colour = "black"), legend.key = element_blank())
stopifnot(length(unique(rowSums(rar))) == 1)
cat("Rarefaction depth:", unique(rowSums(rar)), "| samples:", nrow(raw), "| ASVs:", ncol(raw), "| mean depth:", round(mean(rowSums(raw))), "+/-", round(sd(rowSums(raw))), "\n")

## ---------------- contaminant genera (Salter 2014; Eisenhofer 2019) ----------------
cont_genera <- c("Sphingomonas","Pseudomonas","Leifsonia","Bradyrhizobium","Caulobacter","Candidatus Obscuribacter",
  "Ralstonia","Burkholderia-Caballeronia-Paraburkholderia","Methylobacterium-Methylorubrum","Cutibacterium","Staphylococcus",
  "Streptococcus","Lactobacillus","Acinetobacter","Delftia","Stenotrophomonas","Corynebacterium","Pantoea","Bacillus",
  "Escherichia-Shigella","Novosphingobium","Mesorhizobium","Afipia","Curvibacter","Massilia","Microbacterium","Brevundimonas")
cont_asv <- rownames(tax)[tax$Genus %in% cont_genera]

## ---------------- CLR helper (same procedure for all analyses) ----------------
clr_matrix <- function(X, min_samples = 3) {
  X <- X[, colSums(X > 0) >= min_samples, drop = FALSE]
  P <- as.matrix(cmultRepl(X, method = "CZM", label = 0, output = "prop", z.warning = 1, z.delete = FALSE))
  Y <- as.matrix(clr(acomp(P))); Y[, apply(Y, 2, sd) > 0, drop = FALSE]
}
Y <- clr_matrix(raw)
cat("CLR matrix:", dim(Y), "\n")

perm <- how(nperm = 9999)   # permutation design used by all tests
