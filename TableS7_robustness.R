############################################################
## Korbous - Table S7 (robustness to the inclusion of MWC and to the removal of contaminant-affiliated genera) and Table S7b
## Run from the folder containing 00_common.R and data/
############################################################
source("00_common.R")

## ======================= Robustness (Table S7) =======================
cand <- c("Temperature","Chloride","pH","Calcium","Magnesium","Potassium","Sulfate","Sulfur")
select_rda <- function(Yx, samples) {           # same VIF screening and forward selection as Fig4_RDA.R
  e <- as.data.frame(scale(meta_all[samples, cand])); vars <- cand
  repeat { v <- vif.cca(rda(as.formula(paste("Yx ~", paste(vars, collapse = "+"))), data = e))
    if (max(v) <= 10) break; vars <- setdiff(vars, names(which.max(v))) }
  set.seed(123); sel <- ordiR2step(rda(Yx ~ 1, data = e), scope = formula(rda(as.formula(paste("Yx ~", paste(vars, collapse = "+"))), data = e)),
                                   Pin = 0.05, R2scope = TRUE, permutations = perm, trace = FALSE)
  set.seed(123); p <- anova(sel, permutations = perm)$`Pr(>F)`[1]
  list(terms = paste(attr(terms(sel), "term.labels"), collapse = " + "), r2 = RsquareAdj(sel)$adj.r.squared, p = p)
}
run <- function(label, X, drop_cont = FALSE) {
  if (drop_cont) X <- X[, !colnames(X) %in% cont_asv]
  X <- X[, colSums(X) > 0]
  Yx <- clr_matrix(X); p <- prcomp(Yx); v <- 100 * p$sdev^2 / sum(p$sdev^2)
  common <- intersect(rownames(Y), rownames(Yx))      # Procrustes on the samples shared with the main analysis
  set.seed(123); pr <- protest(prcomp(Y)$x[common, 1:2], p$x[common, 1:2], permutations = 9999)
  s <- select_rda(Yx, rownames(Yx))
  data.frame(Scenario = label, n_samples = nrow(Yx), n_ASV_CLR = ncol(Yx), PC1_pct = round(v[1], 2), PC2_pct = round(v[2], 2),
             Procrustes_r_vs_main = round(pr$t0, 3), Procrustes_p = pr$signif,
             RDA_selected = s$terms, R2adj_pct = round(100 * s$r2, 1), p_model = s$p)
}
S7 <- rbind(run("Main analysis (14 samples, MWC excluded)", raw),
            run("MWC included (15 samples)", raw_all),
            run("Contaminant genera removed", raw, drop_cont = TRUE))
write.csv(S7, "tables/TableS7_robustness.csv", row.names = FALSE); print(S7)

## Table S7b: reads assigned to contaminant genera, all 15 libraries (MWC included for reference)
rel_all <- 100 * rar_all / rowSums(rar_all)
cont_frac <- data.frame(Sample = sample_all, pct_reads_contaminant_genera = round(rowSums(rel_all[, intersect(cont_asv, colnames(rel_all))]), 2),
                        used_in_community_analyses = sample_all != "MWC")
write.csv(cont_frac, "tables/TableS7b_contaminant_fraction.csv", row.names = FALSE)
