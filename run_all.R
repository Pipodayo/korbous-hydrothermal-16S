## Runs every analysis in order; outputs go to figures/ and tables/
for (f in c("Fig2_alpha_diversity.R","Fig3_PCA_heatmap.R","Fig4_RDA.R","Fig5_marine_influence.R","FigS2_UpSet.R","TableS7_robustness.R")) {
  cat("\n==>", f, "\n"); source(f, local = new.env())
}
writeLines(capture.output(sessionInfo()), "tables/sessionInfo.txt")
