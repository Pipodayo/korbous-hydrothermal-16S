############################################################
## Korbous - Fig. S2 (UpSet, presence >= 10 reads) and Table S2 (10-read and 1-read thresholds)
## Run from the folder containing 00_common.R and data/
############################################################
source("00_common.R")

## ======================= UpSet (Fig. S2, Table S2) =======================
zs <- as.character(meta[rownames(rar), "Zone"])
ups <- function(thr, X = rar, z = zs) { pres <- sapply(zones, function(k) colSums(X[z == k, , drop = FALSE] >= thr) > 0)
  comb <- apply(pres * 1, 1, paste0, collapse = ""); d <- data.frame(combination = comb, reads = colSums(X))
  d <- d[d$combination != "0000", ]; o <- aggregate(reads ~ combination, d, function(x) c(n = length(x), r = sum(x)))
  data.frame(combination = o$combination, n_ASV = o$reads[, "n"], pct_reads = 100 * o$reads[, "r"] / sum(X)) }
u10 <- ups(10); u1 <- ups(1)
S2 <- merge(u10, u1, by = "combination", all = TRUE, suffixes = c("_10reads", "_1read")); S2 <- S2[order(-S2$n_ASV_10reads), ]
S2$compartments <- sapply(S2$combination, function(cc) paste(zones[strsplit(cc, "")[[1]] == "1"], collapse = " + "))
write.csv(S2, "tables/TableS2_UpSet.csv", row.names = FALSE)
cat("UpSet core (10 reads):", u10[u10$combination == "1111", "n_ASV"], round(u10[u10$combination == "1111", "pct_reads"], 1), "| ASVs below 10 reads everywhere:", ncol(rar) - sum(u10$n_ASV), "\n")
pres <- sapply(zones, function(k) colSums(rar[zs == k, , drop = FALSE] >= 10) > 0) * 1; pres <- pres[rowSums(pres) > 0, ]
m <- make_comb_mat(pres, mode = "distinct"); pct <- setNames(u10$pct_reads, u10$combination)[comb_name(m)]
cols <- rep("grey30", length(comb_degree(m)))
hu <- UpSet(m, set_order = zones, comb_order = order(-comb_size(m)), comb_col = cols, pt_size = unit(2.5, "mm"), lwd = 1.5,
  top_annotation = upset_top_annotation(m, gp = gpar(fill = cols, col = NA), annotation_name_rot = 90, axis_param = list(gp = gpar(fontsize = 7))),
  right_annotation = upset_right_annotation(m, gp = gpar(fill = zone_cols[zones], col = NA), axis_param = list(gp = gpar(fontsize = 7))),
  bottom_annotation = HeatmapAnnotation(`% of reads` = anno_barplot(pct, gp = gpar(fill = "grey60", col = NA), border = FALSE, height = unit(22, "mm"),
    axis_param = list(gp = gpar(fontsize = 7))), annotation_name_side = "left", annotation_name_rot = 0, annotation_name_gp = gpar(fontsize = 8)),
  row_names_gp = gpar(fontsize = 9))
cairo_pdf("figures/FigS2_UpSet.pdf", width = 174/25.4, height = 120/25.4); draw(hu, padding = unit(c(2, 22, 2, 2), "mm")); dev.off()

