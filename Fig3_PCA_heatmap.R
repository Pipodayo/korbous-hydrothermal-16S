############################################################
## Korbous - Fig. 3 (PCA of CLR data; heatmap of the 30 most abundant ASVs)
## Run from the folder containing 00_common.R and data/
############################################################
source("00_common.R")

## ======================= PCA + HEATMAP (Fig. 3) =======================
pc <- prcomp(Y, center = TRUE, scale. = FALSE); ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
sc <- data.frame(pc$x[, 1:2], Sample = rownames(Y), Compartment = meta[rownames(Y), "Zone"])
write.csv(sc[, c("Sample","PC1","PC2")], "tables/Fig3A_PCA_scores.csv", row.names = FALSE)
ld <- data.frame(ASV = colnames(Y), PC1 = pc$rotation[, 1], PC2 = pc$rotation[, 2])
top <- ld[order(abs(ld$PC1), decreasing = TRUE)[1:2], ]
labels30 <- c(ASV1="Pseudoalteromonas sp.",ASV2="Sulfurimonas sp.",ASV3="Colwellia sp.",ASV4="Sphingomonas sp.",ASV5="Thiofaba sp.",
  ASV6="Colwellia sp.",ASV7="Aggregatilineales-affiliated",ASV8="Sulfurimonas sp.",ASV9="Pseudomonas sp.",ASV10="Chlorobaculum sp.",
  ASV11="Thermus sp.",ASV12="Hydrogenothermus sp.",ASV13="Vibrio sp.",ASV14="Anaerolineaceae-affiliated",ASV15="Sulfurivirga sp.",
  ASV16="Ca. Bipolaricaulis sp.",ASV17="Desulfacinum sp.",ASV18="Sphingomonas sp.",ASV19="Sulfurimonas sp.",ASV20="unclassified Bacteria",
  ASV21="Leifsonia sp.",ASV22="Anaerolineaceae-affiliated",ASV23="Colwelliaceae-affiliated",ASV24="Anaerolineaceae-affiliated",
  ASV26="Sulfurimonas sp.",ASV27="Thermotoga sp.",ASV28="Ca. Thiobios-related",ASV31="Paracoccaceae-affiliated",ASV32="Bathyarchaeia-affiliated",ASV33="SJA-15 clade",ASV36="Thermus sp.",ASV41="Sulfurimonas sp.")
## arrow labels as plotmath: genus names in italics, higher ranks in plain text
top$Label <- sapply(top$ASV, function(a) { t <- labels30[[a]]
  if (grepl(" sp\\.$", t)) sprintf('%s~italic(%s)', a, sub(" sp\\.$", "", t)) else sprintf('%s~"%s"', a, t) })
mult <- 0.4 * min(max(abs(sc$PC1)) / max(abs(top$PC1)), max(abs(sc$PC2)) / max(abs(top$PC2)))
fig3A <- ggplot(sc, aes(PC1, PC2)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey70") + geom_vline(xintercept = 0, linetype = 2, colour = "grey70") +
  geom_segment(data = top, aes(x = 0, y = 0, xend = PC1 * mult, yend = PC2 * mult), inherit.aes = FALSE, colour = "grey25",
               arrow = arrow(length = unit(2, "mm"), type = "closed"), linewidth = 0.45) +
  geom_text(data = top, aes(PC1 * mult, PC2 * mult, label = Label), inherit.aes = FALSE, size = 2.2, parse = TRUE,
            hjust = ifelse(rank(-top$PC2) == 1, 0.48, 0.35), vjust = ifelse(rank(-top$PC2) == 1, 0, 1), nudge_y = ifelse(rank(-top$PC2) == 1, 3.2, -2.5)) +
  geom_point(aes(colour = Compartment, shape = Compartment), size = 2.8) +
  geom_text_repel(aes(label = Sample), size = 2.1, seed = 123, segment.colour = "grey65") +
  scale_colour_manual(values = zone_cols, name = NULL) + scale_shape_manual(values = zone_shapes, name = NULL) +
  labs(x = sprintf("PC1 (%.2f%%)", ve[1]), y = sprintf("PC2 (%.2f%%)", ve[2])) + theme_k() +
  theme(legend.position = "bottom", legend.text = element_text(size = 7)) + guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))

top30 <- names(sort(colMeans(rel), decreasing = TRUE))[1:30]
stopifnot(all(top30 %in% names(labels30)))
mat <- t(rel[, top30]); mat_z <- t(scale(t(mat))); mat_z[is.na(mat_z)] <- 0
lab_expr <- as.expression(lapply(top30, function(a) { t <- labels30[[a]]
  if (grepl("^Ca\\. ", t)) return(bquote(plain(.(paste0(a, " "))) * italic("Ca.") * plain(.(sub("^Ca\\.", "", t)))))
  if (grepl(" sp\\.$", t)) return(bquote(plain(.(paste0(a, " "))) * italic(.(sub(" sp\\.$", "", t))) * plain(" sp.")))
  if (grepl("-affiliated$", t) && !grepl("aceae|ales|ia-|clade", t)) return(bquote(plain(.(paste0(a, " "))) * italic(.(sub("-affiliated$", "", t))) * plain("-affiliated")))
  bquote(plain(.(paste(a, t)))) }))
zv <- meta[sample_order, "Zone"]
ht <- Heatmap(mat_z, name = "Z-score", top_annotation = HeatmapAnnotation(Compartment = zv, col = list(Compartment = zone_cols), show_legend = FALSE, show_annotation_name = FALSE),
  right_annotation = rowAnnotation(`Rel. ab. (%)` = anno_barplot(rowMeans(mat), border = FALSE, gp = gpar(fill = "grey70", col = NA), width = unit(12, "mm"),
                                     axis_param = list(gp = gpar(fontsize = 6))), annotation_name_gp = gpar(fontsize = 7), annotation_name_rot = 0, annotation_name_side = "top"),
  cluster_rows = TRUE, clustering_distance_rows = "euclidean", clustering_method_rows = "complete", cluster_columns = FALSE,
  column_split = zv, column_title = NULL, column_names_rot = 45, column_names_gp = gpar(fontsize = 7),
  row_labels = lab_expr, row_names_gp = gpar(fontsize = 7), row_names_max_width = unit(55, "mm"),
  col = colorRamp2(c(-2, 0, 2), c("#4575b4", "white", "#d73027")),
  heatmap_legend_param = list(title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize = 6)))
cairo_pdf("figures/Fig3_PCA_heatmap.pdf", width = 174/25.4, height = 225/25.4)
grid.newpage(); pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(92, 133), "mm"))))
print(fig3A + labs(tag = "A") + theme(legend.position = "right") + guides(colour = guide_legend(ncol = 1), shape = guide_legend(ncol = 1)),
      vp = viewport(layout.pos.row = 1, layout.pos.col = 1, width = unit(150, "mm")))
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(ht, heatmap_legend_side = "right", newpage = FALSE)
grid.text("B", x = unit(2, "mm"), y = unit(1, "npc") - unit(3, "mm"), gp = gpar(fontface = "plain", fontsize = 10)); popViewport(2); dev.off()

