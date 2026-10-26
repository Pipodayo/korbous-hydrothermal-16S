############################################################
## Korbous - Fig. 5 (A: ASVs dominating the hydrothermal waters; B: ion ratios; C: marine planktonic lineages) and Tables S8-S9
## Run from the folder containing 00_common.R and data/
############################################################
source("00_common.R")

## ======================= WATER ASVs AND MARINE INFLUENCE (Fig. 5, Tables S8-S9) =======================
ions <- c(Sodium = "Na^'+'", Potassium = "K^'+'", Magnesium = "Mg^'2+'", Calcium = "Ca^'2+'")
A <- expand.grid(Sample = c("TWS","MAWS","MBWS"), Ion = names(ions), stringsAsFactors = FALSE)
A$ratio <- mapply(function(s, i) (meta[s, i] / meta[s, "Chloride"]) / (meta_all["MWC", i] / meta_all["MWC", "Chloride"]), A$Sample, A$Ion)
write.csv(A, "tables/Fig5B_ion_ratios.csv", row.names = FALSE)
A$Ion <- factor(A$Ion, levels = names(ions)); A$Sample <- factor(A$Sample, levels = c("TWS","MAWS","MBWS"))
labt <- paste(tax$Order, tax$Family, tax$Genus); names(labt) <- rownames(tax)
grp <- list("SAR11 (Pelagibacterales)" = grepl("Pelagibacterales", labt),
            "SAR86, SAR116, AEGEAN-169, OM43" = grepl("SAR86|SAR116|Puniceispirillales|AEGEAN-169|OM43", labt),
            "NS marine groups (Flavobacteriales)" = grepl("NS2b|NS3|NS4|NS5|NS7|NS9", labt),
            "Synechococcus/Prochlorococcus" = grepl("Synechococcus|Prochlorococcus", labt))
gm <- sapply(grp, identity); rownames(gm) <- names(labt); stopifnot(all(rowSums(gm) <= 1))   # no ASV in two groups
B <- do.call(rbind, lapply(names(grp), function(g) { a <- intersect(names(labt)[grp[[g]]], colnames(rel)); data.frame(Sample = sample_order, Group = g, pct = rowSums(rel[sample_order, a, drop = FALSE])) }))
B$Sample <- factor(B$Sample, levels = sample_order); B$Group <- factor(B$Group, levels = names(grp))
short <- c("Marine reference" = "Marine\nref.", "Hydrothermal water" = "Hydroth.\nwaters", "Submarine hydrothermal sediment" = "Submarine\nhydrothermal\nsediments", "Terrestrial hydrothermal sediment" = "Terrestrial\nhydrothermal\nsediments")
B$Comp <- factor(short[as.character(meta[as.character(B$Sample), "Zone"])], levels = short)
tot <- aggregate(pct ~ Sample + Comp, B, sum)
write.csv(reshape(B[, c("Sample","Group","pct")], idvar = "Sample", timevar = "Group", direction = "wide"), "tables/Fig5C_marine_lineages.csv", row.names = FALSE)
print(tot)
mi_asv <- names(labt)[rowSums(gm) == 1]; mi_asv <- mi_asv[colSums(rar[, mi_asv, drop = FALSE]) > 0]
S9 <- data.frame(ASV = mi_asv, Group = apply(gm[mi_asv, , drop = FALSE], 1, function(r) names(grp)[which(r)]),
                 tax[mi_asv, c("Phylum","Class","Order","Family","Genus")], round(t(rel[, mi_asv, drop = FALSE]), 4), check.names = FALSE)
S9 <- S9[order(S9$Group, -S9$TWS), ]
write.csv(S9, "tables/TableS9_marine_indicator_ASVs.csv", row.names = FALSE)

## ---- Table S8: Bray-Curtis similarity (% of community shared) between samples, rarefied table ----
bc_sim <- round(100 * (1 - as.matrix(vegdist(rel[sample_order, ] / 100, method = "bray"))), 1)
write.csv(bc_sim, "tables/TableS8_BrayCurtis_similarity.csv")
print(bc_sim[c("TWS","MAWS","MBWS"), ])
ink <- "#0b0b0b"; ink2 <- "#52514e"
th <- theme_minimal(base_size = 8) + theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
  panel.grid.major.y = element_line(colour = "#e6e5e0", linewidth = 0.3), axis.text = element_text(colour = ink2), strip.text = element_text(size = 7))
pA <- ggplot(A, aes(Ion, ratio, shape = Sample)) + geom_hline(yintercept = 1, colour = ink2, linetype = "22") +
  geom_point(position = position_dodge(width = 0.55), size = 2.4, colour = ink, fill = "white", stroke = 0.6) +
  scale_shape_manual(values = c(TWS = 24, MAWS = 22, MBWS = 23), name = NULL) + scale_y_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3), limits = c(0.5, 3.6)) +
  scale_x_discrete(labels = parse(text = ions)) + labs(x = NULL, y = "(X/Cl) sample / (X/Cl) seawater") + th + theme(legend.position = c(0.22, 0.9))
pB <- ggplot(B, aes(Sample, pct, fill = Group)) + geom_col(width = 0.72, colour = "#fcfcfb", linewidth = 0.25) +
  geom_text(data = tot[tot$pct >= 1, ], aes(Sample, pct, label = sprintf("%.1f", pct)), inherit.aes = FALSE, vjust = -0.4, size = 2.4) +
  facet_grid(. ~ Comp, scales = "free_x", space = "free_x") + scale_fill_manual(values = c("#1f4e79","#8c6bb1","#8a8a8a","#c9a227"), name = NULL,
    labels = c("SAR11 (Pelagibacterales)", "SAR86, SAR116, AEGEAN-169, OM43", "NS marine groups (Flavobacteriales)", expression(italic(Synechococcus)*"/"*italic(Prochlorococcus)))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) + labs(x = NULL, y = "Marine planktonic lineages\n(% of reads)") +
  th + theme(legend.position = "top", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.key.size = unit(3, "mm")) + guides(fill = guide_legend(nrow = 2))
asv4 <- c("ASV2","ASV8","ASV26","ASV6")
C <- do.call(rbind, lapply(asv4, function(a) data.frame(Sample = sample_order, ASV = paste(a, tax[a, "Genus"]), pct = rel[sample_order, a])))
C$Sample <- factor(C$Sample, levels = sample_order); C$ASV <- factor(C$ASV, levels = rev(unique(C$ASV)))
C$Comp <- factor(short[as.character(meta[as.character(C$Sample), "Zone"])], levels = short)
pC <- ggplot(C, aes(Sample, ASV, fill = pmax(pct, 0.01))) + geom_tile(colour = "#fcfcfb", linewidth = 0.6) +
  geom_text(data = C[C$pct >= 1, ], aes(label = sprintf("%.0f", pct), colour = pct > 15), size = 2.3, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = ink)) + facet_grid(. ~ Comp, scales = "free_x", space = "free_x") +
  scale_fill_gradient(low = "#efedfb", high = "#4a3aa7", trans = "log10", breaks = c(0.01, 0.1, 1, 10), labels = c("<0.01","0.1","1","10"), name = "% of reads") +
  labs(x = NULL, y = NULL) + th + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_y_discrete(labels = function(x) parse(text = sub("^(ASV[0-9]+) (.*)$", "plain('\\1 ')*italic('\\2')", x)))
ggsave("figures/Fig5_water_ASVs_marine_influence.pdf", (pC / ((pA | pB) + plot_layout(widths = c(0.4, 1)))) + plot_layout(heights = c(0.5, 1)) + plot_annotation(tag_levels = "A"),
       width = 174, height = 160, units = "mm", device = cairo_pdf)

