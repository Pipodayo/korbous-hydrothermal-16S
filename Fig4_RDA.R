############################################################
## Korbous - Fig. 4 (RDA, predictors retained by forward selection) and Tables S3-S6 (Spearman, VIF, forward selection, permutation tests)
## Run from the folder containing 00_common.R and data/
############################################################
source("00_common.R")

## ======================= RDA (Fig. 4, Tables S3-S6) =======================
cand <- c("Temperature","Chloride","pH","Calcium","Magnesium","Potassium","Sulfate","Sulfur")
env <- meta[rownames(Y), cand]; env_z <- as.data.frame(scale(env))
write.csv(round(cor(env, method = "spearman"), 2), "tables/TableS3_Spearman_env.csv")
vars <- cand; vlog <- list(); st <- 1
repeat { v <- vif.cca(rda(as.formula(paste("Y ~", paste(vars, collapse = "+"))), data = env_z))
  vlog[[st]] <- data.frame(Step = st, Variable = names(v), VIF = round(v, 1), Removed = names(v) == names(which.max(v)) & max(v) > 10)
  if (max(v) <= 10) break; vars <- setdiff(vars, names(which.max(v))); st <- st + 1 }
write.csv(do.call(rbind, vlog), "tables/TableS4_VIF_screening.csv", row.names = FALSE)
perm <- how(nperm = 9999)
rda_glob <- rda(as.formula(paste("Y ~", paste(vars, collapse = "+"))), data = env_z); rda0 <- rda(Y ~ 1, data = env_z)
set.seed(123); sel <- ordiR2step(rda0, scope = formula(rda_glob), Pin = 0.05, R2scope = TRUE, permutations = perm, trace = FALSE)
S4 <- data.frame(Step = rownames(sel$anova), as.data.frame(sel$anova), check.names = FALSE)
## global test of the model with all predictors retained after VIF screening (Blanchet et al. 2008: selection only if significant)
set.seed(123); tfull <- anova(rda_glob, permutations = perm)
S4[S4$Step == "<All variables>", c("Df", "F", "Pr(>F)")] <- c(tfull$Df[1], tfull$F[1], tfull$`Pr(>F)`[1])
cat(sprintf("Global test of the candidate model: adjR2 %.4f, F %.3f, p %.4f\n", RsquareAdj(rda_glob)$adj.r.squared, tfull$F[1], tfull$`Pr(>F)`[1]))
write.csv(S4, "tables/TableS5_forward_selection.csv", row.names = FALSE)
set.seed(123); tg <- anova(sel, permutations = perm); set.seed(123); tm <- anova(sel, by = "margin", permutations = perm); set.seed(123); ta <- anova(sel, by = "axis", permutations = perm)
S5 <- rbind(data.frame(Test = "Global", Term = rownames(tg), as.data.frame(tg), check.names = FALSE),
            data.frame(Test = "Marginal", Term = rownames(tm), as.data.frame(tm), check.names = FALSE),
            data.frame(Test = "Axis", Term = rownames(ta), as.data.frame(ta), check.names = FALSE))
write.csv(S5, "tables/TableS6_RDA_tests.csv", row.names = FALSE)
r2sel <- RsquareAdj(sel)$adj.r.squared
steps <- S4[grepl("^\\+", S4$Step), ]; r2cum <- steps$R2.adj; terms_sel <- sub("^\\+ ", "", steps$Step)
r2s <- r2cum[1]
axp <- 100 * sel$CCA$eig / sel$tot.chi
cat(sprintf("RDA: %s | adjR2 %.4f (first term %.4f); axes %.2f %.2f; global p %.4f\n", paste(terms_sel, collapse = " + "), r2sel, r2s, axp[1], axp[2], tg$`Pr(>F)`[1]))

site <- as.data.frame(scores(sel, display = "sites", choices = 1:2, scaling = 2)); site$Sample <- rownames(site); site$Compartment <- meta[site$Sample, "Zone"]
bp <- as.data.frame(scores(sel, display = "bp", choices = 1:2, scaling = 2)); vlab <- c(Sulfur = "Dissolved sulfide", pH = "pH", Temperature = "Temperature", Chloride = "Chloride", Calcium = "Calcium",
          Magnesium = "Magnesium", Potassium = "Potassium", Sulfate = "Sulfate")
bp$Label <- vlab[rownames(bp)]
m4 <- 0.52 * min(max(abs(site$RDA1)) / max(abs(bp$RDA1)), max(abs(site$RDA2)) / max(abs(bp$RDA2)))
fig4A <- ggplot(site, aes(RDA1, RDA2)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey72") + geom_vline(xintercept = 0, linetype = 2, colour = "grey72") +
  geom_segment(data = bp, aes(x = 0, y = 0, xend = RDA1 * m4, yend = RDA2 * m4), inherit.aes = FALSE, colour = "grey60", linewidth = 0.6, arrow = arrow(length = unit(2.2, "mm"), type = "closed")) +
  geom_text_repel(data = bp, aes(RDA1 * m4, RDA2 * m4, label = Label), inherit.aes = FALSE, fontface = "italic", size = 3, seed = 123) +
  geom_point(aes(colour = Compartment, shape = Compartment), size = 3) + geom_text_repel(aes(label = Sample), size = 2.4, seed = 123, segment.colour = "grey65") +
  scale_colour_manual(values = zone_cols, name = NULL) + scale_shape_manual(values = zone_shapes, name = NULL) +
  labs(x = sprintf("RDA1 (%.2f%%)", axp[1]), y = sprintf("RDA2 (%.2f%%)", axp[2]), subtitle = bquote(.(paste(vlab[terms_sel], collapse = " + ")) * "; adjusted" ~ R^2 == .(sprintf("%.1f%%", 100 * r2sel)))) +
  theme_k(9) + theme(legend.position = "bottom", legend.text = element_text(size = 7), plot.subtitle = element_text(size = 7.5)) +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
lv <- c(paste0(sub(" ", "\n", vlab[terms_sel[1]]), "\nalone"), paste0(vlab[terms_sel[-1]], "\nadded"), "Unexplained")
sq <- data.frame(F = factor(lv, levels = lv), P = 100 * c(r2cum[1], diff(r2cum), 1 - r2sel))
fig4B <- ggplot(sq, aes(F, P, fill = F)) + geom_col(width = 0.68) + geom_text(aes(label = sprintf("%.1f%%", P)), vjust = -0.5, size = 2.8) +
  scale_fill_manual(values = c("#d95f02","#1b9e77","#7570b3","#e7298a")[seq_len(nrow(sq) - 1)] |> c("#d9d9d9"), guide = "none") + scale_y_continuous(limits = c(0, 80), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Community variation (%)") + theme_k(9) + theme(axis.text.x = element_text(size = 6.3))
ggsave("figures/Fig4_RDA.pdf", (fig4A | fig4B) + plot_layout(widths = c(1.35, 1)) + plot_annotation(tag_levels = "A"), width = 174, height = 115, units = "mm", device = cairo_pdf)

