############################################################
## Korbous - Fig. 2 (alpha diversity by compartment, descriptive), Fig. S1 (rarefaction), Table S1
## Run from the folder containing 00_common.R and data/
############################################################
source("00_common.R")

## ======================= ALPHA DIVERSITY (Fig. 2, Fig. S1, Table S1) =======================
S1f <- read.csv("data/TableS1_alpha_diversity_Fabrice.csv", row.names = 1)[sample_order, ]
alpha <- data.frame(Sample = sample_order, Compartment = meta$Zone,
  Observed_ASVs = rowSums(rar > 0),
  Shannon = round(diversity(rar, "shannon"), 3),
  Gini_Simpson = round(diversity(rar, "simpson"), 3),
  PD_obs = S1f$PD_obs)
stopifnot(all(alpha$Observed_ASVs == S1f$Observed_ASVs), max(abs(alpha$Shannon - S1f$Shannon)) < 0.002)
write.csv(alpha, "tables/TableS1_alpha_diversity.csv", row.names = FALSE)

## Compartments are compared descriptively: samples of a compartment come from one or two sites and are not
## independent replicates (for all four indices, terrestrial sediments do not overlap with any other sample)
idx <- c(Observed_ASVs = "Observed richness", Shannon = "Shannon diversity", Gini_Simpson = "Gini-Simpson index", PD_obs = "Phylogenetic diversity")
al <- alpha |> select(Compartment, Observed_ASVs, Shannon, Gini_Simpson, PD_obs) |>
  pivot_longer(-Compartment, names_to = "Index", values_to = "Value") |>
  mutate(Index = factor(idx[Index], levels = idx))
fig2 <- ggplot(al, aes(Compartment, Value)) +
  geom_boxplot(data = filter(al, Compartment != "Marine reference"), aes(fill = Compartment, colour = Compartment),
               outlier.shape = NA, width = 0.6, alpha = 0.25, linewidth = 0.4) +
  geom_jitter(aes(colour = Compartment), width = 0.12, size = 1.6, alpha = 0.9) +
  facet_wrap(~Index, scales = "free_y", labeller = as_labeller(setNames(paste0(LETTERS[1:4], ". ", idx), idx))) +
  scale_fill_manual(values = zone_cols) + scale_colour_manual(values = zone_cols) +
  scale_x_discrete(limits = zones, labels = function(x) gsub(" ", "\n", x)) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  labs(x = NULL, y = "Alpha diversity") +
  theme_bw(base_size = 8) + theme(legend.position = "none", axis.text.x = element_text(size = 6.5),
    strip.text = element_text(face = "bold"), strip.background = element_rect(fill = "grey95", colour = "grey70"))
ggsave("figures/Fig2_alpha_diversity.pdf", fig2, width = 174, height = 152, units = "mm", device = cairo_pdf)

rc <- do.call(rbind, lapply(sample_order, function(s) { n <- sum(raw[s, ]); st <- unique(c(seq(1, n, by = 500), n))
  data.frame(Sample = s, Depth = st, ASVs = as.numeric(rarefy(raw[s, , drop = FALSE], st))) }))
rc$Compartment <- meta[rc$Sample, "Zone"]
figS1 <- ggplot(rc, aes(Depth, ASVs, group = Sample, colour = Compartment)) + geom_line(linewidth = 0.5) +
  geom_vline(xintercept = unique(rowSums(rar)), linetype = "dashed", colour = "grey30") +
  scale_colour_manual(values = zone_cols, name = NULL) + scale_x_continuous(labels = label_number(big.mark = ",")) +
  labs(x = "Sequencing depth (reads)", y = "Observed ASV richness") + theme_k() + theme(legend.position = "bottom") +
  guides(colour = guide_legend(ncol = 2))
ggsave("figures/FigS1_rarefaction.pdf", figS1, width = 130, height = 100, units = "mm", device = cairo_pdf)

