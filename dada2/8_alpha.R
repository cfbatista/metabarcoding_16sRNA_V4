library(ggplot2)
library(ggpubr)
library(phyloseq)

G1 = c('ARI+', 'SARS+')
G2 = c('ARI-', 'SARS-')
G3 = c('ARI+', 'ARI-')
G4 = c('SARS+', 'SARS-')

ps_G1 = ps %>%
  subset_samples(groups %in% G1)

ps_G2 = ps %>%
  subset_samples(groups %in% G2)

ps_G3 = ps %>%
  subset_samples(groups %in% G3)

ps_G4 = ps %>%
  subset_samples(groups %in% G4)

# ------------------------------- Wilcoxon test

rich = estimate_richness(ps_G4_prune)
pairwise.wilcox.test(rich$Shannon, sample_data(ps_G4_prune)$groups_beta)


# ------------------------------- GRUPO 1
ps_G1_prune <- prune_taxa(taxa_sums(ps_G1) > 0, ps_G1)
tab1 <- microbiome::alpha(ps_G1_prune, index = "shannon");
tab2 <- microbiome::alpha(ps_G1_prune, index = "simpson");
tab3 <- microbiome::alpha(ps_G1_prune, index = "chao1");

ps_G1_prune.meta <- microbiome::meta(ps_G1);
ps_G1_prune.meta$Shannon <- tab1$diversity_shannon
ps_G1_prune.meta$Simpson <- tab2$dominance_simpson
ps_G1_prune.meta$chao1 <- tab3$chao1

bmi_G1 <- levels(ps_G1_prune.meta$groups)
bmi_G1.pairs_group <- combn(seq_along(bmi_G1), 2, simplify = FALSE, FUN = function(i)bmi_G1[i])

plot_shannon_G1 <- ggboxplot(ps_G1_prune.meta, x = "groups", y = "Shannon",
                            color = "black", palette =c("#E7B800", "#9370DB"), fill="groups",
                            facet.by = "age", short.panel.labs = FALSE) + geom_jitter(width = .3, shape = 21, color = "black")
plot_shannon_G1 <- plot_shannon_G1 + stat_compare_means(comparisons = bmi_G1.pairs_group, aes(label = paste0("Wilcoxon, p = ", after_stat(p.format))))
plot_shannon_G1



# ------------------------------- GRUPO 2
ps_G2_prune <- prune_taxa(taxa_sums(ps_G2) > 0, ps_G2)
tab1 <- microbiome::alpha(ps_G2_prune, index = "shannon");
tab2 <- microbiome::alpha(ps_G2_prune, index = "simpson");
tab3 <- microbiome::alpha(ps_G2_prune, index = "chao1");

ps_G2_prune.meta <- microbiome::meta(ps_G2);
ps_G2_prune.meta$Shannon <- tab1$diversity_shannon
ps_G2_prune.meta$Simpson <- tab2$dominance_simpson
ps_G2_prune.meta$chao1 <- tab3$chao1

bmi_G2 <- levels(ps_G2_prune.meta$groups)
bmi_G2.pairs_group <- combn(seq_along(bmi_G2), 2, simplify = FALSE, FUN = function(i)bmi_G2[i])


plot_shannon_G2 <- ggboxplot(ps_G2_prune.meta, x = "groups", y = "Shannon",
                             color = "black", palette =c("#7CCD7C", "palevioletred2"), fill="groups",
                             facet.by = "age", short.panel.labs = FALSE) + geom_jitter(width = .3, shape = 21, color = "black")
plot_shannon_G2 <- plot_shannon_G2 + stat_compare_means(comparisons = bmi_G2.pairs_group, aes(label = paste0("Wilcoxon, p = ", after_stat(p.format))))
plot_shannon_G2



# ------------------------------- GRUPO 3
ps_G3_prune <- prune_taxa(taxa_sums(ps_G3) > 0, ps_G3)
tab1 <- microbiome::alpha(ps_G3_prune, index = "shannon");
tab2 <- microbiome::alpha(ps_G3_prune, index = "simpson");
tab3 <- microbiome::alpha(ps_G3_prune, index = "chao1");

ps_G3_prune.meta <- microbiome::meta(ps_G3);
ps_G3_prune.meta$Shannon <- tab1$diversity_shannon
ps_G3_prune.meta$Simpson <- tab2$dominance_simpson
ps_G3_prune.meta$chao1 <- tab3$chao1

bmi_G3 <- levels(ps_G3_prune.meta$groups)
bmi_G3.pairs_group <- combn(seq_along(bmi_G3), 2, simplify = FALSE, FUN = function(i)bmi_G3[i])


plot_shannon_G3 <- ggboxplot(ps_G3_prune.meta, x = "groups", y = "Shannon",
                             color = "black", palette =c("#E7B800","#7CCD7C"), fill="groups",
                             facet.by = "age", short.panel.labs = FALSE) + geom_jitter(width = .3, shape = 21, color = "black")
plot_shannon_G3 <- plot_shannon_G3 + stat_compare_means(comparisons = bmi_G3.pairs_group, aes(label = paste0("Wilcoxon, p = ", after_stat(p.format))))
plot_shannon_G3


# ------------------------------- GRUPO 4
ps_G4_prune <- prune_taxa(taxa_sums(ps_G4) > 0, ps_G4)
tab1 <- microbiome::alpha(ps_G4_prune, index = "shannon");
tab2 <- microbiome::alpha(ps_G4_prune, index = "simpson");
tab3 <- microbiome::alpha(ps_G4_prune, index = "chao1");

ps_G4_prune.meta <- microbiome::meta(ps_G4);
ps_G4_prune.meta$Shannon <- tab1$diversity_shannon
ps_G4_prune.meta$Simpson <- tab2$dominance_simpson
ps_G4_prune.meta$chao1 <- tab3$chao1

bmi_G4 <- levels(ps_G4_prune.meta$groups)
bmi_G4.pairs_group <- combn(seq_along(bmi_G4), 2, simplify = FALSE, FUN = function(i)bmi_G4[i])


plot_shannon_G4 <- ggboxplot(ps_G4_prune.meta, x = "groups", y = "Shannon",
                             color = "black", palette =c("#9370DB","palevioletred2"), fill="groups",
                             facet.by = "age", short.panel.labs = FALSE) + geom_jitter(width = .3, shape = 21, color = "black")
plot_shannon_G4 <- plot_shannon_G4 + stat_compare_means(comparisons = bmi_G4.pairs_group, aes(label = paste0("Wilcoxon, p = ", after_stat(p.format))))
plot_shannon_G4

