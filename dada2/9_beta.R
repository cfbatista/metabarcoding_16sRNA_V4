library(ggplot2)
library(ggpubr)
library(phyloseq)
library(vegan)
library(tidyverse)

G1 = c('ARI+', 'SARS+')
G2 = c('ARI-', 'SARS-')
G3 = c('ARI+', 'ARI-')
G4 = c('SARS+', 'SARS-')

ps_G1 = ps %>%
  subset_samples(groups %in% G1)
ps_G1 <- prune_taxa(taxa_sums(ps_G1) > 0, ps_G1)
ps_G2 = ps %>%
  subset_samples(groups %in% G2)

ps_G3 = ps %>%
  subset_samples(groups %in% G3)

ps_G4 = ps %>%
  subset_samples(groups %in% G4)


# PERMANOVA

metadata <- as(sample_data(ps), "data.frame")
dist.bc <- phyloseq::distance(ps, method = "bray")
permanova <- adonis2(dist.bc ~ groups_beta, data = metadata, perm=9999)
# --------------------------------------------------------------------

# ------------ Grupo 1 - Genus
ps_G1_genus = tax_glom(ps_G1, 'Genus')

metadata <- as(sample_data(ps_G1_genus), "data.frame")
dist.bc <- phyloseq::distance(ps_G1_genus, method = "bray")
permanova_G1_genus <- adonis2(dist.bc ~ groups_beta, data = metadata, perm=9999)

ord.bray_G1_genus <- ordinate(ps_G1_genus, "MDS", "bray");

ps_G1_genus_plt <- plot_ordination(ps_G1_genus, ord.bray_G1_genus, type='Sample', color="groups_beta") +
  guides(color = guide_legend(title = "Groups")) +
  annotate("text", x = -0.32, y = 0.5, label = "Pr(>F)=0.59") +
  annotate("text", x = 0.45, y = 0.5, label = "Genus") +
  geom_line() + geom_point(size=2)

# ------------ Grupo 1 - Species
ps_G1_specie = tax_glom(ps_G1, 'Species')

metadata <- as(sample_data(ps_G1_specie), "data.frame")
dist.bc <- phyloseq::distance(ps_G1_specie, method = "bray")
permanova_G1_specie <- adonis2(dist.bc ~ groups_beta, data = metadata, perm=9999)

ord.bray_G1_specie <- ordinate(ps_G1_specie, "MDS", "bray");

ps_G1_genus_plt <- plot_ordination(ps_G1_specie, ord.bray_G1_specie, type='Sample', color="groups_beta") +
  guides(color = guide_legend(title = "Groups")) +
  annotate("text", x = -0.32, y = 0.5, label = "Pr(>F)=0.58") +
  annotate("text", x = 0.45, y = 0.5, label = "Genus") +
  geom_polygon(aes(fill=groups_beta), fill=NA) + geom_point(size=2)

# ------------ Grupo 2 - Genus
ps_G2_genus = tax_glom(ps_G2, 'Genus')

metadata <- as(sample_data(ps_G2_genus), "data.frame")
dist.bc <- phyloseq::distance(ps_G2_genus, method = "bray")
permanova_G2_genus <- adonis2(dist.bc ~ groups_beta, data = metadata, perm=9999)

ord.bray_G2_genus <- ordinate(ps_G2_genus, "MDS", "bray");

ps_G2_genus_plt <- plot_ordination(ps_G2_genus, ord.bray_G2_genus, type='Sample', color="groups_beta") +
  guides(color = guide_legend(title = "Groups")) +
  annotate("text", x = -0.32, y = 0.5, label = "Pr(>F)=0.01") +
  annotate("text", x = 0.35, y = 0.5, label = "Genus") +
  geom_path() + geom_point(size=2)

# ------------ Grupo 2 - Species


# ------------ Grupo 3 - Genus
ps_G3_genus = tax_glom(ps_G3, 'Genus')

metadata <- as(sample_data(ps_G3_genus), "data.frame")
dist.bc <- phyloseq::distance(ps_G3_genus, method = "bray")
permanova_G3_genus <- adonis2(dist.bc ~ groups_beta, data = metadata, perm=9999)

ord.bray_G3_genus <- ordinate(ps_G3_genus, "MDS", "bray");

ps_G3_genus_plt <- plot_ordination(ps_G3_genus, ord.bray_G3_genus, type='Sample', color="groups_beta") +
  guides(color = guide_legend(title = "Groups")) +
  annotate("text", x = -0.32, y = 0.5, label = "Pr(>F)=0.09") +
  annotate("text", x = 0.35, y = 0.5, label = "Genus") +
  geom_polygon(aes(fill=groups_beta), fill=NA) + geom_point(size=2)

# ------------ Grupo 3 - Species


# ------------ Grupo 4 - Genus
ps_G4_genus = tax_glom(ps_G4, 'Genus')

metadata <- as(sample_data(ps_G4_genus), "data.frame")
dist.bc <- phyloseq::distance(ps_G4_genus, method = "bray")
permanova_G4_genus <- adonis2(dist.bc ~ groups_beta, data = metadata, perm=9999)

ord.bray_G4_genus <- ordinate(ps_G4_genus, "MDS", "bray");

ps_G4_genus_plt <- plot_ordination(ps_G4_genus, ord.bray_G4_genus, type='Sample', color="groups_beta") +
  guides(color = guide_legend(title = "Groups")) +
  annotate("text", x = -0.20, y = 0.5, label = "Pr(>F)=0.04") +
  annotate("text", x = 0.45, y = 0.5, label = "Genus") +
  geom_polygon(aes(fill=groups_beta), fill=NA) + geom_point(size=2)

# ------------ Grupo 4 - Species



