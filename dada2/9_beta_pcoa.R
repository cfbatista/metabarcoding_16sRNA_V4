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

ps_G2 = ps %>%
  subset_samples(groups %in% G2)

ps_G3 = ps %>%
  subset_samples(groups %in% G3)

ps_G4 = ps %>%
  subset_samples(groups %in% G4)

# ------------ Grupo 1
metadata <- as(sample_data(ps_G1), "data.frame")
dist.bc <- phyloseq::distance(ps_G1, method = "bray")
permanova_G1 <- adonis2(dist.bc ~ groups, data = metadata, perm=9999)

ordu_G1 = ordinate(ps_G1, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps_G1, ord.bray, type="samples", color="groups")+
  stat_ellipse() + 
  annotate("text", x = -2, y = 3, label = "Pr(>F)=0.48") +
  theme(legend.position = c(0.9, 0.2)) +
  theme(legend.background = element_rect(fill="lightgray", 
                                         size=0.5, linetype="solid"))
  


# ------------ Grupo 2
metadata <- as(sample_data(ps_G2), "data.frame")
dist.bc <- phyloseq::distance(ps_G2, method = "bray")
permanova_G2 <- adonis2(dist.bc ~ groups, data = metadata, perm=9999)

ordu_G2 = ordinate(ps_G2, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps_G2, ordu_G2, type="samples", color="groups")+
  stat_ellipse()


# ------------ Grupo 3
metadata <- as(sample_data(ps_G3), "data.frame")
dist.bc <- phyloseq::distance(ps_G3, method = "bray")
permanova_G3 <- adonis2(dist.bc ~ groups, data = metadata, perm=9999)

ordu_G3 = ordinate(ps_G3, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps_G3, ordu_G3, type="samples", color="groups")+
  stat_ellipse()



# ------------ Grupo 4
metadata <- as(sample_data(ps_G4), "data.frame")
dist.bc <- phyloseq::distance(ps_G4, method = "bray")
permanova_G4 <- adonis2(dist.bc ~ groups, data = metadata, perm=9999)

ordu_G4 = ordinate(ps_G4, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps_G4, ordu_G4, type="samples", color="groups")+
  stat_ellipse()

