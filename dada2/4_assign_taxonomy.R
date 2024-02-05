# Assign Taxonomy
taxaSilva <- assignTaxonomy(seqTab.nochim, 
                            'C:/Users/camil/Documents/databases/silva_nr99_v138.1_train_set.fa.gz',  multithread = TRUE)


taxaSilva <- addSpecies(taxaSilva, 'C:/Users/camil/Documents/databases/silva_species_assignment_v138.1.fa.gz')

#tempo execução 40 m
taxaGTDB <- assignTaxonomy(seqTab.nochim, 
                           'C:/Users/camil/Documents/databases/GTDB_bac120_arc53_ssu_r207_fullTaxo.fa.gz',
                           multithread = TRUE)

taxaGTDB <- addSpecies(taxaGTDB, 'C:/Users/camil/Documents/databases/GTDB_dada2_assignment_species.fa.gz')


# tempo execução 10 m
taxaRDP <- assignTaxonomy(seqTab.nochim, 
                          'C:/Users/camil/Documents/databases/rdp_train_set_18.fa.gz',
                          multithread = TRUE)

taxaRDP <- addSpecies(taxaRDP, 'C:/Users/camil/Documents/databases/rdp_species_assignment_18.fa.gz')


#taxaGG <- assignTaxonomy(seqTab.nochim, '~silva_nr99_v138.1_train_set.fa.gz',multithread = TRUE)


taxaSilva.print <- taxaSilva
rownames(taxaSilva.print) <- NULL
head(taxaSilva.print)


taxaRDP.print <- taxaRDP
rownames(taxaRDP.print) <- NULL
head(taxaRDP.print)

taxaGTDB.print <- taxaGTDB
rownames(taxaGTDB.print) <- NULL
head(taxaGTDB.print)
