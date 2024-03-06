library("phyloseq")
library("qiime2R")

#Load sample information
metadata = read_q2metadata("sample-metadata.tsv")

# Garantir que as samples fiquem com o mesmo nome
samples.out <- rownames(seqtab.nochim)
rownames(metadata) <- samples.out

#Conferir os nomes
sample_names(otu_table(seqtab.nochim, taxa_are_rows=FALSE))
sample_names(sample_data(metadata))
sample_names(tax_table(taxaSilva))



ps = phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), 
              sample_data(metadata), 
              tax_table(taxaSilva),
              phy_tree(filtGTR$tree))

ps <- merge_phyloseq(ps,metadata)
ps


# rooted tree
set.seed(711)
phy_tree(ps) <- root(phy_tree(ps), sample(taxa_names(ps), 1), resolve.root = TRUE)
is.rooted(phy_tree(ps))
