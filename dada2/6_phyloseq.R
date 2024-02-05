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
sample_names(tax_table(taxtab))

ps = phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), 
              sample_data(metadata), 
              tax_table(spec_silva))

ps
