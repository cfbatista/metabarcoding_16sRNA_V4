library("DECIPHER")
library("phangorn")

# Align Sequences
sequences <- getSequences(seqtab.nochim)
names(sequences) <- sequences
alignmnet <- AlignSeqs(DNAStringSet(seqs))

#change sequence alignment output unto phydat structure
phang.align <- phyDat(as(alignmnet, "matrix"), type="DNA")
write.phyDat(phang.align, file = "alignment.fasta", format = "fasta")
write.phyDat(phang.align, file = "alignment.aln", format = "phylip")

#Create distance matrix
dm <- dist.ml(phang.align)

#Perform Neighbor joining
treeNJ <- NJ(dm)

#Internal maximum likehood
fit <- pml(treeNJ, data=phang.align)

filtGTR <- update(fit, k=4, inv=0.2)
filtGTR <- optim.pml(filtGTR, model="GTR", optInv = TRUE, optGamma = TRUE, rearrangement = "stochastic", control=pml.control(trace=0))

