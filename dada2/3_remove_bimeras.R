# Remove Quimeras
seqtab.nochim = removeBimeraDenovo(seqtab, method = "consensus", 
                                   multithread = TRUE, 
                                   verbose = TRUE)

# Update the same names to exclude "fastq.gz" in name
rownames(seqtab.nochim) = sample.names
dim(seqtab)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) # Print the proportion of non-chimeras in merged sequence reads


# Track Reads through the pipeline
getN <- function(x) sum(getUniques(x))

track = cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), 
              sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) = c("Trimmed", "Filtered", "denoisedF", "denoisedR", "merged", "nonchim")
track = cbind(data.frame(SampleID = sample.names), track)

head(track)

