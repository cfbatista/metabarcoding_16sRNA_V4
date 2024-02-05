# Sample Inference
dadaFs = dada(filtFs, err = errF, pool = "pseudo", multithread = TRUE)
dadaRs = dada(filtRs, err = errR, pool = "pseudo", multithread = TRUE)

#Merge paired reads
mergers = mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

#Construct sequence table
seqtab = makeSequenceTable(mergers)

table(nchar(getSequences(seqtab)))
head(mergers[[1]])

head(track)
saveRDS(seqtab, "seqtab.rds")
