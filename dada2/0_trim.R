library("dada2")
library("ggplot2")

setwd("C:/Users/camil/Documents/TCC/covid/dada2")

# Config folders and paths
fastq = "raw_data"

head(list.files(fastq), 20)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Inspect Quality of trimmed files
# Plot quality profile of fastq files
ii = 1:length(sample.names)
pdf(paste0(images, "/plotQualityProfile.pdf"), width = 8, height = 8, pointsize = 12)
for(i in ii) {
  message(paste0("[", i ,"/", length(sample.names), "] ", sample.names[i]))
  print(plotQualityProfile(fnFs[i]) + ggtitle("Fwd"))
  print(plotQualityProfile(fnRs[i]) + ggtitle("Rev"))
}
invisible(dev.off())



# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)