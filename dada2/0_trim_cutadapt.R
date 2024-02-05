# Load packages
library("dada2")
library("ggplot2")

setwd("C:/Users/camil/Documents/TCC/covid/dada2")

# Config folders and paths
fastq = "raw_data"      # raw fastq files
trimmed = "trimmed"     # cutadapt trimmed fastq files
filt = "filt"           # dada2 trimmed fastq files
outfiles = "outfiles"   # output files
images = "images"       # output images

if(!dir.exists(trimmed)) dir.create(trimmed)
if(!dir.exists(filt)) dir.create(filt)
if(!dir.exists(outfiles)) dir.create(outfiles)
if(!dir.exists(images)) dir.create(images)

head(list.files(fastq), 20)


# Primmer trtimming -- With Cutadapt (being called by system2)
fns = sort(list.files(fastq, full.names = TRUE))
fnFs = fns[grep("R1_001.fastq.gz", fns)]
fnRs = fns[grep("R2_001.fastq.gz", fns)]

fnFs.cut = file.path(trimmed, basename(fnFs))
fnRs.cut = file.path(trimmed, basename(fnRs))
log.cut = gsub(".R1_001.fastq.gz", ".log", fnFs.cut)
sample.names = gsub(".R1_001.fastq.gz", "", basename(fnFs.cut))

# Define the primer set used to perform PCR
FWD = "GTGCCAGCMGCCGCGGTAA"       # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3753973/
REV = "TAATCTWTGGGVHCATCAGG"   # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3753973/

# Get reverse complement DNA sequences
FWD.RC = dada2::rc(FWD)
REV.RC = dada2::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags = paste("-g", FWD, "-a", REV.RC)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags = paste("-G", REV, "-A", FWD.RC)

# Run cutadapt to remove primers
cutadapt = "C:/Users/camil/AppData/Local/Programs/Python/Python312/Scripts/cutadapt"

for(i in seq_along(fnFs)) {
  print(paste0("[", i ,"/", length(sample.names), "] ", sample.names[i]))
  
  system2(cutadapt,
          stdout = log.cut[i], stderr = log.cut[i], # log file
          args = c(R1.flags, R2.flags,
                   "-n 2",                   # -n 2 required to remove FWD and REV from reads
                   "--match-read-wildcards", # enable IUPAC nucleotide codes (wildcard characters)
                   "--length 300",           # Truncate reads to 300 bp
                   "-m 150",                 # discard reads shorter than LEN (avoid length zero sequences)
                   "--overlap 10",           # min overlap between read and adapter for an adapter to be found
                   "-j 0",                   # auto-detection of CPU cores, only available on Python 3
                   "-o", fnFs.cut[i], "-p", fnRs.cut[i], # trimmed files
                   fnFs[i], fnRs[i])         # input files
  )
}



# Build trimmed file
fns = sort(list.files(trimmed, full.names = TRUE))
fnFs = fns[grep("R1_001.fastq.gz", fns)]     # Update the grep pattern when necessary
fnRs = fns[grep("R2_001.fastq.gz", fns)]     # Update the grep pattern when necessary
sample.names = gsub(".1.fastq.gz", "", basename(fnFs))  # Update the gsub pattern when necessary

# Check objects
fnFs
fnRs

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


# Filter and Trim
# Set paths to the dada2-filterd files
filtFs = file.path(filt, basename(fnFs))
filtRs = file.path(filt, basename(fnRs))

# Perform filtering and trimming
out = filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                    # Need to keep paramters consistent between runs of the same study
                    truncLen = c(240,200), minLen = 200, maxN = 0, truncQ = 2, maxEE = c(2,5),
                    rm.phix = TRUE, compress = TRUE, verbose = TRUE, multithread = TRUE)

out = as.data.frame(out)
rownames(out) = sample.names

head(out, 10)
