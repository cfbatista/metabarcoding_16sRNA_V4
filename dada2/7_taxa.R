library("phyloseq")
library("DESeq2")
library("ggplot2")
library("ggbeeswarm")
library("ggrepel")
library("vegan")
library("tidyverse")
library("data.table")

ps

# Create a total counts data.table
tdt = data.table(setDT(as.data.frame(tax_table(ps))), 
                 TotalCounts = taxa_sums(ps), SV = taxa_names(ps))

tdt


ggplot(tdt, aes(TotalCounts)) + geom_histogram(bins = 50) + theme_bw() + 
  ggtitle("Histogram of Total Counts")

tdt[(TotalCounts == 1), .N]     # singletons
tdt[(TotalCounts == 2), .N]     # doubletons



# taxa cumulative sum
taxcumsum = tdt[, .N, by = TotalCounts]
setkey(taxcumsum, TotalCounts)
taxcumsum[, CumSum := cumsum(N)]

# Plot the cumulative sum of ASVs against the total counts
pCumSum = ggplot(taxcumsum, aes(TotalCounts, CumSum)) + geom_point() + theme_bw() + 
  xlab("Filtering Threshold") + ylab("ASV Filtered")

gridExtra::grid.arrange(pCumSum, pCumSum + xlim(0, 500), 
                        pCumSum + xlim(0, 100), pCumSum + xlim(0, 50), nrow = 2, 
                        top = "ASVs that would be filtered vs. minimum taxa counts threshold")



#Prevalence table
fast_melt = function(physeq,
                     includeSampleVars = character(),
                     omitZero = FALSE){
  require("phyloseq")
  require("data.table")
  # supports "naked" otu_table as `physeq` input.
  otutab = as(otu_table(physeq), "matrix")
  if(!taxa_are_rows(physeq)){otutab <- t(otutab)}
  otudt = data.table(otutab, keep.rownames = TRUE)
  setnames(otudt, "rn", "TaxaID")
  # Enforce character TaxaID key
  otudt[, TaxaIDchar := as.character(TaxaID)]
  otudt[, TaxaID := NULL]
  setnames(otudt, "TaxaIDchar", "TaxaID")
  # Melt count table
  mdt = melt.data.table(otudt, 
                        id.vars = "TaxaID",
                        variable.name = "SampleID",
                        value.name = "count")
  if(omitZero){
    # Omit zeroes and negative numbers
    mdt <- mdt[count > 0]
  }
  # Omit NAs
  mdt <- mdt[!is.na(count)]
  # Calculate relative abundance
  mdt[, RelativeAbundance := count / sum(count), by = SampleID]
  if(!is.null(tax_table(physeq, errorIfNULL = FALSE))){
    # If there is a tax_table, join with it. Otherwise, skip this join.
    taxdt = data.table(as(tax_table(physeq, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
    setnames(taxdt, "rn", "TaxaID")
    # Enforce character TaxaID key
    taxdt[, TaxaIDchar := as.character(TaxaID)]
    taxdt[, TaxaID := NULL]
    setnames(taxdt, "TaxaIDchar", "TaxaID")
    # Join with tax table
    setkey(taxdt, "TaxaID")
    setkey(mdt, "TaxaID")
    mdt <- taxdt[mdt]
  }
  # includeSampleVars = c("DaysSinceExperimentStart", "SampleType")
  # includeSampleVars = character()
  # includeSampleVars = c()
  # includeSampleVars = c("aksjdflkas") 
  wh.svars = which(sample_variables(physeq) %in% includeSampleVars)
  if( length(wh.svars) > 0 ){
    # Only attempt to include sample variables if there is at least one present in object
    sdf = as(sample_data(physeq), "data.frame")[, wh.svars, drop = FALSE]
    sdt = data.table(sdf, keep.rownames = TRUE)
    setnames(sdt, "rn", "SampleID")
    # Join with long table
    setkey(sdt, "SampleID")
    setkey(mdt, "SampleID")
    mdt <- sdt[mdt]
  }
  setkey(mdt, "TaxaID")
  return(mdt)
}
mdt = fast_melt(ps)

mdt


prevdt = mdt[, list(Prevalence = sum(count > 0), TotalCounts = sum(count), 
                    MaxCounts = max(count)), by = TaxaID]

prevdt

ggplot(prevdt, aes(Prevalence)) + geom_histogram(bins = 50) + theme_bw() +
  ggtitle("Histogram of Taxa Prevalence")

prevdt[(Prevalence == 1), .N]   # singletons
prevdt[(Prevalence == 2), .N]   # doubletons

ggplot(prevdt, aes(MaxCounts)) + geom_histogram(bins = 100) + theme_bw() + 
  ggtitle("Histogram of Maximum TotalCounts")

table(prevdt$MaxCounts)[1:50]

prevdt[(MaxCounts == 1)]    # singletons


# taxa cumulative sum
prevcumsum = prevdt[, .N, by = Prevalence]
setkey(prevcumsum, Prevalence)
prevcumsum[, CumSum := cumsum(N)]

# Plot the cumulative sum of ASVs against the prevalence
pPrevCumSum = ggplot(prevcumsum, aes(Prevalence, CumSum)) + geom_point(size = 2, alpha = 0.5) + 
  theme_bw() + xlab("Filtering Threshold") + ylab("ASVs Filtered") + 
  ggtitle("ASVs that would be filtered vs. minimum sample count threshold")

pPrevCumSum

# Prevalence vs. Total Count Scatter plot
ggplot(prevdt, aes(Prevalence, TotalCounts)) + geom_point(size = 2, alpha = 0.5) + 
  scale_y_log10() + theme_bw() + xlab("Prevalence [No. Samples]") + ylab("TotalCounts [Taxa]")


