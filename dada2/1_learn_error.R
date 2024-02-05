#Learn the error rates
# the learnErrors function to performs dereplication and learn the error rates. 

errF = learnErrors(filtFs, multithread = TRUE)
errR = learnErrors(filtRs, multithread = TRUE)

# Visualize the estimated error rates
pdf(paste0(images, "/plotErrors.pdf"), width = 10, height = 10, pointsize = 12)
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
invisible(dev.off())
