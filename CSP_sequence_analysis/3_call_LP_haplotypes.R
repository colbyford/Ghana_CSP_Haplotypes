

## [5] CREATE FINAL HAPLOTYPE TABLE
## Call Haplotype Options
minCov <- 3
detectionLimit <- 1/100
minOccHap <- 2
minCovSample <- 25


## Call Final Haplotypes without SNP filtering allows for length-polymorphismus

postfixLP <- paste0(postfix, "LP")
finalTab <- createFinalHaplotypTable(outputDir = outputDir,
                                     sampleTable = procReads,
                                     markerTable = markerTab,
                                     referenceSequence = NULL,
                                     snpList = NULL,
                                     postfix = postfixLP,
                                     minHaplotypCoverage = minCov,
                                     minReplicate = minOccHap, 
                                     detectability = detectionLimit,
                                     minSampleCoverage = minCovSample)

saveRDS(finalTab, "finalTab_LP.RDS")
write.csv(do.call(rbind, finalTab), "finalTab_LP.csv", row.names = FALSE)


