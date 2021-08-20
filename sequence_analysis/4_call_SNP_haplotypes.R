
postfixSNP <- paste0(postfix, "SNP")

## [4] CREATE SNP LIST
minMMrate <- 0.50
minOccGen <- 3 #2

## Process Each Marker
snpLst <- lapply("csp", function(marker){
  # Calculate mismatch rate
  seqErrLst <- calculateMismatchFrequencies(as.character(procReads[procReads$MarkerID == marker, "ReadFile"]), 
                                            refSeq[marker], 
                                            method ="pairwiseAlignment", # c("pairwiseAlignment","compareDNAString"), 
                                            minCoverage=100L)
  names(seqErrLst) <- procReads[procReads$MarkerID == marker, "SampleID"]
  seqErr <- do.call(cbind, lapply(seqErrLst, function(l){
    l[,"MisMatch"]/l[,"Coverage"]
  }))
  write.table(seqErr, file.path(outputDir, sprintf("mismatchRate_rate_%s%s.txt", marker, postfixSNP)), sep="\t", row.names=F)
  
  # Call SNPs
  potSNP <- callGenotype(seqErr, minMismatchRate=minMMrate, minReplicate=minOccGen)
  snpRef <- unlist(lapply(potSNP, function(snp){
    as.character(subseq(refSeq[marker], start=snp, width=1))
  }))
  snps <- data.frame(Chr=marker, Pos=potSNP, Ref=snpRef, Alt="N", stringsAsFactors=F)
  write.table(snps, file=file.path(outputDir, sprintf("potentialSNPlist_rate%.0f_occ%i_%s%s.txt", 
                                                      minMMrate*100, minOccGen, marker, postfixSNP)), 
              row.names=F, col.names=T, sep="\t", quote=F)
  
  # Plot mismatch rate and SNP calls
  png(file.path(outputDir, sprintf("plotMisMatchRatePerBase_rate%.0f_occ%i_%s%s.png", 
                                   minMMrate*100, minOccGen, marker, postfixSNP)), 
      width=1500 , height=600)
  matplot(seqErr, type="p", pch=16, cex=0.4, col="#00000088", ylim=c(0, 1),
          ylab="Mismatch Rate", xlab="Base Position", main=marker, cex.axis=2, cex.lab=2)
  abline(v=snps[,"Pos"], lty=2, col="grey")
  abline(h=minMMrate, lty=1, col="red")
  if(!is.null(refSNP[marker])){
    axis(1, at=refSNP[[marker]][,"Pos"], line=-0.5, col=0, col.ticks=2, lwd.ticks=2, labels = NA, outer=F)
    axis(3, at=refSNP[[marker]][,"Pos"], line=-0.5, col=0, col.ticks=2, lwd.ticks=2, labels = NA)
  }
  dev.off()
  
  return(snps)
})

names(snpLst) <- "csp"
saveRDS(snpLst, "snpLst_csp.RDS")

## Add in reference SNPs
refSNP <- list(csp=read.delim("meta_data/refSNP_csp.txt"))

snpLst[["csp"]] <- rbind(snpLst[["csp"]],
                         refSNP[["csp"]] %>% 
                           filter(!Pos %in% snpLst[["csp"]]$Pos)) %>%
  arrange(Pos)

# ## fuse new with ref SNPs
# if(!is.null(refSNP[marker])){
#   snpLst <- lapply(markerTab$MarkerID, function(marker){
#     # mIdx <- refSNP[[marker]]$Pos %in% potSNPs[[marker]]$Pos
#     # snps <- rbind(potSNPs[[marker]], refSNP[[marker]][!mIdx,])
#     mIdx <- potSNPs[[marker]]$Pos %in% refSNP[[marker]]$Pos
#     snps <- rbind(refSNP[[marker]], potSNPs[[marker]][!mIdx,])
#     snps <- snps[order(snps$Pos),]
#     rownames(snps) <- NULL
#     return(snps)
#   })
#   names(snpLst) <- markerTab$MarkerID
# }

## Call Final Haplotypes WITH SNP filtering DON'T allow for length-polymorphismus
## [5] CREATE FINAL HAPLOTYPE TABLE
## Call Haplotype Options
minCov <- 3
detectionLimit <- 1/50 #1/100
minOccHap <- 3 #2
minCovSample <- 25


## Call Final Haplotypes
finalTab <- createFinalHaplotypTable(outputDir = outputDir,
                                     sampleTable = procReads[procReads$MarkerID=="csp",],
                                     markerTable = markerTab["csp",],
                                     referenceSequence = refSeq["csp",],
                                     snpList = snpLst["csp"],
                                     postfix = postfixSNP,
                                     minHaplotypCoverage = minCov,
                                     minReplicate = minOccHap, 
                                     detectability = detectionLimit,
                                     minSampleCoverage = minCovSample)

saveRDS(finalTab, "finalTab_csp.RDS")
write.csv(do.call(rbind, finalTab), "finalTab_csp.csv", row.names = FALSE)
