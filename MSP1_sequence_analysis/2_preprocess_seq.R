
## Load in Libraries
library("HaplotypR")
library("ShortRead")
library(dplyr)

## Set Paths to Required Files
primerFile <- "meta_data/marker_msp1.txt"
sampleFile <- "meta_data/demultiplexSampleSummary.txt"

## Set Output Directory
outputDir <- "output"

## [1] DEMULTIPLEXING - Performed Outside of HaplotypR
dePlexSample <- read.table(sampleFile, sep="\t", header = T, stringsAsFactors = F)
is.data.frame(dePlexSample)

## [2] CREATE DEPLEX TABLE
## Create Output Subdirectory 
outDeplexMarker <- file.path(outputDir, "dePlexMarker")
dir.create(outDeplexMarker)

## Process each marker
markerTab <- read.delim(primerFile, stringsAsFactors=F)
rownames(markerTab) <- markerTab$MarkerID
dePlexMarker <- demultiplexByMarker(dePlexSample %>% filter(!SampleID %in% c("Yakabu-1_S13", "Yakabu-2_S12")),
                                    markerTab,
                                    outDeplexMarker,
                                    trimFilenameExt = "R1_001\\.fastq.gz$")

## Filter to fix cross-contamination issues
dePlexMarker <- dePlexMarker %>% filter(numReadOut > 1000) #100

# save summary table
write.table(dePlexMarker, file.path(outputDir, "demultiplexMarkerSummary.txt"), sep="\t", row.names=F, quote=F)


#############################################################
# remove NA samples
dePlexMarker <- dePlexMarker[!is.na(dePlexMarker$FileR1),]

## [3] BIND AMPLICON READS
## Create Output Subdirectory
outProcFiles <- file.path(outputDir, "processedReads")
dir.create(outProcFiles)

# postfix <- "_NGmerge"
# refSeq <- DNAStringSet(markerTab$ReferenceSequence)
# names(refSeq) <- markerTab$MarkerID
# lapply(seq_along(refSeq), function(i){
#   writeFasta(refSeq[i], file.path(outputDir, paste(names(refSeq)[i], postfix, ".fasta", sep="")))
# })

# procReads <- mergeAmpliconReads(fastqFileR1=dePlexMarker$FileR1,
#                                 fastqFileR2=dePlexMarker$FileR2,
#                                 outputDir=outProcFiles,
#                                 method="NGmerge",
#                                 # method="vsearch",
#                                 mergePrefix=postfix)

# procReads <- cbind(dePlexMarker[,c("SampleID", "SampleName","BarcodePair", "MarkerID")], procReads)

# write.table(procReads, file.path(outputDir, sprintf("processedReadSummary%s.txt", postfix)), sep="\t", row.names=F, quote=F)
# #procReads <- read.delim(file.path(outputDir, sprintf("processedReadSummary%s.txt", postfix)), stringsAsFactors = F)

# Trim options
numNtF <- 180 #220 #150 # Adjust to your needs
numNtR <- 180 #220 #150 # Adjust to your needs
postfix <- sprintf("_bind%.0f_%.0f", numNtF, numNtR)

# Adjust reference to trim options and save as fasta file
refSeq <- as.character(markerTab$ReferenceSequence)
refSeq <- DNAStringSet(paste(substr(refSeq, 1,numNtF), substr(refSeq, nchar(refSeq)+1-numNtR, nchar(refSeq)), sep=""))
names(refSeq) <- markerTab$MarkerID
lapply(seq_along(refSeq), function(i){
  writeFasta(refSeq[i], file.path(outputDir, paste(names(refSeq)[i], postfix, ".fasta", sep="")))
})

# Fuse paired read
procReads <- bindAmpliconReads(as.character(dePlexMarker$FileR1),
                               as.character(dePlexMarker$FileR2),
                               outProcFiles, 
                               read1Length=numNtF, read2Length=numNtR)

procReads <- cbind(dePlexMarker[,c("SampleID", "SampleName","BarcodePair", "MarkerID")], procReads)

write.table(procReads, file.path(outputDir, sprintf("processedReadSummary%s.txt", postfix)), sep="\t", row.names=F)