
## Load in Libraries
library("HaplotypR")
library("ShortRead")
library(dplyr)

## Set Paths to Required Files
primerFile <- "meta_data/marker_csp_msp1.txt"
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
# dePlexMarker <- demultiplexByMarker(dePlexSample %>% filter(MarkerID == "csp", SampleID != "Yakabu-2_S12"),
#                                     markerTab %>% filter(MarkerID == "csp"),
#                                     outDeplexMarker,
#                                     trimFilenameExt = "R1_001\\.fastq.gz$")

# # Remove wrong marker file pair, due cross-contamination
# library(stringr)
# dePlexMarker$RealMarker <- do.call(rbind, strsplit(dePlexMarker$SampleID," "))[,2]
# # dePlexMarker$RealMarker <- str_to_lower(dePlexMarker$RealMarker)
# idx <- substring(dePlexMarker$MarkerID,1,4) == dePlexMarker$RealMarker
# dePlexMarker <- dePlexMarker[idx,]

## Filter to fix cross-contamination issues
dePlexMarker <- dePlexMarker %>% filter(numReadOut > 1000) #100

# save summary table
write.table(dePlexMarker, file.path(outputDir, "demultiplexMarkerSummary.txt"), sep="\t", row.names=F, quote=F)


#############################################################
# remove NA samples
dePlexMarker <- dePlexMarker[!is.na(dePlexMarker$FileR1),]

## [3] MERGE AMPLICON READS
## Create Output Subdirectory
outProcFiles <- file.path(outputDir, "processedReads")
dir.create(outProcFiles)

postfix <- "_NGmerge"
refSeq <- DNAStringSet(markerTab$ReferenceSequence)
names(refSeq) <- markerTab$MarkerID
lapply(seq_along(refSeq), function(i){
  writeFasta(refSeq[i], file.path(outputDir, paste(names(refSeq)[i], postfix, ".fasta", sep="")))
})

procReads <- mergeAmpliconReads(fastqFileR1=dePlexMarker$FileR1,
                                fastqFileR2=dePlexMarker$FileR2,
                                outputDir=outProcFiles,
                                method="NGmerge",
                                # method="vsearch",
                                mergePrefix=postfix)

procReads <- cbind(dePlexMarker[,c("SampleID", "SampleName","BarcodePair", "MarkerID")], procReads)

write.table(procReads, file.path(outputDir, sprintf("processedReadSummary%s.txt", postfix)), sep="\t", row.names=F, quote=F)
#procReads <- read.delim(file.path(outputDir, sprintf("processedReadSummary%s.txt", postfix)), stringsAsFactors = F)


##################################################################
## Test Sequence Quality
# reads <- list.files("data/", full.names = T)
# 
# reads <- c("data/Yakabu-2_S12_L001_R1_001.fastq.gz",
#            "data/Yakabu-2_S12_L001_R2_001.fastq.gz")
# qaSummary <- qa(reads, type="fastq")
# # browseURL(report(qaSummary))
# perCycle <- qaSummary[["perCycle"]]
# 
# pdf("plotQC.pdf")
# ShortRead:::.plotCycleQuality(perCycle$quality)
# dev.off()
