library(ape)
library(phangorn)
library(msa)
library(dplyr)
library(tidyr)
library(stringr)

outputDir = "output/"

## [6] Filter Output Fasta Files

###########################
### Filter MSP-1 SNP
msp_snp_fasta_path <- list.files(path = outputDir,
                                 pattern = "(msp).+SNP\\.fasta$")

msp_ref_fasta_path <- "meta_data/msp-1_reference.fasta"

msp_SNP_haplotypes <- readRDS("finalTab_msp.RDS")[["msp1"]][1:5] %>% 
  filter(str_detect(Haplotype, "msp1"))

msp_fasta <- ape::read.FASTA(paste0(outputDir, "/", msp_snp_fasta_path), type = "DNA")
msp_ref_fasta <- ape::read.FASTA(msp_ref_fasta_path, type = "DNA")

#### Get list of Final Haplotypes
msp_haplotypes_filtered <- unique(msp_SNP_haplotypes$Haplotype)

#### Remove Haplotypes and Append Reference Sequence
msp_fasta_filtered <- c(msp_ref_fasta, msp_fasta[msp_haplotypes_filtered])

#### Write out the filtered FASTA File
ape::write.FASTA(msp_fasta_filtered, paste0(outputDir, "/", str_remove(msp_snp_fasta_path, ".fasta"),"_filtered.fasta"))

#### Align Seqs
write_alignment <- function(alignment, filename) {
  sink(filename)
  for(i in seq(1, length(rownames(alignment)))) {
    cat(paste0('>', rownames(alignment)[i]))
    cat('\n')
    seq <- toString(unmasked(alignment)[[i]])
    cat(seq)
    cat('\n')  
  }
  sink(NULL)
}

#### Perform Alignment
msp_to_align <- readDNAStringSet(paste0(outputDir, "/", str_remove(msp_snp_fasta_path, ".fasta"), "_filtered.fasta"))
msp_alignment <- msa(msp_to_align)

write_alignment(msp_alignment, paste0(outputDir, "/", str_remove(msp_snp_fasta_path, ".fasta"), "_filtered_aligned.fasta"))


## [7] Create Haplotype Trees

### MSP
msp_aligned_fasta <- ape::read.FASTA(paste0(outputDir, "/", str_remove(msp_snp_fasta_path, ".fasta"), "_filtered_aligned.fasta"), type = "DNA")
msp_dist <- dist.ml(msp_aligned_fasta, model="JC69")

# msp_UPGMA <- upgma(msp_dist) %>% root(outgroup = "msp-1_reference")
# plot(msp_UPGMA, main="UPGMA")
msp_NJ  <- NJ(msp_dist) %>% root(outgroup = "msp-1_reference")
# plot(msp_NJ, main="NJ")

# msp_optim <- pml(msp_NJ, data = phyDat(msp_aligned_fasta))
msp_optim <- optim.parsimony(msp_NJ, phyDat(msp_aligned_fasta))
# plot(msp_optim)
write.tree(msp_optim, file = paste0(outputDir, "/", str_remove(msp_snp_fasta_path, ".fasta"), "_filtered_aligned_parsimony.tre"))


