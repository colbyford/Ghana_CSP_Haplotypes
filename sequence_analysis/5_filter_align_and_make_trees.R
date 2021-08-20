library(ape)
library(phangorn)
library(msa)
library(dplyr)
library(tidyr)
library(stringr)

outputDir = "output/"

## [6] Filter Output Fasta Files

###########################
### MSP Filter MSP-1 LP and Align
msp_lp_fasta_path <- list.files(path = outputDir,
                                pattern = "(msp1).+LP\\.fasta$")

msp_ref_fasta_path <- "meta_data/msp-1_reference.fasta"

msp_LP_haplotypes <- readRDS("finalTab_LP.RDS")[["msp1"]][1:5] %>% 
  filter(str_detect(Haplotype, "msp1"))

msp_fasta <- ape::read.FASTA(paste0(outputDir, "/", msp_lp_fasta_path), type = "DNA")
msp_ref_fasta <- ape::read.FASTA(msp_ref_fasta_path, type = "DNA")

#### Get list of Final Haplotypes
msp_haplotypes_filtered <- unique(msp_LP_haplotypes$Haplotype)

#### Remove Haplotypes and Append Reference Sequence
msp_fasta_filtered <- c(msp_ref_fasta, msp_fasta[msp_haplotypes_filtered])

#### Write out the filtered FASTA File
ape::write.FASTA(msp_fasta_filtered, paste0(outputDir, "/", str_remove(msp_lp_fasta_path, ".fasta"), "_filtered.fasta"))

#### Align LPs
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
msp_to_align <- readDNAStringSet(paste0(outputDir, "/", str_remove(msp_lp_fasta_path, ".fasta"), "_filtered.fasta"))
msp_alignment <- msa(msp_to_align)

write_alignment(msp_alignment, paste0(outputDir, "/", str_remove(msp_lp_fasta_path, ".fasta"), "_filtered_aligned.fasta"))


###########################
### Filter CSP SNP
csp_snp_fasta_path <- list.files(path = outputDir,
                                 pattern = "(csp).+SNP\\.fasta$")

csp_ref_fasta_path <- "meta_data/csp_reference.fasta"

csp_SNP_haplotypes <- readRDS("finalTab_csp.RDS")[["csp"]][1:5] %>% 
  filter(str_detect(Haplotype, "csp"))

csp_fasta <- ape::read.FASTA(paste0(outputDir, "/", csp_snp_fasta_path), type = "DNA")
csp_ref_fasta <- ape::read.FASTA(csp_ref_fasta_path, type = "DNA")

#### Get list of Final Haplotypes
csp_haplotypes_filtered <- unique(csp_SNP_haplotypes$Haplotype)

#### Remove Haplotypes and Append Reference Sequence
csp_fasta_filtered <- c(csp_ref_fasta, csp_fasta[csp_haplotypes_filtered])

#### Write out the filtered FASTA File
ape::write.FASTA(csp_fasta_filtered, paste0(outputDir, "/", str_remove(csp_snp_fasta_path, ".fasta"),"_filtered.fasta"))


## [7] Create Haplotype Trees

### MSP
msp_aligned_fasta <- ape::read.FASTA(paste0(outputDir, "/", str_remove(msp_lp_fasta_path, ".fasta"), "_filtered_aligned.fasta"), type = "DNA")
msp_dist <- dist.ml(msp_aligned_fasta, model="JC69")

# msp_UPGMA <- upgma(msp_dist) %>% root(outgroup = "msp-1_reference")
# plot(msp_UPGMA, main="UPGMA")
msp_NJ  <- NJ(msp_dist) %>% root(outgroup = "msp-1_reference")
# plot(msp_NJ, main="NJ")

# msp_optim <- pml(msp_NJ, data = phyDat(msp_aligned_fasta))
msp_optim <- optim.parsimony(msp_NJ, phyDat(msp_aligned_fasta))
# plot(msp_optim)
write.tree(msp_optim, file = paste0(outputDir, "/", str_remove(msp_lp_fasta_path, ".fasta"), "_filtered_aligned_parsimony.tre"))

### CSP
csp_dist <- dist.ml(csp_fasta_filtered, model="JC69")
csp_NJ  <- NJ(csp_dist) %>% root(outgroup = "csp_reference")
# plot(csp_NJ, main="NJ")

csp_optim <- optim.parsimony(csp_NJ, phyDat(csp_fasta_filtered))
plot(csp_optim)
write.tree(csp_optim, file = paste0(outputDir, "/", str_remove(csp_snp_fasta_path, ".fasta"), "_filtered_parsimony.tre"))



############
## Create Matrix for Sample-based Tree Generation

# csp_fasta_filtered_aa <- ape::read.FASTA(paste0(outputDir, "/csp_HaplotypeSeq_NGmergeSNP_filtered_translated.fasta"), type = "AA")
# 
# csp_fasta_filtered_aa_df <- ape::as.character.AAbin(csp_fasta_filtered_aa) %>%
#   sapply(paste, collapse="") %>%
#   as.data.frame()
# 
# colnames(csp_fasta_filtered_aa_df) <- c("aa_sequence")
# csp_fasta_filtered_aa_df$Haplotype <- row.names(csp_fasta_filtered_aa_df)
# 
# 
# csp_sample_haplotypes <- read.csv("finalTab_csp.csv") %>%
#   filter(!Haplotype %in% c("Noise", "Indels", "Singelton", "Chimera")) %>% 
#   select(SampleName, Haplotype) %>% 
#   left_join(csp_fasta_filtered_aa_df)
# 
# csp_sample_haplotypes_pvt <- csp_sample_haplotypes %>%
#   pivot_wider(names_from = Haplotype,
#               values_from = aa_sequence,
#               values_fill = rep("-", nchar(csp_sample_haplotypes$aa_sequence[1])) %>% paste0(collapse = ""))
# 
# write.table(csp_sample_haplotypes_pvt,
#             paste0(outputDir, "/csp_HaplotypeSeq_NGmergeSNP_filtered_translated_matrix.txt"),
#             sep = "\t",
#             quote = FALSE,
#             row.names = FALSE)

