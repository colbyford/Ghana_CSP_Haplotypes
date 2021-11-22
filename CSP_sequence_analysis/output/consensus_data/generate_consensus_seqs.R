library(Biostrings)
library(dplyr)
library(stringr)
library(tidyr)
library(DECIPHER)
library(ShortRead)

cwd <- "../dePlexMarker/"

## Get Sample File Information

deplex_samples <- list.files("../dePlexMarker/") %>% 
  data.frame() %>% 
  rename(FileName = ".") %>% 
  mutate(SampleName = str_extract(FileName, "[A-Za-z0-9]+-[0-9]_[A-Z0-9]+_[A-Za-z0-9]+_[a-z0-9]+")) %>% 
  separate(FileName, into = c("SampleID", "meta"), sep = "-", extra = 'merge', remove = FALSE) %>% 
  separate(meta, into = c("rm1", "rm2", "rm3", "gene", "direction"), sep = "_", extra = "drop", remove = TRUE) %>%
  mutate(direction = str_sub(direction,1,1)) %>% 
  select(!starts_with("rm")) %>% 
  mutate(size = file.size(paste0("../dePlexMarker/", FileName))) %>% 
  ## Get F and R files for each genes by picking files with most reads
  group_by(SampleID) %>%
  arrange(size, .by_group = TRUE) %>%
  top_n(4) %>% 
  select(!size)

csp_F_samples <- deplex_samples %>%
  filter(gene == "csp") %>%
  filter(direction == "F") %>% 
  select(!c(direction, gene)) %>% 
  rename("FileF_csp" = "FileName")

csp_R_samples <- deplex_samples %>%
  filter(gene == "csp") %>%
  filter(direction == "R") %>% 
  select(!c(direction, gene)) %>% 
  rename("FileR_csp" = "FileName")

csp_samples <- csp_F_samples %>% full_join(csp_R_samples)
csp_ref <- readFasta("../../meta_data/csp_reference.fasta") %>% sread()

msp1_F_samples <- deplex_samples %>%
  filter(gene == "msp1") %>%
  filter(direction == "F") %>% 
  select(!c(direction, gene)) %>% 
  rename("FileF_msp1" = "FileName")

msp1_R_samples <- deplex_samples %>%
  filter(gene == "msp1") %>%
  filter(direction == "R") %>% 
  select(!c(direction, gene)) %>% 
  rename("FileR_msp1" = "FileName")

msp1_samples <- msp1_F_samples %>% full_join(msp1_R_samples)
msp1_ref <- readFasta("../../meta_data/msp-1_reference.fasta") %>% sread()

## Loop Through Files to Generate Consensus FASTAs
## CSP
for (i in 1:nrow(csp_samples)){
  samplename <- csp_samples$SampleName[i]
  sampleid <- csp_samples$SampleID[i]
  
  cat("Getting CSP consensus sequence for sample:", samplename, "\n")

  csp_F <- readFastq(paste0(cwd, csp_samples$FileF_csp[i])) %>%
    sread() %>% 
    ConsensusSequence(ambiguity = FALSE)
  
  # csp_R <- readFastq(paste0(cwd, csp_samples$FileR_csp[i])) %>%
  #   sread() %>%
  #   reverseComplement() %>%
  #   ConsensusSequence(ambiguity = FALSE)

  csp_ali <- pairwiseAlignment(csp_ref, csp_F) %>% aligned()
  names(csp_ali) <- sampleid
  
  writeFasta(csp_ali, file = paste0(samplename,"_consensus.fasta"))
}


## MSP-1
for (i in 1:nrow(msp1_samples)){
  samplename <- msp1_samples$SampleName[i]
  sampleid <- msp1_samples$SampleID[i]
  
  cat("Getting MSP-1 consensus sequence for sample:", samplename, "\n")
  
  msp1_F <- readFastq(paste0(cwd, msp1_samples$FileF_msp1[i])) %>%
    sread() %>% 
    ConsensusSequence(ambiguity = FALSE)
  
  # msp1_R <- readFastq(paste0(cwd, msp1_samples$FileR_msp1[i])) %>%
  #   sread() %>%
  #   reverseComplement() %>%
  #   ConsensusSequence(ambiguity = FALSE)
  
  msp1_ali <- pairwiseAlignment(msp1_ref, msp1_F) %>% aligned()
  names(msp1_ali) <- sampleid
  
  writeFasta(msp1_ali, file = paste0(samplename,"_consensus.fasta"))

}

