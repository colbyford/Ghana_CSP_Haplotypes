library(Biostrings)
library(dplyr)
library(stringr)
library(tidyr)
library(DECIPHER)
library(ShortRead)
library(msa)

## Get Sample File Information

deplex_samples <- list.files("../dePlexMarker/") %>% 
  data.frame() %>% 
  rename(FileName = ".") %>% 
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

sample_files <- csp_samples %>% full_join(msp1_samples)

cwd <- "../dePlexMarker/"

## Loop Through Files to Generate Consensus FASTAs
for (i in 1:nrow(sample_files)){
  sampleid <- sample_files$SampleID[i]
  
  csp_F <- readFastq(paste0(cwd, sample_files$FileF_csp[i])) %>% sread() %>% msa(method = "Muscle", verbose = TRUE) %>% ConsensusSequence(ambiguity = FALSE)
  csp_R <- readFastq(paste0(cwd, sample_files$FileR_csp[i])) %>% sread() %>% reverseComplement() %>% ConsensusSequence(ambiguity = FALSE)
  
  csp <- append(csp_F, csp_R)
  
  csp_cons <- ConsensusSequence(csp, ambiguity = FALSE)
  
  csp_cons
}

