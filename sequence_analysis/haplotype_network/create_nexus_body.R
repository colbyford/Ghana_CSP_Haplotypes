## Script to prepare the body sequences of an
## input Nexus file for TCs haplotype network generation

library(dplyr)
library(stringr)
library(readr)
library(ape)


sample_regions <- read_csv("haplotype_regions.csv")


## CSP

csp_haps <- readRDS("../finalTab_csp.RDS")[['csp']] %>% 
  select(SampleName, Haplotype) %>% 
  filter(str_starts(Haplotype, "csp-"))
  
# csp_fasta <- "../output/csp_HaplotypeSeq_NGmergeSNP_filtered.fasta"
csp_fasta <- "../output/csp_HaplotypeSeq_NGmergeSNP_filtered_translated.fasta"

# csp_seqs <- read.FASTA(csp_fasta) %>%
csp_seqs <- read.FASTA(csp_fasta, type = "AA") %>%
  as.character() %>% 
  lapply(paste0,collapse="") %>% 
  unlist() %>% 
  data.frame(stringsAsFactors = FALSE) %>% 
  tibble::rownames_to_column() %>% 
  rename("seq"=".") %>% 
  mutate(seq = toupper(seq)) %>% 
  rename("Haplotype" = "rowname")
  

csp_body <- csp_haps %>%
  left_join(csp_seqs, by = "Haplotype") %>% 
  mutate(name = paste0(Haplotype, "_", SampleName)) %>% 
  select(name, seq)

### Write the NEXUS body sequences
# write_delim(csp_body, "csp_nexus_body.txt", delim = "\t", col_names = FALSE)
write_delim(csp_body, "csp_aminoacids_nexus_body.txt", delim = "\t", col_names = FALSE)


### Write the Haplotype group list for tcsBU
csp_groups <- csp_haps %>% 
  left_join(sample_regions, by = "SampleName") %>% 
  mutate(name = paste0(Haplotype, "_", SampleName)) %>% 
  select(name, Region)

write_delim(csp_groups, "csp_groups.csv", delim = ";", col_names = FALSE)



## MSP-1

msp_haps <- readRDS("../../MSP/finalTab_msp.RDS")[['msp1']]%>% 
  select(SampleName, Haplotype) %>% 
  filter(str_starts(Haplotype, "msp1-"))


msp_seqs <- read.FASTA("../../MSP/output/msp1_HaplotypeSeq_bind180_180SNP_filtered_aligned.fasta") %>%
  as.character() %>% 
  lapply(paste0,collapse="") %>% 
  unlist() %>% 
  data.frame(stringsAsFactors = FALSE) %>% 
  tibble::rownames_to_column() %>% 
  rename("seq"=".") %>% 
  mutate(seq = toupper(seq)) %>% 
  rename("Haplotype" = "rowname")

msp_body <- msp_haps %>%
  left_join(msp_seqs, by = "Haplotype") %>% 
  mutate(name = paste0(Haplotype, "_", SampleName)) %>% 
  select(name, seq)

### Write the NEXUS body sequences
write_delim(msp_body, "msp_nexus_body.txt", delim = "\t", col_names = FALSE)

### Write the Haplotype group list for tcsBU
msp_groups <- msp_haps %>% 
  left_join(sample_regions, by = "SampleName") %>% 
  mutate(name = paste0(Haplotype, "_", SampleName)) %>% 
  select(name, Region)

write_delim(msp_groups, "msp_groups.csv", delim = ";", col_names = FALSE)

