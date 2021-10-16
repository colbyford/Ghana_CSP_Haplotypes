library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(stringr)

antibody_data <- read_xls("csp antibody data.xls", sheet = "Samples All") %>%
  filter(! `Sample Name` %in% c('no sample')) %>% 
  mutate(SampleID = tolower(`Sample Name`)) %>% 
  select(SampleID, `Avg Conc`)


CSP_data <- read_csv("../CSP_sequence_analysis/finalTab_csp.csv") %>% 
  select(SampleID, Haplotype, Reads) %>% 
  filter(! Haplotype %in% c("Chimera", "Indels", "Noise", "Singelton")) %>% 
  pivot_wider(names_from = "Haplotype", values_from = "Reads", values_fill = 0) %>% 
  mutate(SampleID = str_split(SampleID, "-", n = 2) %>% sapply(head, 1)) %>% 
  filter(SampleID %in% antibody_data$SampleID) %>% 
  # mutate(TotalReads = rowSums(select_if(., is.numeric), na.rm = TRUE))
  gather(variable, value, -SampleID) %>% 
  group_by(SampleID) %>% 
  mutate(percentage = value/sum(value)) %>% 
  select(-value) %>% 
  spread(variable, percentage) %>% 
  ungroup(SampleID) %>% 
  # mutate(TotalReads = rowSums(select_if(., is.numeric), na.rm = TRUE))
  pivot_longer(cols = starts_with("csp-"), names_to = "Haplotype", values_to = "expression")



HADDOCK_data <- read_xlsx("../CSP_sequence_analysis/protein/HADDOCK_Results.xlsx") %>% 
  filter(Region == "Th2R") %>% 
  select(`CSP Peptide`, `HLA-CSPpep\r\nHADDOCK score`, `HLA-TCR \r\nHADDOCK score`) %>%
  rename(Haplotype = `CSP Peptide`) %>% 
  rename(HLA_CSP_HADDOCK = `HLA-CSPpep\r\nHADDOCK score`) %>% 
  rename(HLA_TCR_HADDOCK = `HLA-TCR \r\nHADDOCK score`)



comparison_data <- CSP_data %>%
  inner_join(antibody_data, by = "SampleID") %>% 
  mutate(weighted_avg_conc = expression * as.numeric(`Avg Conc`)) %>% 
  inner_join(HADDOCK_data, by = "Haplotype") %>% 
  filter(expression > 0)
