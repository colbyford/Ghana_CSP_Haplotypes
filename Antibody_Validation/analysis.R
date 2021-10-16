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
  # ungroup(SampleID) %>% 
  # mutate(TotalReads = rowSums(select_if(., is.numeric), na.rm = TRUE))



gather(variable, value, -id) %>% 
  group_by(id) %>% 
  mutate(percentage = value/sum(value)) %>% 
  select(-value) %>% 
  spread(variable, percentage)


HADDOCK_data <- read_xlsx("../CSP_sequence_analysis/protein/HADDOCK_Results.xlsx")