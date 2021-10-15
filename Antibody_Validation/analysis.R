library(dplyr)
library(readr)
library(readxl)
library(stringr)

antibody_data <- read_xls("csp antibody data.xls", sheet = "Samples All") %>%
  select(`Sample Name`, `Avg Conc`) %>% 
  filter(! `Sample Name` %in% c('no sample'))


CSP_data <- read_csv("../CSP_sequence_analysis/finalTab_csp.csv") %>% 
  select(SampleID, Haplotype, Reads) %>% 
  filter(! Haplotype %in% c("Chimera", "Indels", "Noise", "Singelton"))
