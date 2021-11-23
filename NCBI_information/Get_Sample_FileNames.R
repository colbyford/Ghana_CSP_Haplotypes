
library(dplyr)
library(stringr)
library(tidyr)

cwd <- "../CSP_sequence_analysis/data/"

## Get Sample File Information

samples <- list.files(cwd) %>% 
  data.frame() %>% 
  rename(FileName = ".") %>% 
  separate(FileName, into = c("SampleName", "meta"), sep = "-", extra = 'merge', remove = FALSE) %>% 
  separate(meta, into = c("gene", "rm1", "rm2", "direction"), sep = "_", extra = "drop", remove = TRUE) %>% 
  select(!starts_with("rm")) %>% 
  group_by(SampleName) %>% 
  pivot_wider(id_cols = "SampleName",
              names_from = c("gene", "direction"),
              values_from = "FileName")

write.csv(samples, "SampleFileNames.csv")
