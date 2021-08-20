## Script to perform Kruskal-Wallis test on HADDOCK metrics
## for Peptides containing Y or Q

library(dplyr)
library(readr)
library(readxl)

data <- read_excel("HADDOCK_Results.xlsx", sheet = "Sheet1") %>% 
  filter(Region == "Th2R", Duplicate == FALSE) %>% 
  select(`CSP Peptide`,
         `HLA-TCR \r\nHADDOCK score`,
         `HLA-TCR \r\nVan der Waals energy`,
         `HLA-TCR \r\nElectrostatic energy`,
         `HLA-TCR \r\nDesolvation energy`,
         `HLA-TCR \r\nRestraints violation energy`,
         `Peptide Contains Tyrosine`,
         `Peptide Contains Glutamine`)

metrics <- c("`HLA-TCR \r\nHADDOCK score`",
             "`HLA-TCR \r\nVan der Waals energy`",
             "`HLA-TCR \r\nElectrostatic energy`",
             "`HLA-TCR \r\nDesolvation energy`",
             "`HLA-TCR \r\nRestraints violation energy`")

cont <- c("`Peptide Contains Tyrosine`",
          "`Peptide Contains Glutamine`")

tests <- data.frame()

for (metric in metrics){
  for (aa in cont){
    form <- paste0(metric, " ~ ", aa) %>% as.formula()
    
    kw <- kruskal.test(form, data = data)
    
    test_iter <- data.frame(`HADDOCK_Metric` = metric,
                            `Peptide_Contains` = aa,
                            chi_sq = kw$statistic,
                            p_value = kw$p.value)
    
    tests <- rbind(test_iter, tests)
  }
}

tests <- tests %>%
  mutate(HADDOCK_Metric = gsub("`", "", gsub("HLA-TCR \r\n","",
                               as.character(HADDOCK_Metric)))) %>% 
  mutate(Peptide_Contains = gsub("`", "", gsub("`Peptide Contains ","",
                                 as.character(Peptide_Contains)))) 

write_csv(tests, "kw_peptide_tests.csv")
