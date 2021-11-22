## Calculate Wilcoxon test for MOI for Region Comparison

library(readxl)
library(dplyr)

moi <- read_xlsx("MSP-CSP_mMOI_Comparison.xlsx", sheet="Pivot_raw") %>% select(msp1, Region)

wt <- pairwise.wilcox.test(moi$msp1, moi$Region,
                           p.adjust.method="bonferroni")

