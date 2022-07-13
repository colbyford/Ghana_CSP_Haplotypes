library(dplyr)
library(readr)
library(tidyr)

netmhcii_output_path <- "../netMHCIIpan-4.1/netMHCII_CSP_output_CWIDalleles.fsa.netmhcout"


binders <- read_delim(netmhcii_output_path,
                    delim = " ", comment = "#", skip = 15, trim_ws = TRUE,
                    skip_empty_rows = TRUE,
                    col_types = list(
                      "i",
                      "c",
                      "c",
                      "i",
                      "c",
                      "n",
                      "c",
                      "n",
                      "n",
                      "n",
                      "n",
                      "n",
                      "n",
                      "c"
                    ),
                    col_names = c(
                      "Pos",
                      "MHC",
                      "Peptide",
                      "Of",
                      "Core",
                      "Core_Rel",
                      "Identity",
                      "Score_EL",
                      "pct_Rank_EL",
                      "Exp_Bind",
                      "Score_BA",
                      "Aff_nM",
                      "pct_Rank_BA",
                      "BindLevel"
                    )
                    ) %>% 
  filter(!is.na(Pos))

unique(binders$Identity)


binders_output <- binders %>% 
  # filter(!is.na(BindLevel)) %>% 
  select(
    Pos,
    MHC,
    Peptide,
    Identity,
    Score_EL,
    pct_Rank_EL,
    Score_BA,
    pct_Rank_BA,
    Aff_nM,
    BindLevel
  ) %>% 
  mutate(BindLevel = case_when(
    BindLevel == "<=WB" ~ "WB",
    BindLevel == "<=SB" ~ "SB"
  ))

write_csv(binders_output, paste0(netmhcii_output_path, "_binders.csv"))


## PIVOT Data

binders_output <- read_csv(paste0(netmhcii_output_path, "_binders.csv"))

binders_pvt <- binders_output %>% 
  mutate(peptide_length  = nchar(Peptide),
         ensPos = Pos + peptide_length) %>% 
  ## filter to Th2R Region
  filter(Pos >= 29,
         ensPos <= 45) %>% 
  group_by(MHC, Identity) %>% 
  group_by(MHC, Identity) %>% 
  count(BindLevel) %>% 
  pivot_wider(id_cols = c(MHC, Identity),
              names_from = BindLevel,
              values_from = n,
              values_fill = 0) %>% 
  mutate(
    total = SB + WB + `NA`,
    SB_pct = SB/total,
    WB_pct = WB/total,
    Binder_pct = (SB + WB)/total
  )


write_csv(binders_pvt, paste0(netmhcii_output_path, "_binders_pvt.csv"))
