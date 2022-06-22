library(dplyr)
library(readr)
library(tidyr)

# netmhc_output_path <- "netMHC_CSP_output_HLA-A0101.netmhcout"
netmhc_output_path <- "../netMHCpan-4.1/netMHC_CSP_output_CWIDalleles.fsa.netmhcout"


binders <- read_delim(netmhc_output_path,
                    delim = " ", comment = "#", skip = 50, trim_ws = TRUE,
                    skip_empty_rows = TRUE,
                    col_types = list(
                      "i",
                      "c",
                      "c",
                      "c",
                      "i",
                      "i",
                      "i",
                      "i",
                      "i",
                      "c",
                      "c",
                      "n",
                      "n",
                      "n",
                      "n",
                      "n",
                      "c",
                      "c"
                    ),
                    col_names = c(
                      "Pos",
                      "MHC",
                      "Peptide",
                      "Core",
                      "Of",
                      "Gp",
                      "Gl",
                      "Ip",
                      "Il",
                      "Icore",
                      "Identity",
                      "Score_EL",
                      "pct_Rank_EL",
                      "Score_BA",
                      "pct_Rank_BA",
                      "Aff_nM",
                      "Binder",
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
  )

write_csv(binders_output, paste0(netmhc_output_path, "_binders.csv"))


binders_pvt <- binders_output %>% 
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


write_csv(binders_pvt, paste0(netmhc_output_path, "_binders_pvt.csv"))
