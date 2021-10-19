library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(stringr)

antibody_data <- read_xls("csp antibody data.xls", sheet = "Samples All") %>%
  filter(! `Sample Name` %in% c('no sample')) %>% 
  mutate(SampleID = tolower(`Sample Name`)) %>% 
  select(SampleID, `Avg Conc`)


## DESeq2 Prep
library(DESeq2)
CSP_data_orig <- read_csv("../CSP_sequence_analysis/finalTab_csp.csv") %>% 
  select(SampleID, Haplotype, Reads) %>% 
  filter(! Haplotype %in% c("Chimera", "Indels", "Noise", "Singelton")) %>% ## Only True Haplotypes
  # mutate(Reads = Reads / exp(mean(log(Reads[Reads > 0])))) %>%  ## Geometric Mean Ratio 
  pivot_wider(names_from = "Haplotype", values_from = "Reads", values_fill = 1) %>% 
  mutate(SampleID = str_split(SampleID, "-", n = 2) %>% sapply(head, 1)) %>% 
  filter(SampleID %in% antibody_data$SampleID)

CSP_readcounts <- CSP_data_orig %>% 
  select(-SampleID) %>% 
  t() %>% data.frame()

colnames(CSP_readcounts) <- CSP_data_orig$SampleID

CSP_sampleinfo <- data.frame(row.names = names(CSP_readcounts)) %>% 
  mutate(gene = "CSP")


DESeq.ds <- DESeqDataSetFromMatrix(countData = CSP_readcounts,
                                   colData   = CSP_sampleinfo,
                                   design    = ~ 1)

colSums(counts(DESeq.ds)) %>% barplot

DESeq.ds <- estimateSizeFactors(DESeq.ds) # calculate SFs, add them to object
plot( sizeFactors(DESeq.ds), colSums(counts(DESeq.ds)), # assess them
      ylab = "library sizes", xlab = "size factors", cex = .6 )


log.norm.counts <- log2(counts(DESeq.ds, normalized = TRUE))

##########

reads_data <- read_tsv("../CSP_sequence_analysis/output/demultiplexMarkerSummary.txt") %>% 
  filter(MarkerID == "csp") %>% 
  select(SampleID, numReadOut) %>% 
  mutate(SampleID = str_split(SampleID, "-", n = 2) %>% sapply(head, 1))



# CSP_data <- read_csv("../CSP_sequence_analysis/finalTab_csp.csv") %>% 
#   select(SampleID, Haplotype, Reads) %>% 
#   filter(! Haplotype %in% c("Chimera", "Indels", "Noise", "Singelton")) %>% ## Only True Haplotypes
#   mutate(SampleID = str_split(SampleID, "-", n = 2) %>% sapply(head, 1)) %>% 
#   filter(SampleID %in% antibody_data$SampleID) %>% 
#   inner_join(reads_data, by = "SampleID") %>% 
#   mutate(expression = Reads/numReadOut) %>% 
#   select(SampleID, Haplotype, expression)

 
CSP_data <- read_csv("../CSP_sequence_analysis/finalTab_csp.csv") %>%
  select(SampleID, Haplotype, Reads) %>%
  filter(! Haplotype %in% c("Chimera", "Indels", "Noise", "Singelton")) %>% ## Only True Haplotypes
  # mutate(Reads = Reads / exp(mean(log(Reads[Reads > 0])))) %>%  ## Geometric Mean Ratio
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
  pivot_longer(cols = starts_with("csp-"), names_to = "Haplotype", values_to = "expression") %>%
  filter(expression > 0)



HADDOCK_data <- read_xlsx("../CSP_sequence_analysis/protein/HADDOCK_Results.xlsx") %>% 
  filter(Region == "Th2R") %>% 
  select(`CSP Peptide`, `HLA-CSPpep\r\nHADDOCK score`, `HLA-TCR \r\nHADDOCK score`) %>%
  rename(Haplotype = `CSP Peptide`) %>% 
  rename(HLA_CSP_HADDOCK = `HLA-CSPpep\r\nHADDOCK score`) %>% 
  rename(HLA_TCR_HADDOCK = `HLA-TCR \r\nHADDOCK score`)


######

comparison_data <- CSP_data %>%
  inner_join(antibody_data, by = "SampleID") %>% 
  mutate(weighted_avg_conc = expression * as.numeric(`Avg Conc`)) %>% 
  inner_join(HADDOCK_data, by = "Haplotype") %>%
  mutate(weighted_HLA_CSP_HADDOCK = expression * HLA_CSP_HADDOCK,
         weighted_HLA_TCR_HADDOCK = expression * HLA_TCR_HADDOCK) %>% 
  group_by(SampleID) %>%
  mutate(weighted_avg_HLA_CSP_HADDOCK = sum(weighted_HLA_CSP_HADDOCK),
         weighted_avg_HLA_TCR_HADDOCK = sum(weighted_HLA_TCR_HADDOCK))  
  # filter(expression == max(expression))

write_csv(comparison_data, "antibody_conc_HADDOCK_comparison.csv")

# 
# HLA_CSP_model <- lm(HLA_CSP_HADDOCK ~ weighted_avg_conc, data = comparison_data)
# HLA_TCR_model <- lm(HLA_TCR_HADDOCK ~ weighted_avg_conc, data = comparison_data)
