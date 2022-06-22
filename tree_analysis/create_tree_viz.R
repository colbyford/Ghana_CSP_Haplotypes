
# BiocManager::install("msa")
# library(msa)
# 
# csp_seqs <- readDNAStringSet("NCBI_PfCSP_sequences.fasta")
# 
# csp_alignment <- msa(csp_seqs)


library(ggtree)
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)

tree <- read.tree("RAxML_GUI_ModelTest_PfCSP_aligned_aa_wHaplotypes_modified.tree")
tree$tip.label <- replace(tree$tip.label, tree$tip.label=="csp_reference", "csp-reference")

metadata <- read_excel("plasmodiumAfrica.xlsx")

tip_labels <- data.frame(label = tree$tip.label) %>% 
  separate(label, c("accession"), sep = "_", extra = "drop") %>% 
  mutate(accession = ifelse(str_detect(accession, "^csp-"), accession, paste0(accession,".1")))
  
tree$tip.label <- tip_labels$accession


tip_meta <- tip_labels %>%
  left_join(metadata, by = c("accession" = "Record ID")) %>%
  replace_na(list(Country = "HAPLOTYPE", Match = "Haplotype")) %>% 
  mutate(hap_accession = ifelse(Match == "Haplotype", accession, ""))


missing_meta_tips <- tip_meta %>% filter(Match == "Missing") %>% select(accession)

tree_reduced <- tree %>% treeio::drop.tip(missing_meta_tips$accession) %>% treeio::root(outgroup = "csp-reference", edgelabel = TRUE)

# only_hap_tips <- tip_meta %>%
#   filter(accession %in% tree_reduced$tip.label) %>% 
#   mutate(accession = ifelse(Match == "Haplotype", accession, ""))
# 
# tree_reduced$tip.label <- only_hap_tips$accession


ggtree(tree_reduced, layout="rectangular", ladderize = FALSE, branch.length = "none") %<+% 
  tip_meta +
  # geom_text(aes(label=node), hjust=.3) +
  geom_tippoint(aes(color=Match)) +
  geom_tiplab(aes(label=hap_accession), hjust = -.1)
  # geom_tiplab(hjust = -.1)


p1 <- ggtree(tree_reduced, layout="circular", ladderize = TRUE, branch.length = "none") %<+% 
  tip_meta +
  # geom_text(aes(label=node), hjust=.3) +
  geom_tippoint(aes(color=Match)) +
  geom_tiplab(aes(label=hap_accession), hjust = -.1) +
  scale_color_manual(values=c("chocolate2", "purple4", "green4"))
  

p2 <- p1 %>% 
  ggtree::collapse(node=743) %>%
  ggtree::collapse(node=1931) %>% 
  ggtree::collapse(node=2215) %>% 
  ggtree::collapse(node=2453) + 
  geom_point2(aes(subset=(node==743)), shape=23, size=2, fill='grey') +
  geom_point2(aes(subset=(node==1931)), shape=23, size=2, fill='grey') +
  geom_point2(aes(subset=(node==2215)), shape=23, size=2, fill='grey') +
  geom_point2(aes(subset=(node==2413)), shape=23, size=2, fill='grey')

matches <- list(
  Elsewhere = p1$data$label[which(p1$data$Match == "Elsewhere")],
  Africa = p1$data$label[which(p1$data$Match == "Africa")],
  Haplotype = p1$data$label[which(p1$data$Match == "Haplotype")]
)

p3 <- groupOTU(p2, matches, 'match') +
  aes(color=match)

p3


ggplotly(p1)
