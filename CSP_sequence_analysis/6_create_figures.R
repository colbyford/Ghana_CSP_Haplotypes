# install.packages("ggmsa")

library(ggmsa)
library(ggplot2)

csp_aa <- Biostrings::readAAStringSet("output/csp_HaplotypeSeq_NGmergeSNP_filtered_translated.fasta")
# csp_aa <- system.file("extdata", "sample.fasta", package = "ggmsa")

# ggmsa(csp_aa[1], color = "Chemistry_AA", seq_name = TRUE)

## Using Chemistry AA Colors
ggmsa(csp_aa,
      # start = 30,
      # end = 83,
      # start = 1,
      # end = 111,
      start = 29,
      end = 95,
      # color = "Chemistry_AA",
      color = "Clustal",
      seq_name = TRUE,
      char_width = 0.5,
      posHighligthed = c(30:46, 60:83),
      disagreement = TRUE,
      font = "DroidSansMono") + 
  scale_x_continuous(breaks = c(29, 30, 46, 60, 83, 95),
                     # breaks = c(1, 30, 46, 60, 83, 111),
                     sec.axis = dup_axis(
                       breaks = c(38, 71),
                       labels = c("Th2R", "Th3R")
                     ))
 
  # scale_x_continuous(breaks = c(0, 30, 46, 60, 83, 110),
  #                    label = c("", "30\nTh2R", "46", "60\nTh3R", "83", "110"))

## Hack for Clustal Colors
p <- ggmsa(csp_aa,
           # start = 30,
           # end = 83,
           start = 1,
           end = 111,
           # color = "Chemistry_AA",
           color = "Clustal",
           seq_name = TRUE,
           char_width = 0.5,
           # posHighligthed = c(30:46, 60:83),
           disagreement = TRUE,
           font = "DroidSansMono") + 
  scale_x_continuous(breaks = c(1, 30, 46, 60, 83, 111),
                     sec.axis = dup_axis(
                       breaks = c(38, 71),
                       labels = c("Th2R", "Th3R")
                     ))

p$layers[[1]]$data$color <- ifelse(p$layers[[1]]$data$position %in% c(30:46, 60:83),
                                   p$layers[[1]]$data$color,
                                   "#ffffff")
p

## Tree
library(ggtree)
# csp_tree <- ape::read.tree("tree_building/RAxML_bipartitions.csp_HaplotypeSeq_NGmergeSNP_filtered_translated_matrix.out")
csp_tree <- ape::read.tree("tree_building/RAxML_bipartitions.csp_HaplotypeSeq_NGmergeSNP_filtered_translated_matrix - fixed_taxa.out")

metadata <- read.csv("tree_building/sample_metadata.csv", stringsAsFactors = FALSE, header = TRUE)

colors <- c(North='red',
            # Central='palegreen3',
            Central='olivedrab3',
            South='khaki')

ggtree(csp_tree, branch.length='none', layout='circular') %<+%
  metadata + 
  # geom_tiplab() +
  geom_tippoint(aes(color=Region), size = 3) +
  scale_color_manual(values=colors) + 
  geom_tiplab2(align=T, linetype=NA, size=3, offset=3, hjust=0.5)


