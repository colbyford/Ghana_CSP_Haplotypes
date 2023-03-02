# Genetic variations of _Plasmodium falciparum_ circumsporozoite protein and the impact on interactions with human immunoproteins and malaria vaccine efficacy

<p align="right">Cheikh Cambel Dieng, Colby T. Ford, Anita Lerch, Dickson Doniou, Jennifer Huynh, Kovidh Vegesna, <br>Jun-tao Guo, Daniel Janies, Liwang Cui, Linda Amoah, Yaw Afrane, and Eugenia Lo</p>

Published in _Infection, Genetics and Evolution_: [![DOI](https://zenodo.org/badge/DOI/10.1016/j.meegid.2023.105418.svg)](https://doi.org/10.1016/j.meegid.2023.105418)


## Abstract

In October 2021, the world's first malaria vaccine RTS,S was endorsed by WHO for broad use in children, despite its low efficacy. This study examined polyclonal infections and the associations of parasite genetic variations with binding affinity to human leukocyte antigen (HLA). Multiplicity of infection was determined by amplicon deep sequencing of _PfMSP1_. Genetic variations in _PfCSP_ were examined across 88 samples from Ghana and analyzed together with 1655 _PfCSP_ sequences from other African and non-African isolates. Binding interactions of _PfCSP_ peptide variants and HLA were predicted using NetChop and HADDOCK. High polyclonality was detected among infections, with each infection harboring multiple non-3D7 _PfCSP_ variants. Twenty-seven _PfCSP_ haplotypes were detected in the Ghanaian samples, and they broadly represented _PfCSP_ diversity across Africa. The number of genetic differences between 3D7 and non-3D7 _PfCSP_ variants does not influence binding to HLA. However, CSP peptide length after proteolytic degradation significantly affects its molecular weight and binding affinity to HLA. Despite the high diversity of HLA, the majority of the HLAI and II alleles interacted/bound with all Ghana CSP peptides. Multiple non-3D7 strains among _P. falciparum_ infections could impact the effectiveness of RTS,S. Longer peptides of the Th2R/Th3R CSP regions should be considered in future versions of RTS,S.

<p align="middle"><img src="https://raw.githubusercontent.com/colbyford/Ghana_CSP_Haplotypes/main/CSP_sequence_analysis/protein/CD4_Th2R/renderings/HLA-TCR_csp-reference.png" width="500"></p>


## Process

<p align="middle"><img src="https://raw.githubusercontent.com/colbyford/Ghana_CSP_Haplotypes/main/CSP_sequence_analysis/figures/process_v2.png" width="500"></p>

## Summarized Table of Contents

- [CSP_sequence_analysis/](CSP_sequence_analysis/)
  - Final List of CSP Haplotypes (SNP-based): [finalTab_csp.csv](/CSP_sequence_analysis/finalTab_csp.csv)
  - Haplotype Creation Script: [MalariaHaplotyping_NGmerge_script.R](/CSP_sequence_analysis/MalariaHaplotyping_NGmerge_script.R)
  - Haplotype Network Resources: [haplotype_network/](/CSP_sequence_analysis/haplotype_network/)
  - Protein Modeling Resources: [protein/](/CSP_sequence_analysis/protein/)
    - HADDOCK Docking Results: [HADDOCK_Results.xlsx](/CSP_sequence_analysis/protein/HADDOCK_Results.xlsx)
    - CD4+/Th2R Protein Docking Files: [CD4_Th2R/](/CSP_sequence_analysis/protein/CD4_Th2R/)
    - CD8+/Th3R Protein Docking Files: [CD8_Th3R/](/CSP_sequence_analysis/protein/CD8_Th23/)
- [MSP1_sequence_analysis/](MSP1_sequence_analysis/)
  - Final List of MSP1 Haplotypes: [finalTab_msp.csv](/MSP1_sequence_analysis/finalTab_msp.csv)
  - Haplotype Creation Script: [MalariaHaplotyping_NGmerge_script.R](/MSP1_sequence_analysis/MalariaHaplotyping_NGmerge_script.R)
- [CSP-MSP1_Comparison/](CSP-MSP1_Comparison/)
  -  Multiplicity of Infection by Sample: [MSP-CSP_mMOI_Comparison.xlsx](/CSP-MSP1_Comparison/MSP-CSP_mMOI_Comparison.xlsx)

## Data
All Amplicon-seq paired-end FASTQ files of CSP and MSP-1 genes are available on NCBI SRA: https://www.ncbi.nlm.nih.gov/sra/PRJNA783000

## Resources

- HaplotypR pacakge: https://github.com/lerch-a/HaplotypR
- NetChop 3.1 Server: http://www.cbs.dtu.dk/services/NetChop/
- PeptideBuilder package: https://github.com/clauswilke/peptidebuilder
- HADDOCK Server: https://wenmr.science.uu.nl/haddock2.4/
- PyMol: https://pymol.org/
- PDB Structures:
  - HLA-DR: [`6V1A`](https://www.rcsb.org/structure/6V1A)
  - HLA-A2*01: [`6TRN`](https://www.rcsb.org/structure/6TRN)
  - TCR T594: [`6PY2`](https://www.rcsb.org/structure/6py2)
