# Protein structural analysis of Ghanian _Plasmodium falciparum_ haplotypes with human immunoproteins explains variations in vaccine efficacy

[![DOI](https://zenodo.org/badge/DOI/10.1101/2021.12.29.474460.svg)](https://www.biorxiv.org/content/10.1101/2021.12.29.474460)

<p align="right">Cheikh Cambel Dieng, Colby T. Ford, Anita Lerch, Dickson Doniou, Jennifer Huynh, <br>Kovidh Vegesna, Jun-tao Guo, Daniel Janies, Linda Amoah, Yaw Afrane, and Eugenia Lo</p>


## Abstract

Today, the RTS,S malaria vaccine only provides partial protection (~36%) against _Plasmodium falciparum_ in certain sub-Saharan African countries such as Ghana, though this vaccine previously showed efficacy well above 60\% in previous clinical trials. Our study analyzes the diversity of the _PfCSP_ gene across 88 samples from Ghana and we have derived 27 prevalent haplotypes. We further investigated these haplotypes by understanding how variations in this gene affect protein binding interactions with human immunological proteins (a human leukocyte antigen and a T-cell receptor). By predicting and quantifying these interactions, we can then better understand variations in vaccine efficacy and posit sequence options for next-generation malaria vaccine design.

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
