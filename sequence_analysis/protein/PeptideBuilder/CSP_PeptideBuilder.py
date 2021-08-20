# -*- coding: utf-8 -*-
"""
Peptide Builder
for PfCSP Th2R/Th3R

Colby T. Ford, Ph.D.
"""

## Load in Libraries

import PeptideBuilder, Bio.PDB
from PeptideBuilder import Geometry
from Bio import SeqIO
import pandas as pd

## Input Sequences

# seq_str = "DIEIQN"
# seq = list(seq_str)

seq_records = SeqIO.parse("../../output/csp_HaplotypeSeq_NGmergeSNP_filtered_translated.fasta", "fasta")

cleave_data = pd.read_csv("NetChop3.0_20S_Thresh0.5.csv")

# peptide_res = [31,34,35,41,42,43]
# peptide_idx = [res - 1 for res in peptide_res]

## Th2R Peptide Builder
for seq_record in seq_records:
    seq_name = seq_record.id
    print(seq_name)
    
    seq_full = seq_record.seq
    seq_full = list(str(seq_full).replace('X', ''))
    
    peptide_range = [i for i in range(28, 44+1)] 
    
       
    peptide_res_df = cleave_data[(cleave_data.pos.isin(peptide_range)) &
                                 (cleave_data.C.isin(['S'])) &
                                 (cleave_data.Ident == seq_name)]
    
    peptide_res = peptide_res_df['pos'].tolist()
    
    peptide_idx = [res - 1 for res in peptide_res]
       
    
    seq = [seq_full[i] for i in peptide_idx]
    
    print(seq)
    
    ## Initialize Structure

    structure = PeptideBuilder.initialize_res(seq[0])
    seq.pop(0)
    
    ## Build Structure
    
    for res in seq:
        # print(res)
    
        geo = Geometry.geometry(res)
        structure = PeptideBuilder.add_residue(structure, geo)
    
    ## Add terminal oxygen (OXT) to the final glycine
    
    PeptideBuilder.add_terminal_OXT(structure)
    
    ## Output Structure
    
    output = Bio.PDB.PDBIO()
    output.set_structure(structure)
    
    print(seq_name + "_Th2R_peptide.pdb")
    print("\n")
    
    output.save(seq_name + "_Th2R_peptide.pdb")

## Th3R Peptide Builder
for seq_record in seq_records:
    seq_name = seq_record.id
    print(seq_name)
    
    seq_full = seq_record.seq
    seq_full = list(str(seq_full).replace('X', ''))
    
    peptide_range = [i for i in range(59, 81+1)] 
    
       
    peptide_res_df = cleave_data[(cleave_data.pos.isin(peptide_range)) &
                                 (cleave_data.C.isin(['S'])) &
                                 (cleave_data.Ident == seq_name)]
    
    peptide_res = peptide_res_df['pos'].tolist()
    
    peptide_idx = [res - 1 for res in peptide_res]
       
    
    seq = [seq_full[i] for i in peptide_idx]
    
    print(seq)
    
    ## Initialize Structure

    structure = PeptideBuilder.initialize_res(seq[0])
    seq.pop(0)
    
    ## Build Structure
    
    for res in seq:
        # print(res)
    
        geo = Geometry.geometry(res)
        structure = PeptideBuilder.add_residue(structure, geo)
    
    ## Add terminal oxygen (OXT) to the final glycine
    
    PeptideBuilder.add_terminal_OXT(structure)
    
    ## Output Structure
    
    output = Bio.PDB.PDBIO()
    output.set_structure(structure)
    
    print(seq_name + "_Th3R_peptide.pdb")
    print("\n")
    
    output.save(seq_name + "_Th3R_peptide.pdb")



