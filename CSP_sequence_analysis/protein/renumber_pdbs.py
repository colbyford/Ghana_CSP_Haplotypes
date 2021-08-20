# -*- coding: utf-8 -*-
"""
Renumber PDB structures of the HLA-CSPpep Complex using PyMol

@author: Colby
"""
import glob, os, re
from pymol import cmd

for file in glob.glob("./CD*_Th*R/pdb_files/HLA-CSPpep_csp-*/cluster*.pdb"):

    output_file = file.replace(".pdb", "_renumbered.pdb")
    print(output_file)

    ## PyMol Part
    cmd.load(file)

    cmd.alter('chain B', 'resi=int(resi)+500')

    cmd.save(output_file, selection='all', format='pdb')

    cmd.reinitialize('everything')
