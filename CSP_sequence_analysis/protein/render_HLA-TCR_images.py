# -*- coding: utf-8 -*-
"""
Render HLA-TCR PDB images using PyMol

@author: Colby
"""
import glob, os, re
from pymol import cmd

##################################
## Th2R Region - CD4+

## Get HLA reference structure (for alignment)
hla_ref_path = "./CD4_Th2R/pdb_files/HLA-DRB1_Class2/6v1a_onlyMHC_renumbered.pdb"
hla_ref_name = "6v1a_onlyMHC_renumbered"

for file in glob.glob("./CD4_Th2R/pdb_files/HLA-TCR_csp-*/cluster*.pdb"):

    output_file = "./CD4_Th2R/renderings/" + re.search('HLA-TCR_csp-[a-z0-9]+', file).group(0) + ".png"

    hap_num = re.search('-[a-z0-9]+', output_file).group(0).replace('-','')
    
    obj_name = re.search('cluster[a-z0-9\_]+', file).group(0)

    print(f"Number: {hap_num} \tIn: {file} (called {obj_name}) \tOut: {output_file}")
    
    ## PyMol Part

    ## Load HLA reference structure
    cmd.load(hla_ref_path)

    ## Load HLA-TCR structure
    cmd.load(file)

    ## Align HLA-TCR structure to HLA reference structure
    cmd.align(obj_name, hla_ref_name)

    ## Remove HLA reference structure
    cmd.delete(hla_ref_name)

    ## Color TCR green and HLA cyan
    cmd.color('green', 'chain B')
    cmd.color('cyan', 'chain A')
    
    ## Color CSP peptide magenta and turn into spheres
    cmd.select('csppep', 'resi 501-511')
    cmd.color('magenta', 'csppep')
    cmd.show('spheres', 'csppep')

    ## Rotate structure into similar orientation
    cmd.rotate('y', -105)
    cmd.rotate('x', 105)
    cmd.rotate('z', -10)

    cmd.zoom(obj_name, complete = 1)
    
    cmd.set('opaque_background', 0)
    
    cmd.png(output_file, width = 3640, height = 2160, dpi = 600, ray=1)
    
    cmd.reinitialize('everything')



##################################
## Th3R Region - CD8+

## Get HLA reference structure (for alignment)
hla_ref_path = "./CD8_Th3R/pdb_files/HLA-A_Class1/6trn_onlyMHC_renumbered.pdb"
hla_ref_name = "6trn_onlyMHC_renumbered"

for file in glob.glob("./CD8_Th3R/pdb_files/HLA-TCR_csp-*/cluster*.pdb"):

    output_file = "./CD8_Th3R/renderings/" + re.search('HLA-TCR_csp-[a-z0-9]+', file).group(0) + ".png"

    hap_num = re.search('-[a-z0-9]+', output_file).group(0).replace('-','')
    
    obj_name = re.search('cluster[a-z0-9\_]+', file).group(0)

    print(f"Number: {hap_num} \tIn: {file} (called {obj_name}) \tOut: {output_file}")
    
    ## PyMol Part

    ## Load HLA reference structure
    cmd.load(hla_ref_path)

    ## Load HLA-TCR structure
    cmd.load(file)

    ## Align HLA-TCR structure to HLA reference structure
    cmd.align(obj_name, hla_ref_name)

    ## Remove HLA reference structure
    cmd.delete(hla_ref_name)

    ## Color TCR green and HLA cyan
    cmd.color('green', 'chain B')
    cmd.color('cyan', 'chain A')
    
    ## Color CSP peptide magenta and turn into spheres
    cmd.select('csppep', 'resi 501-511')
    cmd.color('magenta', 'csppep')
    cmd.show('spheres', 'csppep')

    ## Rotate structure into similar orientation
    cmd.rotate('z', 30)

    cmd.zoom(obj_name, complete = 1)
    
    cmd.set('opaque_background', 0)
    
    cmd.png(output_file, width = 3640, height = 2160, dpi = 600, ray=1)
    
    cmd.reinitialize('everything')