# -*- coding: utf-8 -*-
"""
Generate PyMol Rendering of the HLA-CSPpep-TCR complex

@author: Colby
"""
import glob, os, re
from pymol import cmd

for file in glob.glob("../interaction*/HLA-TCR/cluster*.pdb"):

    output_file = re.search('interaction_csp-[a-z0-9]+', file).group(0) + ".png"
    hap_num = re.search('-[a-z0-9]+', output_file).group(0).replace('-','')
    
    # print(hap_num)
    
    ## PyMol Part
    cmd.load(file)
    
    ## Color TCR green and HLA cyan
    cmd.color('green', 'chain B')
    cmd.color('cyan', 'chain A')
    
    ## Color CSP peptide magenta and turn into spheres
    cmd.select('csppep', 'resi 501-511')
    cmd.color('magenta', 'csppep')
    cmd.show('spheres', 'csppep')
    
    
    ## Rotate structures into similar orientation
    if hap_num in ['12', '19','22','32','35','68','94','117']:
        cmd.rotate('y', 180)
        
    if hap_num in ['7', '8', '15', '47', '117']:
        cmd.rotate('x', -90)
        
    if hap_num in ['16', '35', '68', '94']:
        cmd.rotate('x', 90)
        
    if hap_num in ['reference', '10', '12', '19']:
        cmd.rotate('x', 180)
        
    if hap_num in ['reference', '1', '19', '35', '94']:
        cmd.rotate('x', -30)
        
    if hap_num in ['15', '117']:
        cmd.rotate('x', 15)
        
    if hap_num in ['12', '15', '117']:
        cmd.rotate('x', 45)
        cmd.rotate('y', -20)
    
    cmd.rotate('y', 20)
    
    cmd.set('opaque_background', 0)
    
    cmd.png(output_file, 0, 0, -1, ray=1)
    
    cmd.reinitialize('everything')
