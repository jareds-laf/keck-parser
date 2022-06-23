# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 11:50:59 2022

@author: Jared
"""

import numpy as np
import matplotlib.pyplot as plt 
from astropy.table import Table, Column, vstack
import os
import astropy.io.ascii as at
import linecache as lc


# Path of files
input_path = os.path.expanduser(r"G:/Shared drives/DouglasGroup/models/BHAC15")
output_path = os.path.expanduser(r"G:/Shared drives/DouglasGroup/Jared Sofair 2022/CSV Files/BHAC Tables")


# Try numpy loadtext
input_file_tmass = os.path.join(input_path, r"BHAC15_iso.2mass")
input_file_gaia = os.path.join(input_path, r"BHAC15_iso.GAIA")

bhac_tmass = os.path.join(output_path, r"BHAC_2MASS.csv")
bhac_gaia = os.path.join(output_path, r"BHAC_GAIA.csv")


#%% 2MASS
# Reformat the table
age_tmass = "0.8000"
data_tmass = Table()
with open(input_file_tmass, 'r') as f_tmass:
    l = f_tmass.readline()
    while l != "":
        
        if age_tmass in l:
            # Add age column to headers
            l = f_tmass.readline()
            l = f_tmass.readline()
            colnames = l[1:].split()
            colnames.append("Age")
            data_tmass = Table(names=colnames)
            
            # Add age column to data and format data as astropy table
            l = f_tmass.readline()
            l = f_tmass.readline()
            
            while "!" not in l:
                lsplit = l.split()
                lsplit.append(age_tmass)
                data_tmass.add_row(lsplit)
                l = f_tmass.readline()
                
        else:
            l = f_tmass.readline()
    
    # Write data to csv :D
    data_tmass.write(bhac_tmass, overwrite=True)

#%% GAIA
# Reformat the table
age_gaia = "0.8000"
data_gaia = Table()
with open(input_file_gaia, 'r') as f_gaia:
    l = f_gaia.readline()
    while l != "":
        
        if age_gaia in l:
            # Add age column to headers
            l = f_gaia.readline()
            l = f_gaia.readline()
            colnames = l[2:].split()
            colnames.append("Age")
            data_gaia = Table(names=colnames)
            
            # Add age column to data and format data as astropy table
            l = f_gaia.readline()
            l = f_gaia.readline()
            
            while "!" not in l:
                lsplit = l.split()
                lsplit.append(age_gaia)
                data_gaia.add_row(lsplit)
                l = f_gaia.readline()
                
        else:
            l = f_gaia.readline()
    
    # Write data to csv :D
    data_gaia.write(bhac_gaia, overwrite=True)
#%%

# Can also put it in Excel and reformat it there


# Read in tables
# tmass_initial = at.read(tmass_file, comment='!', data_start=759, data_end=782, format='no_header', guess=False, fast_reader=False, names=names_tmass)#, names = names_tmass) # Can't start var name with a # :(
# gaia = at.read(os.path.join(bhac_path, r'BHAC15_iso.GAIA'), data_start=821, data_end=844, format='basic', delimiter=' ', guess=False, fast_reader=False, names = names_gaia)
# print(tmass_initial)
    
    # lc might open and close the file everytime you open it? if so, just use f.readline