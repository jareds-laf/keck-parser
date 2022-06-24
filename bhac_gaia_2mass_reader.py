# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 11:50:59 2022

@author: Jared
"""

import numpy as np
import matplotlib.pyplot as plt 
from astropy.table import Table, Column, vstack, join
import os
import astropy.io.ascii as at
import linecache as lc
from numpy import genfromtxt


# Path of files
input_path = os.path.expanduser(r"G:/Shared drives/DouglasGroup/models/BHAC15")
output_path = os.path.expanduser(r"G:/Shared drives/DouglasGroup/Jared Sofair 2022/CSV Files/BHAC Tables")


# Try numpy loadtext
input_file_tmass = os.path.join(input_path, r"BHAC15_iso.2mass")
input_file_gaia = os.path.join(input_path, r"BHAC15_iso.GAIA")

bhac_tmass = os.path.join(output_path, r"BHAC_2MASS.csv")
bhac_gaia = os.path.join(output_path, r"BHAC_GAIA.csv")


#%% Read 2MASS BHAC
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

#%% Read GAIA BHAC
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
#%% Get absG mag and BP-RP from praesepe_merged (meant to be run once just to create the appropriate table)
# csv_path = r"C:\Users\Jared\Documents\GitHub\data-parser\CSV Files"
# pm = Table.read(os.path.join(csv_path, r'praesepe_merged.csv'))
# targets_gaia = os.path.expanduser(os.path.join(csv_path, r'targets_gaia_dr2.csv'))
# targets_2mass = os.path.expanduser(os.path.join(csv_path, r'targets_2mass.csv'))
# targets_names = Table.read(os.path.expanduser(os.path.join(csv_path, r'targets_names.csv')))

# table_gaia = Table.read(targets_gaia, names=["Name"], data_start=0)
# table_2mass = Table.read(targets_2mass, names=["Name"], data_start=0)

# pm_row_index_list = []
# apparentG_list = []
# BPminRP_list = []
# absG_list = []

# names_gaia = table_gaia["Name"][1:].data
# names_2mass = table_2mass["Name"][1:].data

#%%% Get list of indeces for each target star using 2MASS names
# targets_mag_color = Table()

# j=0
# for star in names_2mass:
#     pm_row_index_list.append(np.where(pm["name"]==star))
#     # print(pm_row_index_list[j][0])

#     pos = pm_row_index_list[j][0]
    
#     # Get list of apparentG mags for each target
#     apparentG_list.append(pm["G"][pos])

#     # Get BP-RP for each target
#     BPminRP_list.append(float(pm["BP"][pos]-pm["RP"][pos]))
    
#     j+=1


#%%% Convert apparent G mag to absolute G mag
# k=0
# for value in apparentG_list:
#     pos = pm_row_index_list[k][0]
#     absG_mag = value - 5*np.log10(pm['D'][pos]) + 5 # absmag = appmag - 5*log(D) + 5
#     absG_list.append(absG_mag[0])
#     k+=1

# targets_mag_color.add_column(names_2mass, name="desig_2mass")
# targets_mag_color.add_column(BPminRP_list, name="BP-RP")
# targets_mag_color.add_column(absG_list, name="absG")

# targets_abr = join(targets_names, targets_mag_color, join_type="outer")


# targets_abr.write(os.path.expanduser(os.path.join(csv_path, r"targets_abr.csv")), overwrite=True)

# targets_mag_color.write(os.path.expanduser(os.path.join(csv_path, r"targets_mag_color.csv")), overwrite=True)