# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 16:11:19 2022

@author: Jared
"""

from keck_parser import get_data_masking, get_data_psf, plot_star, export_star
import os
import numpy as np
import linecache as lc
import glob
from astropy.table import Table, join
import astropy.io.ascii as at
import matplotlib.pyplot as plt
import pandas as pd
import re

#%% Former keck_parser.py actions
#%%% Plot all targets
# targets = pd.read_excel(r'G:/Shared drives/DouglasGroup/data/Copy of Keck Targets.xlsx', index_col=0)
# try:
#     fig, ax = plt.subplots()
#     for name, obsdate in targets.iterrows():
#         print("Name:", name, "Obsdate:", obsdate[14])
#         if not pd.isnull(obsdate[14]):
#             star = name.replace(" ", "_")
#             # print(star)
#             plot_star(star, ax)
#         else:
#             print(f"Plot not generated for {name} (no observation date)")
# except KeyboardInterrupt:
#     plt.title("\u0394 Kp Magnitude vs. Separation", fontsize = 12)
#     ax.set_xlabel("Separation Values (mas)")
#     ax.set_ylabel("\u0394 Kp Magnitudes")
#     if ax.get_ylim()[0] < ax.get_ylim()[1]:
#         ax.invert_yaxis()
#     plt.grid(visible = True)
#     plt.show()
#     plt.close("all")
# else:
#     plt.title("\u0394 Kp Magnitude vs. Separation", fontsize = 12)
#     ax.set_xlabel("Separation Values (mas)")
#     ax.set_ylabel("\u0394 Kp Magnitudes")
#     if ax.get_ylim()[0] < ax.get_ylim()[1]:
#         ax.invert_yaxis()
#     plt.grid(visible = True)
#     plt.savefig(os.path.expanduser(r'C:/Users/Jared/Documents/GitHub/data-parser/plots/keck_contrast_curves.pdf'))
#     plt.savefig(os.path.expanduser(r'C:/Users/Jared/Documents/GitHub/data-parser/plots/keck_contrast_curves.png'))
#     plt.show()
#     plt.close("all")

# print(get_data_masking("JS355"))

#%%% Export all targets to its own file that MOLUSC can take
# list_plot = []
# list_exp_masking = []
# list_exp_psf = []

# targets = pd.read_excel(r'G:/Shared drives/DouglasGroup/data/Copy of Keck Targets.xlsx', index_col=0)
# for name, obsdate in targets.iterrows():
#     print("Name:", name, "Obsdate:", obsdate[14])
#     if not pd.isnull(obsdate[14]): # Only try for stars we have data with
        
#         if get_data_masking(name) is not None: # Has masking data
#             star = name.replace(" ", "_")    
#             # export_star(name)
#             print(star, "\n")
#             list_exp_masking.append(star)
#         else: # No masking data :(
#             star = name.replace(" ", "_")
#             print(f"No masking data for {name}")
#             print(star, "\n")
        
#         if get_data_psf(name) is not None: # Has psf data
#             # export_star(name)  
#             print(star, "\n")
#             list_exp_psf.append(star)
#         else: # No psf data :(
#             star = name.replace(" ", "_")
#             print(f" No PSF data for {name}")    
#             print(star, "\n")
            
#     else: # We don't have data!
#         print(f"File not generated for {name} (no observation date)\n")

# print("# of plots:", len(list_plot))
# print("# of exports:", len(list_exp_psf), len(list_exp_masking))
# print("PSF exports:", list_exp_psf)
#%% Create targets_abr.csv
csv_path_github = os.path.expanduser(r"C:\Users\Jared\Documents\GitHub\data-parser\CSV Files")
csv_path_drive = r"G:/Shared drives/DouglasGroup/Jared Sofair 2022/CSV Files"
simbad = Table.read(os.path.expanduser(r"G:/Shared drives/DouglasGroup/data/keck_psf_detections_praesepe/Keck_Targets_2018B_Simbad.csv"))

pm = Table.read(os.path.join(csv_path_github, r'praesepe_merged.csv'))

names = ['AD_1427', 'AD_2354', 'AD_2595', 'AD_3663', 'EPIC211885995', 'HSHJ510', 'JS117', 'JS169', 'JS178', 'JS181', 'JS191', 'JS301', 'JS317', 'JS318', 'JS352', 'JS355', 'JS373', 'JS394', 'JS405', 'JS432', 'JS513', 'JS533', 'JS545', 'JS552', 'JS582', 'JS620', 'JS649', 'KW569', 'AD_0738', 'AD_3411', 'JS_46', 'JS113', 'JS230', 'JS364', 'JS430', 'JS452', 'KW564', 'AD_1660', 'EPIC211998192', 'EPIC212011416', 'EPIC212127087', 'HSHJ300', 'JS119', 'JS148', 'JS174', 'JS187', 'JS200', 'JS232', 'JS246', 'JS415', 'JS468', 'JS488', 'JS505', 'JS541', 'JS550', 'JS689']

desig_gaia_dr2_all = []
desig_2mass_all = []
simbad_row_index_list = []
pm_index_list = []
targets_abr = Table()

#%%% Match name formatting of pm
j=0
for item in names:
    if item.find("_") != -1:
        item = item.replace("_", " ")
        names[j] = item
    j+=1

# Add name column to targets_abr table after sorting alphabetically
names.sort()
targets_abr.add_column(names, name="name")

#%%% Get simbad indeces
j=0
for star in targets_abr["name"]:
    # Get simbad index
    simbad_row_index_list.append(np.where(simbad["Name"]==star))
    # Next line is necessary because the output style of np.where is horrendous
    simbad_row_index_list[j] = simbad_row_index_list[j][0][0]
    j+=1

# Add simbad indeces to targets_abr
targets_abr.add_column(simbad_row_index_list, name="simbad_index")

# Make sure that stars and simbad indeces are properly paired
j=0
for ind in targets_abr["simbad_index"]:
    star_abr = targets_abr["name"][j]
    star_simbad = simbad["Name"][ind]

    if star_abr != star_simbad:
        print(f"ERROR, stars do not match: {star_abr}, {star_simbad}")
    j+=1

#%%% Get gaia DR2 and 2MASS designations
j=0
for ind in targets_abr["simbad_index"]:
    # Add star designation to appropriate list
    desig_gaia_dr2_all.append(simbad['GaiaDR2'][ind])
    desig_2mass_all.append(simbad['2MASS'][ind])
    
    j+=1

# Add Gaia DR2 and 2MASS designations to targets_abr
targets_abr.add_column(desig_gaia_dr2_all, name="desig_gaia_dr2")
targets_abr.add_column(desig_2mass_all, name="desig_2mass")

# Check if Gaia DR2 and 2MASS designations were properly matched with the target star
j=0
for ind in targets_abr["simbad_index"]:
    star_simbad_g = simbad["GaiaDR2"][ind]
    star_simbad_t = simbad["2MASS"][ind]
    
    star_g = targets_abr["desig_gaia_dr2"][j]
    star_t = targets_abr["desig_2mass"][j]

    if star_g != star_simbad_g:
        print(f"ERROR, Gaia names do not match: {star_g}, {star_simbad_g}")
    
    if star_t != star_simbad_t:
        print(f"ERROR, 2MASS names do not match: {star_t}, {star_simbad_t}")
    # else:
    #     print(f"-------------------------{star_t}")
    j+=1

#%%% Get pm index list using 2MASS names
j=0
for star in targets_abr["desig_2mass"]:
    pm_index_list.append(np.where(pm["name"]==star))
    pm_index_list[j] = pm_index_list[j][0][0]
    
    j+=1

# Add pm index list to targets_abr
targets_abr.add_column(pm_index_list, name="pm_index")

# Check if pm index was properly matched to star using 2MASS and Gaia deisgnations
j=0
for ind in targets_abr["pm_index"]:
    star_g = targets_abr["desig_gaia_dr2"][j]
    star_t = targets_abr["desig_2mass"][j]
    star_pm_g = pm["DR2_desig"][ind]
    star_pm_t = pm["name"][ind]
    # print(star_g, star_pm_g, star_t, star_pm_t, j)
    if star_g != star_pm_g:
        print(f"ERROR, Gaia names do not match: {star_g}, {star_pm_g}")
    
    if star_t != star_pm_t:
        print(f"ERROR, 2MASS names do not match: {star_t}, {star_pm_t}")
    j+=1

targets_abr.write(os.path.join(csv_path_github, "targets_abr2.csv"), overwrite=True)


#%%%% Export the list of targets
# np.savetxt(os.path.expanduser(r'C:/Users/Jared/Documents/GitHub/data-parser/CSV Files/targets_gaia_dr2.csv'), desig_gaia_dr2_all, header="desig_gaia_dr2", fmt='%s', delimiter=',')
# np.savetxt(os.path.expanduser(r'C:/Users/Jared/Documents/GitHub/data-parser/CSV Files/targets_2mass.csv'), desig_2mass_all,  header="desig_2mass", fmt='%s', delimiter=',')
# np.savetxt(os.path.expanduser(r'C:/Users/Jared/Documents/GitHub/data-parser/CSV Files/targets_keck.csv'), all_exports,  header="name", fmt='%s', delimiter=',')

# df = pd.DataFrame({"name" : all_exports, "desig_2mass" : desig_2mass_all, "desig_gaia_dr2" : desig_gaia_dr2_all})
# df.to_csv(os.path.expanduser(r'C:/Users/Jared/Documents/GitHub/data-parser/CSV Files/targets_names.csv'), index=False)




#%%% Get absG mag, absK mag, and BP-RP from praesepe_merged
# csv_path = r"C:\Users\Jared\Documents\GitHub\data-parser\CSV Files"
# targets_abr = Table.read(os.path.join(csv_path, r'targets_abr.csv'))
# targets_names = Table.read(os.path.expanduser(os.path.join(csv_path, r'targets_names.csv')))

# targets_gaia = os.path.expanduser(os.path.join(csv_path, r'targets_gaia_dr2.csv'))
# targets_2mass = os.path.expanduser(os.path.join(csv_path, r'targets_2mass.csv'))
# table_gaia = Table.read(targets_gaia, names=["Name"], data_start=0)
# table_2mass = Table.read(targets_2mass, names=["Name"], data_start=0)

# pm_row_index_list = []
# apparentG_list = []
# apparentK_list = []
# BPminRP_list = []
# absG_list = []
# absK_list = []
# pos_list = []

# names_gaia = table_gaia["Name"][1:].data
# names_2mass = table_2mass["Name"][1:].data

#%%% Get list of indeces for each target star using 2MASS names
# targets_mag_color = Table()

# j=0
# for star in targets_abr["desig_2mass"]:
#     pm_row_index_list.append(np.where(pm["name"]==star))
#     pos = pm_row_index_list[j][0]
#     print(pm["name"][pos])
#     pos_list.append(pos)
    
    # pos_list[j] = pos_list[j][0]
    # print(pos_list[j], pm["name"][pos])
    # print(pos_list[j][0])
    
    # pos_list.append(pm_row_index_list[j][0])
    # print(pos)
    
    # print(star)

    
#     # Get list of apparentG mags for each target
#     apparentG_list.append(pm["G"][pos])

#     # Get list of apparentK mags for each target
#     apparentK_list.append(pm["K"][pos])

#     # Get BP-RP for each target
#     BPminRP_list.append(float(pm["BP"][pos]-pm["RP"][pos]))
    
    # j+=1
    
# targets_abr.add_column(pos_list, name="pm_index")

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

#%%% Convert apparent K mag to absolute K mag
# print(targets_abr["desig_2mass"])
# print(pm["name"][92])

# n=0
# for ind in targets_abr["pm_index"]:
#     apparentK = pm["K"][ind]
#     # print(apparentK)
#     absK_mag = apparentK - 5*np.log10(pm['D'][ind]) + 5 # absmag = appmag - 5*log(D) + 5
#     # print(absK_mag)
#     absK_list.append(absK_mag)
#     n+=1
# print(absK_list)

# targets_mag_color.add_column(absG_list, name="absK")


# targets_abr = join(targets_names, targets_mag_color, join_type="outer")

# targets_abr.add_column(absK_list, name="absK")
# targets_abr.write(os.path.expanduser(os.path.join(csv_path, r"targets_abr.csv")), overwrite=True)

# targets_mag_color.write(os.path.expanduser(os.path.join(csv_path, r"targets_mag_color.csv")), overwrite=True)