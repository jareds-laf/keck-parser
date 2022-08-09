# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 16:11:19 2022

@author: Jared
"""

from keck_parser import get_data_masking, get_data_psf, plot_star, export_star
import os
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
import tkinter
from tkinter import filedialog

#%% Former keck_parser.py actions
#%%% Plot all targets
targets = pd.read_excel(r'G:/Shared drives/DouglasGroup/data/Copy of Keck Targets.xlsx', index_col=0)
try:
    fig, ax = plt.subplots()
    for name, obsdate in targets.iterrows():
        print("Name:", name, "Obsdate:", obsdate[14])
        if not pd.isnull(obsdate[14]):
            star = name.replace(" ", "_")
            # print(star)
            plot_star(star, ax)
        else:
            print(f"Plot not generated for {name} (no observation date)")
except KeyboardInterrupt:
    plt.title("\u0394 Kp Magnitude vs. Separation", fontsize = 12)
    ax.set_xlabel("Separation Values (mas)")
    ax.set_ylabel("\u0394 Kp Magnitudes")
    if ax.get_ylim()[0] < ax.get_ylim()[1]:
        ax.invert_yaxis()
    plt.grid(visible = True)
    ax.set_xscale('log')
    plt.show()
    plt.close("all")
else:
    plt.title("\u0394 Kp Magnitude vs. Separation", fontsize = 12)
    ax.set_xlabel("Separation Values (mas)")
    ax.set_ylabel("\u0394 Kp Magnitudes")
    if ax.get_ylim()[0] < ax.get_ylim()[1]:
        ax.invert_yaxis()
    plt.grid(visible = True)
    ax.set_xscale('log')
    plt.savefig(os.path.expanduser(r'C:/Users/Jared/Documents/GitHub/data-parser/plots/keck_contrast_curves.pdf'))
    plt.savefig(os.path.expanduser(r'C:/Users/Jared/Documents/GitHub/data-parser/plots/keck_contrast_curves.png'))
    plt.show()
    plt.close("all")

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

#%%% Export the list of targets (originally used to create the inaccurate targets_abr table)
# np.savetxt(os.path.expanduser(r'C:/Users/Jared/Documents/GitHub/data-parser/CSV Files/targets_gaia_dr2.csv'), desig_gaia_dr2_all, header="desig_gaia_dr2", fmt='%s', delimiter=',')
# np.savetxt(os.path.expanduser(r'C:/Users/Jared/Documents/GitHub/data-parser/CSV Files/targets_2mass.csv'), desig_2mass_all,  header="desig_2mass", fmt='%s', delimiter=',')
# np.savetxt(os.path.expanduser(r'C:/Users/Jared/Documents/GitHub/data-parser/CSV Files/targets_keck.csv'), all_exports,  header="name", fmt='%s', delimiter=',')

# df = pd.DataFrame({"name" : all_exports, "desig_2mass" : desig_2mass_all, "desig_gaia_dr2" : desig_gaia_dr2_all})
# df.to_csv(os.path.expanduser(r'C:/Users/Jared/Documents/GitHub/data-parser/CSV Files/targets_names.csv'), index=False)


#%% Create targets_abr.csv
# Project to make this faster: use array operations instead of for loops here
csv_path_github = os.path.expanduser(r"C:/Users/Jared/Documents/GitHub/data-parser/CSV Files")
csv_path_drive = os.path.expanduser(r"G:/Shared drives/DouglasGroup/Jared Sofair 2022/CSV Files")
simbad = Table.read(os.path.expanduser(r"G:/Shared drives/DouglasGroup/data/keck_psf_detections_praesepe/Keck_Targets_2018B_Simbad.csv"))
bhac_gaia = Table.read(os.path.join(csv_path_drive, r"BHAC Tables/BHAC_GAIA.csv"))
bhac_2mass = Table.read(os.path.join(csv_path_drive, r"BHAC Tables/BHAC_2MASS.csv"))

pm = Table.read(os.path.join(csv_path_github, r'praesepe_merged.csv'))

# HSHJ510 removed since it is not a member of Praesepe
names = ['AD_1427', 'AD_2354', 'AD_2595', 'AD_3663', 'EPIC211885995', 'JS117', 'JS169', 'JS178', 'JS181', 'JS191', 'JS301', 'JS317', 'JS318', 'JS352', 'JS355', 'JS373', 'JS394', 'JS405', 'JS432', 'JS513', 'JS533', 'JS545', 'JS552', 'JS582', 'JS620', 'JS649', 'KW569', 'AD_0738', 'AD_3411', 'JS_46', 'JS113', 'JS230', 'JS364', 'JS430', 'JS452', 'KW564', 'AD_1660', 'EPIC211998192', 'EPIC212011416', 'EPIC212127087', 'HSHJ300', 'JS119', 'JS148', 'JS174', 'JS187', 'JS200', 'JS232', 'JS246', 'JS415', 'JS468', 'JS488', 'JS505', 'JS541', 'JS550', 'JS689']

desig_gaia_dr2_all = []
desig_2mass_all = []
simbad_row_index_list = []
pm_index_list = []
absK_list = []
apparentK_list = []
absG_list = []
BPminRP_list = []
ra_list = []
de_list = []
testing = []
ra_j2k_list = []
de_j2k_list = []
targets_abr = Table()
lamost_search = Table()


#%%% Match name formatting of pm
for j, item in enumerate(names):
        if item.find("_") != -1:
            item = item.replace("_", " ")
            names[j] = item

        
# Add name column to targets_abr table after sorting alphabetically
names.sort()

# testing = np.where(names.find("_") != 1, names = item.replace("_", " "))
# print(testing.sort())
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
targets_abr.add_columns([desig_gaia_dr2_all, desig_2mass_all], names=["desig_gaia_dr2", "desig_2mass"])

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
    if star_g != star_pm_g:
        print(f"ERROR, Gaia names do not match: {star_g}, {star_pm_g}")
    
    if star_t != star_pm_t:
        print(f"ERROR, 2MASS names do not match: {star_t}, {star_pm_t}")
    j+=1


#%%% Calculate absolute K and absolute G magnitudes using 2MASS names
j=0
for ind in targets_abr["pm_index"]:
    # apparent K --> absolute K
    apparentK = pm["K"][ind]
    apparentK_list.append(apparentK) # Used later to check if mass and stars were properly paired
    absK = apparentK - 5*np.log10(pm["D"][ind]) + 5 # absmag = appmag - 5*log(D) + 5
    absK_list.append(absK)
    
    # apprent G --> absolute G
    apparentG = pm["G"][ind]
    absG = apparentG - 5*np.log10(pm["D"][ind]) + 5 # absmag = appmag - 5*log(D) + 5
    absG_list.append(absG)
    j+=1

# Add absK and absG to targets_abr
targets_abr.add_columns([absK_list, absG_list], names=["absK", "absG"])
targets_abr.add_column(apparentK_list, name="apparentK")

#%%% Calculate BP-RP color index
j=0
for ind in targets_abr["pm_index"]:
    BPminRP = pm["BP"][ind] - pm["RP"][ind]
    # print(BPminRP, ind)
    BPminRP_list.append(BPminRP)
    # j+=1

targets_abr.add_column(BPminRP_list, name="BP-RP")
    
    

#%%% Calculate masses
# This section is taken from praesepe_cmd_mass.py (slightly modified here)
#%%%% Fit data to BHAC models and use interpolation to get masses
x1 = bhac_gaia["G_BP"]-bhac_gaia["G_RP"]
x2 = bhac_gaia["G"]
y = bhac_gaia["M/Ms"]

# Get mass using color
calc_mass_gaia_color = interp1d(x1, y)
mass_color = calc_mass_gaia_color(targets_abr["BP-RP"])
# print(mass_color)

# Get mass using absolute magnitude
calc_mass_gaia_absmag = interp1d(x2, y)
mass_absmag = calc_mass_gaia_absmag(targets_abr["absG"])
# print(mass_absmag)

x_K = bhac_2mass["Mk"]
y_K = bhac_2mass["M/Ms"]

calc_mass_2mass = interp1d(x_K, y_K)
mass_K = calc_mass_2mass(targets_abr["absK"])

targets_abr.add_column(mass_K, name="M/Ms")

#%%%% Plotting masses obtained with color and those obtained with absolute magnitude to see which is better
# fig, ax = plt.subplots()
# plt.title("Masses Using Absolute K Magnitude")
# ax.set_xlabel('Mass Absolute K Magnitude')
# ax.set_ylabel('Gaia Data')
# ax.plot(mass_color, mass_K, '.', label="BP-RP")
# ax.plot(mass_absmag, mass_K, '.', label="absG")
# linear = [0, 0.65]
# ax.plot(linear, linear, '-')
# plt.legend()
# plt.show()
# plt.close()

# From the plot, it is clear that getting mass with absolute G magnitude is more accurate than with BP-RP color index
# (at least with our sample)

#%%%% Plotting apparentK mag vs. mass to make sure the values are paired with the correct targets
# fig, ax = plt.subplots()
# plt.title("Mass vs. Apparent K Magnitude")
# ax.set_xlabel("Apparent K Magnitude")
# ax.set_ylabel("Mass")
# ax.plot(apparentK_list, mass_K , '.')
# plt.show()
# plt.close()
#%%% Get coordinates from Keck Targets
#08h49m26.76s -- +18d31m19.5s
# Get ra and dec info
for j, ind in enumerate(targets_abr["simbad_index"]):
    # Get ra
    h = simbad["RAh"][ind]
    m = simbad["RAm"][ind]
    s = np.round(simbad["RAs"][ind], 2)
    ra = f"0{h}h{m}m{s}s"
    # print(f"{simbad['Name'][ind]}: 0{h}h{m}m{s}s {j}")
    ra_list.append(ra)
    
    # Get dec
    d = simbad["DEd"][ind]
    m2 = simbad["DEm"][ind]
    s2 = np.round(simbad["DEs"][ind], 1)
    de = f"{d}d{m2}m{s2}s"
    # print(f"{simbad['Name'][ind]}: {de} {j}")
    de_list.append(de)
    
# print(simbad["RAh"][np.where(targets_abr["name"]==)])


# Append ra and de columns to targets_abr
targets_abr.add_columns([ra_list, de_list], names=["ra", "de"])

# This checking doesn't work, but I can verify that 8=8 is a true statement, as is 21=21. This code thinks they are false.
# Make sure data pairing was successful
# for j, ind in enumerate(targets_abr["simbad_index"]):
#     m1 = targets_abr["ra"][j][1:2]
#     m2 = targets_abr["de"][j][0:2]
    
#     m1_simbad = simbad["RAh"][ind]
#     m2_simbad = simbad["DEd"][ind]
    
#     # print(simbad['Name'][ind], m1, m2, j)
#     # print(simbad['Name'][ind], m1_simbad, m2_simbad, j)
#     # print(m1)

     
#     if m1 != m1_simbad:
#         print(f"BAD RA {j}")
#     if m2 != m2_simbad:
#         print(f"BAD DEC {j}")
        
# targets_abr.write(os.path.join(csv_path_github, "targets_abr.csv"), overwrite=True)

#%% Get RA_J2000 and DECJ_2000 for LAMOST searches

ra_j2k_list = pm["RA_J2000"][targets_abr["pm_index"]]
de_j2k_list = pm["DEC_J2000"][targets_abr["pm_index"]]

# ra_j2k_list.add_row("blank")
# de_j2k_list.add_row("blank")

lamost_search.add_columns([ra_j2k_list, de_j2k_list], names=["RA_J2000", "DEC_J2000"])
# print(lamost_search.columns[:])

sep_rad = np.full(len(lamost_search), 2.0)
lamost_search.add_column(sep_rad, name="sep")
# print(len(sep_rad))
# print(len(ra_j2k_list))
# print(len(names))
# print(len(targets_abr["pm_index"]))

# lamost_search.show_in_browser()
lamost_search.write(os.path.join(csv_path_github, "lamost_search.csv"), overwrite=True)


# Check if ra_j2k and dec_j2k were properly assigned to each target
for j, i in enumerate(targets_abr["pm_index"]):
    # print(i, j)
    star_ra = lamost_search["RA_J2000"][j]
    pm_ra = pm["RA_J2000"][i]
    
    star_de = lamost_search["DEC_J2000"][j]
    pm_de = pm["DEC_J2000"][i]
    
    if star_ra != pm_ra:
        print(f"ERROR, ra_j2ks do not match: {star_ra}, {pm_ra}")
    
    if star_de != pm_de:
        print(f"ERROR, de_j2ks do not match: {star_de}, {pm_de}")

lamost_search.add_column(names, name="name")
lamost_search.write(os.path.join(csv_path_github, "lamost_search_names.csv"), overwrite=True)

#%% Unsuccessful array operations testing

# testing = np.where(names==simbad["Name"])
# RAh = simbad["RAh"][np.where(simbad("Name")==targets_abr("name"))]
# print(np.where(simbad("Name")==targets_abr("name")))
# ra_list.append(RAh)

# what do i want to do?
# search simbad for the target star
# grab the ra and dec info (3 cols each)
# put it into a format that batch_creator can read

# targets_abr.remove_row()
# targets_array = targets_abr.as_array()


# For array operations, feeding it the whole column of (for example) simbad indeces will access the simbad data in the same order that you feed in the indeces

#%% Tinkering with Tkinter to choose files with file explorer
# tkinter.Tk().withdraw()
# folder_path = filedialog.askdirectory()
