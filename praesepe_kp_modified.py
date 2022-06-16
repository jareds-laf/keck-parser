#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 10:17:45 2021

@author: emilybt
"""
''' Input row number here'''
#row_index = 3
name1 = "JS355"
# name1_2mass = "2MASS J08403789+2020178"
# name2 = "KW569"
# name2_2mass = "2MASS J08393715+1948580"


import numpy as np
import os
import linecache
from astropy.table import Table, Column
import astropy.io.ascii as at
import matplotlib.pyplot as plt

import pandas as pd

from scipy.interpolate import interp1d
from astropy.table import join, Table


names_col1 = ["Name", "Epoch", "Filter", "N_obs", "t_int", "150", "200", "250", "300", "400", "500", "700", "1000", "1500", "2000", "PI"]

keck_psf_dir = os.path.expanduser(r'G:/Shared drives/DouglasGroup/data/keck_psf_detections_praesepe')
praesepe_data = at.read(os.path.join(keck_psf_dir, 'detlimtable_praesepe.txt'), data_start=0, delimiter="&", names= names_col1, fill_values=[("...", np.nan)])

# keck_masking_dir = os.path.expanduser(r'G:/Shared drives/DouglasGroup/data/keck_masking_detections_praesepe')
# star = at.read(os.path.join(keck_masking_dir, 'outsaves_181203/JS355_0227binary.txt'))


praesepe_new = (praesepe_data["PI"] == "Douglas") & (praesepe_data["Filter"] == "Kp") #Kp is a filter in near infrared
our_praesepe_data = praesepe_data[praesepe_new]

#making numpy array of the delta mag values
our_praesepe_data_new = np.array([our_praesepe_data["150"], our_praesepe_data["200"], our_praesepe_data["250"], our_praesepe_data["300"], our_praesepe_data["400"], our_praesepe_data["500"], our_praesepe_data["700"], our_praesepe_data["1000"], our_praesepe_data["1500"], our_praesepe_data["2000"]])

#print(our_praesepe_data_new)

#Seperation values in mas from column headers
seperation_values = np.array((150, 200, 250, 300, 400, 500, 700, 1000, 1500, 2000))

#print(our_praesepe_data_new[:, 1].astype(float))

row_index1 = np.where(our_praesepe_data["Name"]==name1)[0]
# row_index2 = np.where(our_praesepe_data["Name"]==name2)[0]

#Selects data from corresponding row
col1 = our_praesepe_data_new[:, row_index1].astype(float)
col2 = seperation_values
#data_to_plot = np.column_stack((col1, col2))
#print(data_to_plot)

#next star on the list
# col3 = our_praesepe_data_new[:, (row_index2)].astype(float)


# Read in detlimtable_praesepe.txt
names_col = ["Epoch", "15", "30", "60", "120", "200", "280"]

keck_masking_dir = r'G:/Shared drives/DouglasGroup/data/keck_masking_detections_praesepe/'
filedir = os.path.join(keck_masking_dir, r'outsaves_181203/JS355_0227binary.txt')

# Establish what the data is
raw_data = linecache.getline(filedir, 16)
data_as_string = raw_data[raw_data.find(r'&')+1:raw_data.find(r'\\')] # Skips the "99% only" and runs until the "\\" at the end of the line
data_as_arr = data_as_string.split('&')

# Convert data from str to float
data = np.asarray(data_as_arr, dtype="float")

# Create the table from the data
t = Table(names = names_col)
t.add_row(data)

sep_vals = np.asarray(t.colnames[1:],dtype="float")


kp_mags = data[1:]
# print(sep_vals)
# print()
# print(kp_mags)

plt.close('all')

#This block makes K mag vs Seperation plots

#Making the plot

#Making legend with stars' names
#leg1 = str(our_praesepe_data[row_index1][0])
#leg2 = str(our_praesepe_data[(row_index2)][0])

fig1, ax1 = plt.subplots()
ax1.plot(col2, col1, color = 'c', marker = 'o', label = name1)
ax1.plot(sep_vals, kp_mags, marker = 'o', color='red', label = name1)

# ax1.plot(sep_vals, kp_mags, marker = 'o', color='red', label = name1)

# plt.figure()
# plt.plot(col2, col1, color = 'c', marker = 'o', label = name1)
# plt.plot(col2, col3, color = 'r', marker = 'o', label = name2)

# ax1.xaxis.set_ticks(np.arrange[15,2000,15])

# ax1.set_autoscale_on(True)

plt.title('\u0394 Kp Magnitude vs. Seperation', fontsize = 12)
ax1.set_xlabel('Seperation Values (mas)')
ax1.set_ylabel('\u0394 Kp Magnitudes')  #\u0394 is a Delta symbol
ax1.invert_yaxis()
# plt.xlim(0,200)
# plt.xlim([15,2000])
# plt.xticks([15, 30, 60, 120, 150, 200, 250, 280, 300, 400, 500, 700, 1000, 1500, 2000])

# plt.xscale("log")

plt.legend()

ax1.plot()

plt.grid()
plt.show()


'''
#Crossmatching Stars
keck_targets_matched = at.read(os.path.join(keck_psf_dir, 'Keck_Targets_2018B_Simbad.csv'), delimiter=",")

#Combining targets by name
combined_keck_targets_matched = join(our_praesepe_data, keck_targets_matched, keys="Name")


#From GvsColor_07_02_21.py
#This part creates the age model
spot_dir = os.path.expanduser("~/Google_Drive/COLLEGE/Emily_Summer_2021/SPOTS_Somers2020/")
spots = at.read(os.path.join(spot_dir, "f000.isoc"),header_start=2,comment="##")


log10_age = np.log10(790e6)
model_ind = np.argmin(abs(spots["logAge"]-log10_age))
model_age = spots["logAge"][model_ind]
model_age_Myr = 10**(model_age-6)
model = (spots["logAge"]==model_age) & (spots["G_mag"]>0)

#Don't ask why I chose to use data frames here, I don't know :|
df_G_mag = pd.DataFrame(spots["G_mag"][model])
df_BP_minus_RP = pd.DataFrame((spots["BP_mag"][model]-spots["RP_mag"][model]))
model_df = pd.concat([df_G_mag, df_BP_minus_RP], axis=1)

df_Mass = pd.DataFrame(spots["Mass"][model])
df_Mass = np.asarray(df_Mass).squeeze()
df_BP_minus_RP = np.asarray(df_BP_minus_RP).squeeze()


#Interpolation function of mass from color based off model
interp_func1 = interp1d(df_BP_minus_RP, df_Mass, bounds_error=False, fill_value=np.nan)

#importing preasepe_laest_AN file to get Bpsw and Rps
praesepe_latest_dir = os.path.expanduser('~/Google_Drive/COLLEGE/Emily_Summer_2021/')
praesepe_latest_an = at.read(os.path.join(praesepe_latest_dir, 'Praesepe_latest_AN.csv'), delimiter=",")
#renaming so columns have the same name
praesepe_latest_an.rename_column("name", "2MASS")

crossmatched_colors = join(combined_keck_targets_matched, praesepe_latest_an, keys ="2MASS")

new_BP_RP = (crossmatched_colors["BP"] - crossmatched_colors["RP"])
named_BP_RP = [new_BP_RP, crossmatched_colors["2MASS"]]
#print(named_BP_RP)

#getting masses from colors
interpolated_masses = interp_func1(named_BP_RP[0])
named_interpolated_masses = Table({"2MASS Names": named_BP_RP[1], "BP-RP": named_BP_RP[0], "Interpolated Masses": interpolated_masses})

#print(named_interpolated_masses)


row_index3 = np.where(named_interpolated_masses["2MASS Names"]==name1_2mass)[0]
row_index4 = np.where(named_interpolated_masses["2MASS Names"]==name2_2mass)[0]
#print(row_index3, row_index4)

mass1 = interpolated_masses[row_index3].astype(float)
mass2 = interpolated_masses[row_index4].astype(float)
#print(mass1, mass2)

second_model = at.read(os.path.join(spot_dir, "sonora_bobcat_nc+0.0_co1.0_mass_0.8Gyr_interpolated.csv"),header_start=0,comment="##")
add_on_model = second_model["M/Msun_phot","Ks","J"]

#mass_ref = interpolated_masses[(len(interpolated_masses))-1]

#importing function to get mass "ratio" from ref mass from another file
#from kmag_jmag_func import k_j_mags


#making dataframes into arrays
j_mag_array = np.append(spots["J_mag"][model], add_on_model["J"])
k_mag_array = np.append(spots["K_mag"][model], add_on_model["Ks"])
mass_array = np.append((spots["Mass"][model]), add_on_model["M/Msun_phot"])  #Depends on ref star mass

#Combines the 3 arrays into one array
#arr = np.stack((mass_array,j_mag_array, k_mag_array), axis=1)
arr = Table({"J": j_mag_array, "K": k_mag_array, "Mass": mass_array})
#Selects all rows where mass is less than or equal to the refrence star
#arr_new = np.vstack((arr[arr[:,0] <=1]))

arr.sort("Mass")
arr_new= arr[arr["Mass"] <= 1]

#Makes one-column arrays
new_j_mag = arr_new["J"]
new_k_mag = arr_new["K"]
new_mass = arr_new["Mass"]
   
#Getting refrence star K and J values
k_ref = new_k_mag[(len(new_k_mag))-1]
j_ref = new_j_mag[(len(new_j_mag))-1]

#Finding difference in k and j magnitudes
delta_k_mag = (new_k_mag - k_ref)
delta_j_mag = (new_j_mag - j_ref)

#print(delta_k_mag, new_mass)
interp_func2 = interp1d(delta_k_mag, new_mass, bounds_error=False, fill_value=np.nan)
#print(interp_func2(1.5))

#interp1 returns mass from color
#interp2 returns mass ratio from delta k magnitudes

#making them into ratios
mass_ratios1 = (interp_func2(col1))/mass1
mass_ratios2 = (interp_func2(col3))/mass2

#print(mass_ratios)

plt.close('all')

#Making the plot

#Makign legend with stars' names
#leg1 = str(our_praesepe_data[row_index1][0])
#leg2 = str(our_praesepe_data[(row_index2)][0])

plt.figure()
plt.plot(col2, mass_ratios1, color = 'c', marker = 'o', label = name1)
plt.plot(col2, mass_ratios2, color = 'r', marker = 'o', label = name2)

plt.title('Mass Ratio vs. Seperation', fontsize = 12)
plt.xlabel('Seperation Values (mas)')
plt.ylabel(r'Mass Ratio $ \frac{M_2}{M_1} $')

plt.legend()

plt.grid()
plt.show()

'''