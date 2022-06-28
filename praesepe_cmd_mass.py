# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 11:36:18 2022

@author: Jared
"""
import os
import numpy as np
import matplotlib.pyplot as plt 
from astropy.table import Table, Column
from scipy.interpolate import interp1d

# Read in necessasry csv files
csv_path_drive = os.path.expanduser(r"G:/Shared drives/DouglasGroup/Jared Sofair 2022/CSV Files")
csv_path_github = os.path.expanduser(r"C:/Users/Jared/Documents/GitHub/data-parser/CSV Files")
pm = Table.read(os.path.join(csv_path_drive, r'Praesepe_Merged.csv'))
bhac_gaia = Table.read(os.path.join(csv_path_drive, r"BHAC Tables/BHAC_GAIA.csv"))
bhac_2mass = Table.read(os.path.join(csv_path_drive, r"BHAC Tables/BHAC_2MASS.csv"))
targets = Table.read(os.path.join(csv_path_github, r"targets_abr.csv"))

# Creating color index columns
GminusK_pm = Column(data=pm['G']-pm['K'])
BPminusRP_pm = Column(data=pm['BP']-pm['RP'])
pm.add_column(GminusK_pm, name = 'G-K')
pm.add_column(BPminusRP_pm, name = 'BP-RP')

BPminusRP_gaia = Column(data=bhac_gaia['G_BP']-bhac_gaia['G_RP'])
bhac_gaia.add_column(BPminusRP_gaia, name = 'BP-RP')

# Creating column with absolute magnitudes
absG_pm = pm['G'] - 5*np.log10(pm['D']) + 5    # absmag = appmag - 5*log(D) + 5
pm.add_column(absG_pm, name='absG')

# Writing pm to a .csv file to ensure that I did things properly
# pm.write(r'G:/Shared drives/DouglasGroup/Jared Sofair 2022/CSV Files/TestingPM.csv', overwrite = True) # Make sure you close this .csv file before running!

#%% Removing irrelevant stars

# Bad ranges
low1 = pm['BP-RP'] <= 2.25
high1 = pm['absG'] >= 10.0

low2 = pm['BP-RP'] <= 2.75
high2 = pm['absG'] >= 12

low3 = pm['BP-RP'] >= 1.0
high3 = pm['absG'] <= 4

low4 = pm['BP-RP'] >= 0.0
high4 = pm['absG'] <= 0.75

low5 = bhac_gaia['BP-RP'] >= 3.0
high5 = bhac_gaia['G'] >= 12.9

bad_pm = (high1 & low1) | (high2 & low2) | (high3 & low3) | (low4 & high4)
# print(np.where(bad_pm)[0])
pm.remove_rows(np.where(bad_pm)[0])

bad_gaia = high5 & low5

bhac_gaia.remove_rows(np.where(bad_gaia)[0])

#%% Creating bins by color index

# Setting bins
# color_bin1 = np.trunc(pm['G-K'] / 0.5) 
# color_bin2 = np.trunc(pm['BP-RP'] / 0.5)

# pm_grouped_color1 = pm.group_by(color_bin1)
# pm_grouped_color2 = pm.group_by(color_bin2)

# pm_binned_color1 = pm_grouped_color1.groups.aggregate(np.mean)
# pm_binned_color2 = pm_grouped_color2.groups.aggregate(np.mean)

# # Setting a particular bin to a different color
# low_bin1 = pm['absG'] >= pm_binned_color1['absG'][4]
# high_bin1 = pm['absG'] <= pm_binned_color1['absG'][5]

# color_bin1 = low_bin1 & high_bin1

#%% Plotting G vs. BP-RP
# fig, ax = plt.subplots()
# ax.plot(pm['BP-RP'], pm['absG'], '.', color='#695ffa')
# # ax.plot(pm['BP-RP'][bad_pm], pm['absG'][bad_pm], 'o', color='red') # Uncomment this and comment line 65 to see removed data

# ax.plot(bhac_gaia['BP-RP'], bhac_gaia['G'], 'o', color='pink')
# # ax.plot(bhac_gaia['BP-RP'][bad_gaia], bhac_gaia['G'][bad_gaia], 'o', color='red') # Uncomment this and comment line 69 to see removed data



# # ax.plot(pm_binned_color2['BP-RP'], pm_binned_color2['absG'], 'o', color='orange')
# # ax.plot(pm['BP-RP'][color_bin2], pm['absG'][color_bin2], '.', color='green')


# plt.title("G vs. (BP-RP)")
# ax.set_xlabel('(BP-RP)')
# ax.set_ylabel('Absolute G Magnitude')
# ax.invert_yaxis()