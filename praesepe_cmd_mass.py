# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 11:36:18 2022

@author: Jared
"""
import numpy as np
import matplotlib.pyplot as plt 
from astropy.table import Table, Column
import os

# Read in necessasry csv files
csv_path = r"G:/Shared drives/DouglasGroup/Jared Sofair 2022/CSV Files"
pm = Table.read(os.path.join(csv_path, r'Praesepe_Merged.csv'))
bhac_gaia = Table.read(os.path.join(csv_path, r"BHAC Tables/BHAC_GAIA.csv"))
bhac_2mass = Table.read(os.path.join(csv_path, r"BHAC Tables/BHAC_2MASS.csv"))


# Creating color index columns
GminusK_pm = Column(data=(pm['G']-pm['K']))
BPminusRP_pm = Column(data=(pm['BP']-pm['RP']))
pm.add_column(GminusK_pm, name = 'G-K')
pm.add_column(BPminusRP_pm, name = 'BP-RP')


absG_gaia = Column(data=bhac_gaia['G']) # Currently appg, not absg
BPminusRP_gaia = Column(data=bhac_gaia['G_BP']-bhac_gaia['G_RP'])
bhac_gaia.add_column(BPminusRP_gaia, name = 'BP-RP')


# absK_2mass = Column(data=bhac_2mass['Mk']) # Currently appk, not absk
# JminusK_2mass = Column(data=bhac_2mass['Mj']-bhac_2mass['Mk'])
# bhac_2mass.add_column(JminusK_2mass, name = 'J-K')

# Creating column with absolute magnitudes
absG_pm = pm['G'] - 5*np.log(pm['D']) + 5    # absmag = appmag - 5*log(D) + 5
pm.add_column(absG_pm, name='absG')

absG_gaia = bhac_gaia['G']
bhac_gaia.add_column(absG_gaia, name = 'absG')

#%% Removing non-main-sequence stars

# Bad ranges
low1 = pm['BP-RP'] <= 0.0
high1 = pm['absG'] >= -5.0

low2 = pm['BP-RP'] <= 2.75
high2 = pm['absG'] >= -1.0

low3 = pm['BP-RP'] >= 0.6
high3 = pm['absG'] <= -13

low4 = pm['BP-RP'] >= 2.25
high4 = pm['absG'] <= -7.5

low5 = pm['G-K'] <= 6.0
high5 = pm['absG'] >= -1.0

low6 = pm['G-K'] >= 3.0
high6 = pm['absG'] <= -12.0

low7 = pm['G-K'] <= 1.0
high7 = pm['absG'] <= -14.0

low8 = pm['G-K'] >= 2.6
high8 = pm['absG'] <= -7.8

low9 = pm['BP-RP'] >= 1.0
high9 = pm['absG'] <= -14.0

low10 = pm['BP-RP'] >= 1.0
high10 = pm['absG'] <= -14.0

low11 = pm['BP-RP'] <= 2.6
high11 = pm['absG'] >= -4.0

bad = (high1 & low1) | (high2 & low2) | (high3 & low3) | (low5 & high5) | (low6 & high6) | (low7 & high7) | (low8 & high8) | (low9 & high9) | (low10 & high10) | (low11 & high11)
# print(np.where(bad)[0])
pm.remove_rows(np.where(bad)[0])

# Creating bins by color index

#%% Setting bins
color_bin1 = np.trunc(pm['G-K'] / 0.5) 
color_bin2 = np.trunc(pm['BP-RP'] / 0.5)

pm_grouped_color1 = pm.group_by(color_bin1)
pm_grouped_color2 = pm.group_by(color_bin2)

pm_binned_color1 = pm_grouped_color1.groups.aggregate(np.mean)
pm_binned_color2 = pm_grouped_color2.groups.aggregate(np.mean)

# Setting a bin to a different color
low_bin1 = pm['absG'] >= pm_binned_color1['absG'][4]
high_bin1 = pm['absG'] <= pm_binned_color1['absG'][5]

color_bin1 = low_bin1 & high_bin1
#%% Plotting G vs. G-K

# # Creating 1st plot
# fig1, ax1 = plt.subplots()
# ax1.plot(pm['G-K'], pm['absG'], '.', color='#5fa2fa')
# ax1.plot(pm['G-K'][color_bin1], pm['absG'][color_bin1], '.', color='green')

# ax1.plot(pm_binned_color1['G-K'], pm_binned_color1['absG'], 'o', color='orange')
# # ax1.plot(pm_binned_color1['G-K'][color_bin1], pm_binned_color1['absG'][color_bin1], 'o', color='red')
# # ax1.plot(pm['G-K'][bad], pm['absG'][bad], 'o', color='red') # Uncomment to see removed data

# plt.title("G vs. (G-K)")
# ax1.set_xlabel('(G-K)')
# ax1.set_ylabel('Absolute G Magnitude')
# ax1.invert_yaxis()



#%% Plotting G vs. BP-RP

# Creating 2nd plot
fig2, ax2 = plt.subplots()
ax2.plot(pm['BP-RP'], pm['absG'], '.', color='#695ffa')

ax2.plot(pm_binned_color2['BP-RP'], pm_binned_color2['absG'], 'o', color='orange')
ax2.plot(pm['BP-RP'][color_bin1], pm['absG'][color_bin1], '.', color='green')
ax2.plot(bhac_gaia['BP-RP'], bhac_gaia['absG'], '.', color='pink')
# ax2.plot(pm['BP-RP'][bad], pm['absG'][bad], 'o', color='red') # Uncomment to see removed data

plt.title("G vs. (BP-RP)")
ax2.set_xlabel('(BP-RP)')
ax2.set_ylabel('Absolute G Magnitude')
ax2.invert_yaxis()



# Writing pm to a .csv file to ensure that I did everything properly
pm.write(r'G:\Shared drives\DouglasGroup\Jared Sofair 2022\CSV Files\TestingPM.csv', overwrite = True) # Make sure you close this .csv file before running!