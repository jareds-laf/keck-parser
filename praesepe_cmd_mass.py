# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 11:36:18 2022

@author: Jared
"""

import numpy as np
import matplotlib.pyplot as plt 
from astropy.table import Table, Column

# Read Praesepe_Merged.csv
pm = Table.read(r'C:\Users\Jared\Desktop\EXCEL Summer 2022\Catalogs\Praesepe_Merged.csv')

# Creating G-K and BP-RP columns
GminusK = Column(data=(pm['G']-pm['K']))
BPminusRP = Column(data=(pm['BP']-pm['RP']))

pm.add_column(GminusK, name = 'G-K')
pm.add_column(BPminusRP, name = 'BP-RP')

# Creating column with absolute magnitudes
absG = pm['G'] - 5*np.log(pm['D']) + 5    # absmag = appmag - 5*log(D) + 5
pm.add_column(absG, name='absG')

# Removing non-main-sequence stars

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

bad = (high1 & low1) | (high2 & low2) | (high3 & low3) | (low5 & high5)
# print(np.where(bad)[0])
pm.remove_rows(np.where(bad)[0])

# Creating bins by color index

# Setting bins
color_bin1 = np.trunc(pm['G-K'] / 0.5) 
color_bin2 = np.trunc(pm['BP-RP'] / 0.5)

pm_grouped_color1 = pm.group_by(color_bin1)
pm_grouped_color2 = pm.group_by(color_bin2)

pm_binned_color1 = pm_grouped_color1.groups.aggregate(np.mean)
pm_binned_color2 = pm_grouped_color2.groups.aggregate(np.mean)

# Setting a bin to a different color
low_red1 = pm['absG'] >= pm_binned_color1['absG'][4]
high_red1 = pm['absG'] <= pm_binned_color1['absG'][5]

color_r = low_red1 & high_red1
# Plotting G vs. G-K

# Creating 1st plot
fig1, ax1 = plt.subplots()
ax1.plot(pm['G-K'], pm['absG'], '.', color='#5fa2fa')
ax1.plot(pm['G-K'][color_r], pm['absG'][color_r], '.', color='red')

ax1.plot(pm_binned_color1['G-K'], pm_binned_color1['absG'], 'o', color='orange')
# ax1.plot(pm_binned_color1['G-K'][color_r], pm_binned_color1['absG'][color_r], 'o', color='red')
# ax1.plot(pm['G-K'][bad], pm['absG'][bad], 'o', color='red') # Uncomment to see removed data

plt.title("G vs. (G-K)")
ax1.set_xlabel('(G-K)')
ax1.set_ylabel('Absolute G Magnitude')
ax1.invert_yaxis()



# Plotting G vs. BP-RP

    # Creating 2nd plot
fig2, ax2 = plt.subplots()
ax2.plot(pm['BP-RP'], pm['absG'], '.', color='#695ffa')

ax2.plot(pm_binned_color2['BP-RP'], pm_binned_color2['absG'], 'o', color='orange')
ax2.plot(pm['BP-RP'][color_r], pm['absG'][color_r], '.', color='red')

# ax2.plot(pm['BP-RP'][bad], pm['absG'][bad], 'o', color='red') # Uncomment to see removed data

plt.title("G vs. (BP-RP)")
ax2.set_xlabel('(BP-RP)')
ax2.set_ylabel('Absolute G Magnitude')
ax2.invert_yaxis()



# Writing pm to a .csv file to ensure that I did everything properly
pm.write(r'G:\Shared drives\DouglasGroup\Jared Sofair 2022\CSV Files\TestingPM.csv', overwrite = True) # Make sure you close this .csv file before running!