# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 14:07:13 2022

@author: Jared
"""

import os
import numpy as np
import linecache as lc
import glob
from astropy.table import Table
import astropy.io.ascii as at
import matplotlib.pyplot as plt
import pandas as pd
import re

# def get_date(starname): # Get the date for
#     # Get file
#     keck_masking_dir = os.path.expanduser(r'G:/Shared drives/DouglasGroup/data/keck_masking_detections_praesepe/*')
    
#     filelist = glob.glob(os.path.join(keck_masking_dir, f"{starname}*"))
    
#     if len(filelist) == 1: # Proper operation
#         filename = filelist[0]
#     elif len(filelist) == 0: # If no star with the input name is found
#         # print(f"Easton, we have a problem! No such star {starname} detected.")
#         return None
#     elif len(filelist) > 1: # If there is more than one star with the input name
#         print(f"Easton, we have a problem! Multiple masking entries detected for {starname}")
#         print(filelist)
#         return None

def get_data_masking(starname): # Get kp mags for a given star (masking dataset)
    
    # Get file
    keck_masking_dir = os.path.expanduser(r'G:/Shared drives/DouglasGroup/data/keck_masking_detections_praesepe/*')

    filelist = glob.glob(os.path.join(keck_masking_dir, f"{starname}*"))
    
    if len(filelist) == 1: # Proper operation
        filename = filelist[0]
    elif len(filelist) == 0: # If no star with the input name is found
        # print(f"Easton, we have a problem! No such star {starname} detected.")
        return None
    elif len(filelist) > 1: # If there is more than one star with the input name
        print(f"Easton, we have a problem! Multiple masking entries detected for {starname}")
        print(filelist)
        return None

    '''As far as I can tell, the L99 data is always on the 2nd to last line
       The following block would find it in case that ended up being untrue
       Only uncomment this block if my assumption is proven untrue'''
    # Grab the line index where the L99 data is    
    # with open(filedir, 'r') as f:
    #     looking = True
    #     length = len(f.readlines())
    #     while looking:
    #         line = f.readlines()[length, length-1]
    #         if line.find("99% only") != -1:
    #             looking = False
    
    # Get the length of the file
    with open(filename, 'r') as f:
        length = len(f.readlines())
    
    # Turn the data into an array
    try:
        raw_data = lc.getline(filename, length-1)
        data_as_string = raw_data[raw_data.find(r'&')+1:raw_data.find(r'\\')] # Skips the "99% only" and runs until the "\\" at the end of the line
        data_as_arr = data_as_string.split('&')
        data = np.asarray(data_as_arr, dtype="float")
    except ValueError:
        # print(f"No masking data was found for the star {starname}.")
        return None

    # Create the table from the data
    # names_col_masking = ["Epoch", "15", "30", "60", "120", "200", "280"] # Column names for data from keck_masking_detections_praesepe ICE
    
    # sep_vals = np.asarray(t_masking.colnames[1:],dtype="float") (have here ICE)
    # kp_mags_masking = data[1:]

    return data

def get_data_psf(starname): # Get kp mag vs. sep vals for a given star (psf dataste)
    
    # Correcting the name so each star can be properly located
    starname_psf = starname
    if "_" in f"{starname_psf}":
        starname_psf = starname.replace("_", " ")
    elif starname_psf[0:4] == "EPIC" and starname_psf[0:5] != "EPIC ":
        starname_psf = re.sub("EPIC", "EPIC ", starname_psf)

    names_col = ["Name", "Epoch", "Filter", "N_obs", "t_int", "150", "200", "250", "300", "400", "500", "700", "1000", "1500", "2000", "PI"] # Column names for psf data
    
    # Get data path
    keck_psf_dir = os.path.expanduser(r'G:/Shared drives/DouglasGroup/data/keck_psf_detections_praesepe/')
    praesepe_data = at.read(os.path.join(keck_psf_dir, 'detlimtable_praesepe.txt'), data_start=0, delimiter="&", names = names_col, fill_values=[("...", np.nan)])
    
    praesepe_new = (praesepe_data["PI"] == "Douglas") & (praesepe_data["Filter"] == "Kp") # Kp is a filter in near infrared
    our_praesepe_data = praesepe_data[praesepe_new]
    
    # Making numpy array of the delta mag values
    our_praesepe_data_new = np.array([our_praesepe_data["Epoch"], our_praesepe_data["150"], our_praesepe_data["200"], our_praesepe_data["250"], our_praesepe_data["300"], our_praesepe_data["400"], our_praesepe_data["500"], our_praesepe_data["700"], our_praesepe_data["1000"], our_praesepe_data["1500"], our_praesepe_data["2000"]])
    
    row_index_list = np.where(our_praesepe_data["Name"]==starname_psf)
    
    if np.size(row_index_list) == 0: # No entry with the input name
        print(f"No PSF data was found for star {starname_psf}.")
        return None
    elif np.size(row_index_list) > 1: # More than one entry with the input name
        if np.all(get_data_masking(starname) != None): # Check if it has masking data
            data_masking = get_data_masking(starname)
            for entry in range(np.size(row_index_list[0])):
                if data_masking[0] == row_index_list[0][entry]: # Check which entry date matches the masking data date
                    print(row_index_list[0][entry])
        else:
            print(f"Easton, we have a problem! Multiple PSF entries detected for {starname} :(")
            print(row_index_list[0])
            i=0
            for entry in range(np.size(row_index_list[0])):
                i+=1
                print(i, np.asarray(our_praesepe_data_new[:, row_index_list[0][entry]].astype(float)))        
                maximum = np.asarray(our_praesepe_data_new[:, row_index_list[0][0]].astype(float))
                current = np.asarray(our_praesepe_data_new[:, row_index_list[0][entry]].astype(float))
                if maximum[0] < current[0]:
                    maximum = np.asarray(our_praesepe_data_new[:, row_index_list[0][entry]].astype(float))
            print("Highest mag:", maximum)

    else:
        kp_mags_psf = np.asarray(our_praesepe_data_new[:, row_index_list[0]].astype(float))
        kp_mags_psf = np.reshape(kp_mags_psf, len(kp_mags_psf))
        print("All seems good!")
        return kp_mags_psf[1:]

def plot_star(starname, ax=None): # Plot psf and masking curves for a given star
    if ax is None: 
        # Create the plot for the input star
        fig, ax = plt.subplots()
        plt.title("\u0394 Kp Magnitude vs. Seperation", fontsize = 12)
        ax.set_xlabel("Seperation Values (mas)")
        ax.set_ylabel("\u0394 Kp Magnitudes")  #\u0394 is a Delta symbol

    
    # Plot data
    sep_vals_psf = np.array([150, 200, 250, 300, 400, 500, 700, 1000, 1500, 2000])
    sep_vals_masking = np.array([15, 30, 60, 120, 200, 280])
    
    kp_mags_masking = get_data_masking(starname)[1:]
   
    starname_psf = f"{starname}"
    
    if "_" in f"{starname_psf}":
        starname_psf = starname_psf.replace("_", " ")
    elif starname[0:4] == "EPIC" and starname[0:5] != "EPIC ":
        starname_psf = re.sub("EPIC", "EPIC ", starname_psf)
    kp_mags_psf = get_data_psf(f"{starname_psf}")

    
    sep_vals_all = np.concatenate((sep_vals_masking[0 : np.size(sep_vals_masking) - 2], sep_vals_psf))
    
    
    if np.all(kp_mags_psf == None): # If no psf data has been found, do not create a plot
        print(f"No PSF data found; no plot has been generated for {starname}.")
        return None
    elif np.all(kp_mags_masking == None): # If no masking data has been found, do not create a plot
        print(f"No masking data found; no plot has been generated for {starname}.")
        return None
    else:
        # Join masking data (except last 2 points) and psf data into one array and plot
        # Doing so since psf data is better at the overlapping points
        kp_mags_all = np.concatenate((kp_mags_masking[0 : np.size(kp_mags_masking) - 2], kp_mags_psf))
        plt.step(x=sep_vals_all, y=kp_mags_all, color = "#9cffb6", alpha = 0.75,)
        plt.step(x=sep_vals_all, y=kp_mags_all, color = "#9cffb6", alpha = 0.75)

if __name__ == "__main__":
    masking = get_data_masking("JS355")
    masking = get_data_masking("AD_0738")
    print(masking[1:])

    # print()

    # psf = get_data_psf("EPIC211885995")
    # psf = get_data_psf("EPIC211998192")
    # psf = get_data_psf("JS355")
    psf = get_data_psf("AD_0738")
    print(psf)

    # print()

    # plot_star("HSHJ300")
    # plot_star("JS355")
    # plot_star("EPIC211998192")
    # plot_star("AD_0738")

    # List of star names

    # targets = pd.read_excel(r'G:/Shared drives/DouglasGroup/data/Copy of Keck Targets.xlsx', index_col=0)
    # try:
    #     fig, ax = plt.subplots()
    #     for name, obsdate in targets.iterrows():
    #         if not pd.isnull(obsdate[0]):
    #             star = name.replace(" ", "_")
    #             print(star)
    #             plot_star(star, ax)
    # except KeyboardInterrupt:
    #     plt.title("\u0394 Kp Magnitude vs. Seperation", fontsize = 12)
    #     ax.set_xlabel("Seperation Values (mas)")
    #     ax.set_ylabel("\u0394 Kp Magnitudes")
    #     if ax.get_ylim()[0] < ax.get_ylim()[1]:
    #         ax.invert_yaxis()
    #     plt.grid(visible = True)
    #     plt.show()
    #     plt.close("all")
    # else:
    #     plt.title("\u0394 Kp Magnitude vs. Seperation", fontsize = 12)
    #     ax.set_xlabel("Seperation Values (mas)")
    #     ax.set_ylabel("\u0394 Kp Magnitudes")
    #     if ax.get_ylim()[0] < ax.get_ylim()[1]:
    #         ax.invert_yaxis()
    #     plt.grid(visible = True)
    #     plt.show()
    #     plt.close("all")


    # Check why so many PSF plots aren't being generated!
    # Make sure plot_star is working after letting get_data return date as well
    # Export the data to a file that MOLUSC can take
        # Create a separate file for each star. Put them in a separate folder with a good naming convention!