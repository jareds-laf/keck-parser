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


def get_data_masking(starname): # Get kp mags for a given star (masking dataset)
    
    # Get file
    keck_masking_dir = os.path.expanduser(r'G:/Shared drives/DouglasGroup/data/keck_masking_detections_praesepe/*')

    filelist = glob.glob(os.path.join(keck_masking_dir, f"{starname}*"))
    
    if len(filelist) == 1: # Proper operation
        filename = filelist[0]
    elif len(filelist) == 0: # If no star with the input name is found
        print(f"Easton, we have a problem! No such star {starname} detected.")
        return None
    elif len(filelist) > 1: # If there is more than one star with the input name
        print(f"Easton, we have a problem! More than one such star {starname} detected")
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
        print(f"No masking data was found for the star {starname}.")
        return None

    # Create the table from the data
    # names_col_masking = ["Epoch", "15", "30", "60", "120", "200", "280"] # Column names for data from keck_masking_detections_praesepe ICE
    
    # sep_vals = np.asarray(t_masking.colnames[1:],dtype="float") (have here ICE)
    kp_mags_masking = data[1:]

    return kp_mags_masking

def get_data_psf(starname): # Get kp mag vs. sep vals for a given star (psf dataste)
    
    if "_" in f"{starname}":
        starname = starname.replace("_", " ")
    else:
        starname = f"{starname}"

    names_col = ["Name", "Epoch", "Filter", "N_obs", "t_int", "150", "200", "250", "300", "400", "500", "700", "1000", "1500", "2000", "PI"] # Column names for psf data
    
    # Get data path
    keck_psf_dir = os.path.expanduser(r'G:/Shared drives/DouglasGroup/data/keck_psf_detections_praesepe/')
    praesepe_data = at.read(os.path.join(keck_psf_dir, 'detlimtable_praesepe.txt'), data_start=0, delimiter="&", names = names_col, fill_values=[("...", np.nan)])
    
    praesepe_new = (praesepe_data["PI"] == "Douglas") & (praesepe_data["Filter"] == "Kp") # Kp is a filter in near infrared
    our_praesepe_data = praesepe_data[praesepe_new]
    
    # Making numpy array of the delta mag values
    our_praesepe_data_new = np.array([our_praesepe_data["150"], our_praesepe_data["200"], our_praesepe_data["250"], our_praesepe_data["300"], our_praesepe_data["400"], our_praesepe_data["500"], our_praesepe_data["700"], our_praesepe_data["1000"], our_praesepe_data["1500"], our_praesepe_data["2000"]])

    #Seperation values in mas from column headers (ICE)
    # sep_vals = np.array((150, 200, 250, 300, 400, 500, 700, 1000, 1500, 2000))
    
    row_index = np.where(our_praesepe_data["Name"]==starname)[0]
    
    # print(np.where(our_praesepe_data["Name"]==starname)[0])
    
    if len(row_index) == 0:
        print(f"No PSF data was found for star {starname}.")
        return None
    else:
        kp_mags_psf = our_praesepe_data_new[:, row_index].astype(float)
        return kp_mags_psf

def plot_star(starname, ax=None): # Plot psf and masking curves for a given star
    if ax is None: 
        # Create the plot for the input star
        fig, ax = plt.subplots()
        plt.title("\u0394 Kp Magnitude vs. Seperation", fontsize = 12)
        ax.set_xlabel("Seperation Values (mas)")
        ax.set_ylabel("\u0394 Kp Magnitudes")  #\u0394 is a Delta symbol
    
    # Plot psf data
    sep_vals_psf = np.array((150, 200, 250, 300, 400, 500, 700, 1000, 1500, 2000))

    kp_mags_psf = get_data_psf(f"{starname}")     
    
    if np.all(kp_mags_psf == None): # If no psf data has been found, do not create a plot
        print(f"No PSF plot has been generated for {starname}.")
        return None
    else:
        # ax.plot(sep_vals_psf, kp_mags_psf, color = 'blue', marker = 'o')#, label = "PSF")
        plt.step(x=sep_vals_psf, y=kp_mags_psf, marker = 'o')

    # Plot masking data
    sep_vals_masking = np.array([15, 30, 60, 120, 200, 280])
    kp_mags_masking = get_data_masking(starname)
    
    if np.all(kp_mags_masking == None): # If no masking data has been found, do not create a plot
        print(f"No masking plot has been generated for {starname}.")
        return None
    else:
        # ax.plot(sep_vals_masking, kp_mags_masking, marker = 'o', color='red')#, label = "Masking")
        plt.step(x=sep_vals_masking, y=kp_mags_masking, marker = 'o')
        ax.invert_yaxis()
        ax.grid()




if __name__ == "__main__":
    # masking = get_data_masking("HSHJ300")
    # print(masking)

    # print()

    # psf = get_data_psf("AD_2595")
    # print(psf)

    # print()

    # plot_star("HSHJ300")
    # plot_star("JS355")

    # List of star names

    targets = pd.read_excel(r'G:/Shared drives/DouglasGroup/data/Copy of Keck Targets.xlsx', index_col=0)
    try:
        fig, ax = plt.subplots()
        for name, obsdate in targets.iterrows():
            if not pd.isnull(obsdate[0]):
                star = name.replace(" ", "_")
                print(star)
                plot_star(star, ax)
    except KeyboardInterrupt:
        plt.title("\u0394 Kp Magnitude vs. Seperation", fontsize = 12)
        ax.set_xlabel("Seperation Values (mas)")
        ax.set_ylabel("\u0394 Kp Magnitudes")
        ax.invert_yaxis()
        plt.grid()
        plt.show()
        plt.close()
    else:
        plt.title("\u0394 Kp Magnitude vs. Seperation", fontsize = 12)
        ax.set_xlabel("Seperation Values (mas)")
        ax.set_ylabel("\u0394 Kp Magnitudes")
        ax.invert_yaxis()
        plt.grid()
        plt.show()
        plt.close()