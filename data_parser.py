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

def get_data_masking(starname): # Get kp mags for a given star (masking dataset)
    
    # Get file
    keck_masking_dir = os.path.expanduser(r'G:/Shared drives/DouglasGroup/data/keck_masking_detections_praesepe/*')

    filelist = glob.glob(os.path.join(keck_masking_dir, f"{starname}*"))
    
    if len(filelist) == 1: # Proper operation! Only one entry found in masking data :)
        filename = filelist[0]
    elif len(filelist) == 0: # If no star with the input name is found
        print(f"No masking data detected for {starname}.")
        return None
    elif len(filelist) > 1: # If there is more than one star with the input name
        print(f"Easton, we have a problem! Multiple masking entries detected for {starname} :(")
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
        data_as_string = raw_data[raw_data.find(r'&')+1 : raw_data.find(r'\\')] # Skips the "99% only" and runs until the "\\" at the end of the line
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
    starname_psf = f"{starname}" # Names separated due to different name formats between datasets

    if "_" in f"{starname_psf}":
        starname_psf = starname.replace("_", " ")
    elif (starname_psf[0:4] == "EPIC" and starname_psf[0:5] != "EPIC "):
        starname_psf = re.sub("EPIC", "EPIC ", starname_psf)

    # print(f"Names: {starname} (masking), and {starname_psf} (psf)\n")

    names_col = ["Name", "Epoch", "Filter", "N_obs", "t_int", "150", "200", "250", "300", "400", "500", "700", "1000", "1500", "2000", "PI"] # Column names for psf data

    # Get data path
    keck_psf_dir = os.path.expanduser(r'G:/Shared drives/DouglasGroup/data/keck_psf_detections_praesepe/')
    praesepe_data = at.read(os.path.join(keck_psf_dir, 'detlimtable_praesepe.txt'), data_start=0, delimiter="&", names = names_col, fill_values=[("...", np.nan)])
    
    praesepe_new = (praesepe_data["PI"] == "Douglas") & (praesepe_data["Filter"] == "Kp") # Kp is a filter in near infrared
    our_praesepe_data = praesepe_data[praesepe_new]
    
    # Making numpy array of the date (index = 0) and mag values (indeces 1 --> 10)
    our_praesepe_data_new = np.array([our_praesepe_data["Epoch"], our_praesepe_data["150"], our_praesepe_data["200"], our_praesepe_data["250"], our_praesepe_data["300"], our_praesepe_data["400"], our_praesepe_data["500"], our_praesepe_data["700"], our_praesepe_data["1000"], our_praesepe_data["1500"], our_praesepe_data["2000"]])
    
    row_index_list = np.where(our_praesepe_data["Name"]==starname_psf)
    
    # print("Row index list: ", row_index_list[0])
    # print("PSF date:",our_praesepe_data_new[0][5])
    
    if np.size(row_index_list) == 0: # No entry with the input name
        print(f"No PSF data was found for star {starname_psf}.")
        return None
    
    elif np.size(row_index_list) > 1: # More than one entry with the input name
    
        if get_data_masking(starname) is not None: # Check if the input star has masking data
            data_masking = get_data_masking(starname)
            # print("Masking data:", data_masking)
            # print("Masking data date:", data_masking[0])
            # print("PSF data dates:", our_praesepe_data_new[0][row_index_list[0]])
            for entry in range(np.size(row_index_list[0])): # Loop through all entries of the star that has duplicates
                current_entry = our_praesepe_data_new[0][row_index_list[0][entry]] # Date for the current entry being tested
                # print(f"Entry #{entry}:")
                # print("Potential match date:", our_praesepe_data_new[0][row_index_list[0][entry]])
                if abs(data_masking[0] - current_entry) < 0.5: # Check which entry date matches the masking data date. Detects if the dates were taken the same night
                    # print(f"Match found! Entry #{entry}")    
                    # print("Matching date:", current_entry, "\n")
                    kp_mags_psf = np.asarray(our_praesepe_data_new[:, row_index_list[0][entry]].astype(float))
                    kp_mags_psf = np.reshape(kp_mags_psf, len(kp_mags_psf))
                    # print(kp_mags_psf)
                    return kp_mags_psf
                    # break
                else: # Current entry doesn't match, continue with loop
                    # print(f"Entry #{entry} was not observed on the same night as {starname}\n")
                    continue

        else: # No masking data!
            print(f"Easton, we have a problem! Multiple PSF entries detected and no masking data found for {starname} :(")
            print("Row indeces detected:", row_index_list[0])
            i=0 # i is used to differentiate the different entries
            
            # Finding which entry has a higher magnitude contrast (since we could see fainter companions in these cases)
            for entry in range(np.size(row_index_list[0])): # Loop through all entries of the star that has duplicates
                i+=1
                print(i, np.asarray(our_praesepe_data_new[:, row_index_list[0][entry]][1:].astype(float)))
                maximum = np.asarray(our_praesepe_data_new[:, row_index_list[0][entry]].astype(float))
                current = np.asarray(our_praesepe_data_new[:, row_index_list[0][entry]].astype(float))
                if maximum[0] < current[0]:
                    maximum = np.asarray(our_praesepe_data_new[:, row_index_list[0][entry]].astype(float))
            print("Entry with highest mag:", maximum[1:])
            return maximum

    else: # Only one entry for the input name (proper functionality)! :)
        kp_mags_psf = np.asarray(our_praesepe_data_new[:, row_index_list[0]].astype(float))
        kp_mags_psf = np.reshape(kp_mags_psf, len(kp_mags_psf))
        # print("All seems good!")
        return kp_mags_psf

def plot_star(starname, ax=None): # Plot psf and masking curves for a given star
    if ax is None: 
        # Create the plot for the input star
        fig, ax = plt.subplots()
        plt.title(f"\u0394 Kp Magnitude vs. Seperation", fontsize = 12)
        ax.set_xlabel("Seperation Values (mas)")
        ax.set_ylabel("\u0394 Kp Magnitudes")  #\u0394 is a Delta symbol


    # Plot data
    sep_vals_psf = np.array([150, 200, 250, 300, 400, 500, 700, 1000, 1500, 2000])
    sep_vals_masking = np.array([15, 30, 60, 120, 200, 280])

    # Join separation values
    sep_vals_all = np.concatenate((sep_vals_masking[0 : np.size(sep_vals_masking) - 2], sep_vals_psf))
    
    # Make sure there is data from both sources
    if get_data_psf(f"{starname}") is None: # If no psf data has been found, do not create a plot
        print(f"No PSF data found; no plot has been generated for {starname}.\n")
        return None
    elif get_data_masking(f"{starname}") is None: # If no masking data has been found, do not create a plot
        print(f"No masking data found; no plot has been generated for {starname}.\n")
        return None
    else:
        # Join masking data (except last 2 points) and psf data into one array and plot
        # Doing so since psf data is better at the overlapping points
        kp_mags_psf = get_data_psf(starname)[1:]
        kp_mags_masking = get_data_masking(starname)[1:]

        
        kp_mags_all = np.concatenate((kp_mags_masking[0 : np.size(kp_mags_masking) - 2], kp_mags_psf))
        plt.step(x=sep_vals_all, y=kp_mags_all, color = "blue", alpha = 0.75)

def export_star(starname): # Export the data to a file that MOLUSC can take. Lots of overlap with the plot_star function

    # Get separation values
    sep_vals_psf = np.array([150, 200, 250, 300, 400, 500, 700, 1000, 1500, 2000])
    sep_vals_masking = np.array([15, 30, 60, 120, 200, 280])
    sep_vals_all = np.concatenate((sep_vals_masking[0 : np.size(sep_vals_masking) - 2], sep_vals_psf))
    
    # Make sure there is data from both sources
    if (get_data_psf(f"{starname}") is None) and (get_data_masking(f"{starname}") is None): # If no data has been found, do not try to export
        print(f"No data found; no file has been generated for {starname}.\n")
        return None
    elif (get_data_psf(f"{starname}") is not None) and (get_data_masking(f"{starname}") is not None): # Have both sets of data
        # Join masking data (except last 2 points) and psf data into one array and plot
        # Doing so since psf data is better at the overlapping points
        kp_mags_psf = get_data_psf(starname)[1:]
        kp_mags_masking = get_data_masking(starname)[1:]
        kp_mags_all = np.concatenate((kp_mags_masking[0 : np.size(kp_mags_masking) - 2], kp_mags_psf))
        
        # Create a table from the data
        data = Table(data = [sep_vals_all, kp_mags_all], names = ["Sep", "Contrast"])
        
        # Export the data to a file, so long as there is not already a file for the star
        export_dir = os.path.expanduser(r'G:/Shared drives/DouglasGroup/Jared Sofair 2022/MOLUSC Inputs/Data Parser Tables')

        filelist = glob.glob(os.path.join(export_dir, f"{starname}.txt"))
        
        if len(filelist) >= 1: # File found!
            print(f"Error: file detected for {starname}\n")
            return None
        elif len(filelist) == 0: # No file has been previously generated
            data.write(os.path.join(export_dir, f"{starname}.txt"), format = 'ascii.basic', delimiter = ' ', overwrite=True)
            print(f"File generated for {starname}!\n")
            # return data
            
    elif (get_data_psf(f"{starname}") is None) and (get_data_masking(f"{starname}") is not None): # Only masking
        kp_mags_masking = get_data_masking(starname)[1:]
        data = Table(data = [sep_vals_masking, kp_mags_masking], names = ["Sep", "Contrast"])
        
        export_dir = os.path.expanduser(r'G:/Shared drives/DouglasGroup/Jared Sofair 2022/MOLUSC Inputs/Data Parser Tables')
        
        filelist = glob.glob(os.path.join(export_dir, f"{starname}.txt"))

        if len(filelist) >= 1: # File found!
            print(f"Error: file detected for {starname}\n")
            return None
        elif len(filelist) == 0: # No file has been previously generated
            data.write(os.path.join(export_dir, f"{starname}.txt"), format = 'ascii.basic', delimiter = ' ', overwrite=True)
            print(f"File generated for {starname}!\n")

    elif (get_data_psf(f"{starname}") is not None) and (get_data_masking(f"{starname}") is None): # Only psf
        kp_mags_psf = get_data_psf(starname)[1:]
        data = Table(data = [sep_vals_psf, kp_mags_psf], names = ["Sep", "Contrast"])
        
        export_dir = os.path.expanduser(r'G:/Shared drives/DouglasGroup/Jared Sofair 2022/MOLUSC Inputs/Data Parser Tables')
        
        filelist = glob.glob(os.path.join(export_dir, f"{starname}.txt"))
   
        if len(filelist) >= 1: # File found!
            print(f"Error: file detected for {starname}\n")
            return None
        elif len(filelist) == 0: # No file has been previously generated
            data.write(os.path.join(export_dir, f"{starname}.txt"), format = 'ascii.basic', delimiter = ' ', overwrite=True)
            print(f"File generated for {starname}!\n")


if __name__ == "__main__":
    # JS355 is the "control" -- it tends to work without problems
    # AD_0738 is experimental -- it's name changes between data sets due to the difference in formatting between them
    # EPIC211998192 is experimental -- it's name format also changes between data sets
    # HSHJ300 is experimental -- it has a masking file but has no masking data; it has PSF data for Kp and J
    
    
    # masking = get_data_masking("JS355")
    # print(masking)

    # masking = get_data_masking("AD_0738")
    # print(masking)

    # masking = get_data_masking("EPIC211998192")
    # print(masking)

    # masking = get_data_masking("HSHJ300")
    # print(masking)

    # print()

    # psf = get_data_psf("JS355")
    # print(psf)

    # psf = get_data_psf("AD_0738")
    # print(psf)

    # psf = get_data_psf("EPIC211998192")
    # print(psf)

    # psf = get_data_psf("HSHJ300")
    # print(psf)

    # print()

    # plot_star("JS355")
    # plot_star("AD_0738")
    # plot_star("EPIC211998192")
    # plot_star("HSHJ300")

    # print()
    
    # export_star("JS355")
    # print(export)

    # print()
 
    # Plot all targets:
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


    # Export all targets to its own file that MOLUSC can take
    list_plot = []
    list_exp_masking = []
    list_exp_psf = []
    
    targets = pd.read_excel(r'G:/Shared drives/DouglasGroup/data/Copy of Keck Targets.xlsx', index_col=0)
    for name, obsdate in targets.iterrows():
        print("Name:", name, "Obsdate:", obsdate[14])
        if not pd.isnull(obsdate[14]):
            
            if get_data_masking(name) is not None: # Has masking data
                export_star(name)
                star = name.replace(" ", "_")
                print(star, "\n")
                list_exp_masking.append(star)
            else: # No masking data :(
                star = name.replace(" ", "_")
                print(f"No masking data for {name}")
                print(star, "\n")
            
            if get_data_psf(name) is not None: # Has psf data
                export_star(name)  
                star = name.replace(" ", "_")
                print(star, "\n")
                list_exp_psf.append(star)
            else: # No psf data :(
                star = name.replace(" ", "_")
                print(f" No PSF data for {name}")    
                print(star, "\n")
        else:
            print(f"File not generated for {name} (no observation date)\n")

    # print("# of plots:", len(list_plot))
    print("# of exports:", len(list_exp_psf), len(list_exp_masking))


# To do next:
# Output data with J filter as well