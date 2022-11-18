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
import mmap

csv_path_github = os.path.expanduser(r"C:/Users/Jared/Documents/GitHub/data-parser/CSV Files")
csv_path_drive = os.path.expanduser(r"G:/Shared drives/DouglasGroup/data/WIYN")
rv_output_path = os.path.expanduser(r"G:/Shared drives/DouglasGroup/Jared Sofair 2022/MOLUSC/Data Parser tables/MOLUSC RVs")

targets_abr = Table.read(os.path.join(csv_path_github, r'targets_abr.csv'))
wiyn = Table.read(os.path.join(csv_path_drive, r'WIYN_RVs_matchedKeck.csv'))
targets_wiyn = wiyn['Name']  

# Get RV, RVerr, and HJD for each target that we have multiple RVs for
def output_rv_data(starname):
    # Row indeces where the target shows up in wiyn
    target_indeces = np.where(targets_wiyn == starname)[0]
    # print(targets_wiyn[np.where(targets_wiyn == starname)])

    # Grab RV, RVerr, and HJD at each index
    rv = np.array([])
    rverr = np.array([])
    jd = np.array([])
    
    for i in target_indeces:
        rv = np.append(rv, wiyn["rv"][i])
        rverr = np.append(rverr, wiyn["err"][i])
        
        # Also converts HJD to JD using MJD to JD conversion (since we don't need to be too precise :))
        jd = np.append(jd, wiyn["hjd"][i]+2400000.5)
    
    # print(rv)
    # print(rverr)
    # print(jd)
    
    # Export rv data into a text file that MOLUSC can read!
    if starname.find("_") != -1:
        starname = starname.replace(" ", "_")
    rv_data = Table(data = [jd, rv, rverr], names = ["JD", "RV", "RVerr"])
    rv_data.write(os.path.join(rv_output_path, f"{starname}.txt"), format = 'ascii.basic', delimiter = ' ', overwrite=True)

# Is a given target a binary based on the masking data results?
def masking_binary(starname):
    if " " in f"{starname}":
        starname = starname.replace(" ", "_")
    #     print(starname)
    # else:
    #     print(starname)
    
    # Get file
    keck_masking_dir = os.path.expanduser(r'G:/Shared drives/DouglasGroup/data/keck_masking_detections_praesepe/*')

    filelist = glob.glob(os.path.join(keck_masking_dir, f"{starname}*"))

     # Proper operation! Only one entry found in masking data :)
    if len(filelist) == 1:
        filename = filelist[0]
     # If no star with the input name is found
    elif len(filelist) == 0: 
        # print(f"No masking data detected for {starname}.")
        return None
     # If there is more than one star with the input name
    elif len(filelist) > 1:
        print(f"Easton, we have a problem! Multiple masking entries detected for {starname} :(")
        print(filelist)
        return None
    
    with open(filename, 'rb', 0) as file:
        # I actually don't know if I need to do vvv
        # TODO: Get the length of the file - maybe can use mmap features to make more efficient?
        # length = len(file.readlines())

        # Open the file as a memory mapped file (for efficiency)
        f = mmap.mmap(file.fileno(), 0, access=mmap.ACCESS_READ)
        if f.find(b'Significance') != -1:
            # Get line number of the Significance
            sig_ind = f.find(b'Significance')
            
            # TODO: Seek to and get Significance value [WIP]
            f.seek(sig_ind)
            sig_line = f.readline().decode('utf-8').replace(' ', '')
            
            # Significance value happens after the : and ends before the =
            a = sig_line.index(':')+1
            b = sig_line.index('=')
            sig = sig_line[a:b]
            if float(sig) > 10:
                print(f"{starname}: {sig}")
        # else:
        #     print("No significance data found!")

# with open(STAT_FILE, "r+b") as f:
#     map_file = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)
#     for line in iter(map_file.readline, b""):
        # whatever


    # Turn the data into an array
    # try:
    #     # raw_data = lc.getline(filename)
    #     print(raw_data)
    #     data_as_string = raw_data[raw_data.find(r'Significance') : raw_data.find(r')')]
    #     data_as_arr = data_as_string.split('&')
    #     data = np.asarray(data_as_arr, dtype="float")
    #     print(data)
    # except ValueError:
    #     print(f"No significance was found for the star {starname}.")
    #     return None

    # Create the table from the data
    # names_col_masking = ["Epoch", "15", "30", "60", "120", "200", "280"] # Column names for data from keck_masking_detections_praesepe ICE
    
    # sep_vals = np.asarray(t_masking.colnames[1:],dtype="float") (have here ICE)
    # kp_mags_masking = data[1:]
    return
    # return data

def get_data_masking(starname): # Get kp mags for a given star (masking dataset)
    if " " in f"{starname}":
        starname = starname.replace(" ", "_")
        
    # Get file
    keck_masking_dir = os.path.expanduser(r'G:/Shared drives/DouglasGroup/data/keck_masking_detections_praesepe/*')

    filelist = glob.glob(os.path.join(keck_masking_dir, f"{starname}*"))
    
    if len(filelist) == 1: # Proper operation! Only one entry found in masking data :)
        filename = filelist[0]
    elif len(filelist) == 0: # If no star with the input name is found
        # print(f"No masking data detected for {starname}.")
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
        plt.step(x=sep_vals_all, y=kp_mags_all, color = "#a10f05", alpha = 0.37)

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
        export_dir = os.path.expanduser(r'G:/Shared drives/DouglasGroup/Jared Sofair 2022/MOLUSC/Data Parser Tables')
        starname = starname.replace(" ", "_")
        filelist = glob.glob(os.path.join(export_dir, f"{starname}.txt"))
        
        if len(filelist) >= 1: # File found!
            print(f"{starname} file generation skipped; file detected\n")
            return None
        elif len(filelist) == 0: # No file has been previously generated
            data.write(os.path.join(export_dir, f"{starname}.txt"), format = 'ascii.basic', delimiter = ' ', overwrite=True)
            print(f"File generated for {starname}!\n")
            # return data
            
    elif (get_data_psf(f"{starname}") is None) and (get_data_masking(f"{starname}") is not None): # Only masking
        kp_mags_masking = get_data_masking(starname)[1:]
        data = Table(data = [sep_vals_masking, kp_mags_masking], names = ["Sep", "Contrast"])
        
        export_dir = os.path.expanduser(r'G:/Shared drives/DouglasGroup/Jared Sofair 2022/MOLUSC/Data Parser Tables')
        starname = starname.replace(" ", "_")
        filelist = glob.glob(os.path.join(export_dir, f"{starname}.txt"))

        if len(filelist) >= 1: # File found!
            print(f"{starname} file generation skipped; file detected\n")
            return None
        elif len(filelist) == 0: # No file has been previously generated
            data.write(os.path.join(export_dir, f"{starname}.txt"), format = 'ascii.basic', delimiter = ' ', overwrite=True)
            print(f"File generated for {starname}!\n")

    elif (get_data_psf(f"{starname}") is not None) and (get_data_masking(f"{starname}") is None): # Only psf
        kp_mags_psf = get_data_psf(starname)[1:]
        data = Table(data = [sep_vals_psf, kp_mags_psf], names = ["Sep", "Contrast"])
        
        export_dir = os.path.expanduser(r'G:/Shared drives/DouglasGroup/Jared Sofair 2022/MOLUSC/Data Parser Tables')
        starname = starname.replace(" ", "_")
        filelist = glob.glob(os.path.join(export_dir, f"{starname}.txt"))
   
        if len(filelist) >= 1: # File found!
            print(f"{starname} file generation skipped; file detected\n")
            return None
        elif len(filelist) == 0: # No file has been previously generated
            data.write(os.path.join(export_dir, f"{starname}.txt"), format = 'ascii.basic', delimiter = ' ', overwrite=True)
            print(f"File generated for {starname}!\n")

if __name__ == "__main__":    
    print(":)")
#%% Basic function testing
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

    # masking = get_data_masking("HIP14807")
    # print(f"hip masking: {masking}")


    # print()

    # psf = get_data_psf("JS355")
    # print(psf)

    # psf = get_data_psf("AD_0738")
    # print(psf)

    # psf = get_data_psf("EPIC211998192")
    # print(psf)

    # psf = get_data_psf("HSHJ300")
    # print(psf)

    # psf = get_data_psf("HIP14807")
    # print(f"hip psf: {psf}")

    # print()

    # plot_star("JS355")
    # plot_star("AD_0738")
    # plot_star("EPIC211998192")
    # plot_star("HSHJ300")

    # print()
    
    # export_star("JS355")
    # print(export)

    # print()
    
    # masking_binary("JS230")
    # masking_binary("AD_3663")
        
    # Get a list of confirmed binaries from masking data
    # targets_abr = Table.read(os.path.join(csv_path_github, r'targets_abr.csv'))
    # for name in targets_abr.iterrows('name'):
    #     # get_data_masking(name[0])
    #     # plot_star(name[0])
    #     masking_binary(name[0])
        # print(name[0])
    
    
    # print(targets_abr)
    

#%% Get a list of targets we have RVs for (they come up multiple times in WIYN_RVs_matchedKeck.csv)
    
    # List of all targets which appear multiple times in wiyn
    # targets_rv = np.array([])
    
    # for i, name in enumerate(targets_abr.iterrows('name')):
    #     # Count the number of times a target appears in wiyn
    #     count = np.count_nonzero(targets_wiyn == targets_abr['name'][i])
    #     # print(f"{name[0]}: {count}")
        
    #     # If it is in wiyn multiple times, put it in targets_rv
    #     if count > 1:
    #         targets_rv = np.append(targets_rv, name[0])
                        
    # print(targets_rv)
#%% Use the list of targets we have RVs for to generate files that we can feed MOLUSC
# for i, name in enumerate(targets_rv):
#     output_rv_data(name)
#     print(i)
    
#%% TODO: Take targets_abr and add rv, rverr, hjd (converted to mjd) and output them to txt file readable by MOLUSC

        
#%% To do:
# Output data with J filter as well