# data-parser [WIP]

## 1. keck_parser.py (formerly data_parser.py) [WIP]

## 2. praesepe_cmd_mass.py [WIP]

## 3. bhac_gaia_2mass_reader.py [WIP]

## 4. module_runner.py [WIP]
Collection of scripts and scrap work used at some point in the lifetime of this repo.
### Former keck_parser.py actions
Scripts previously in keck_parser.py to:
  plot all of our targets; 
  export the contrast and separation to be read by MOLUSC; 
  and get the Keck, Gaia DR2, and 2MASS designations for our targets;
### Get absG mag and BP-RP from praesepe_merged
Script to get absolute G magnitude and BP-RP color from praesepe_merged (pm) table. Done so by:
  getting pm indeces for target stars,
  converting apparent G magnitude to absolute G magnitude,
  getting apparent K magnitude from pm,
  and converting apparent K magnitude to absolute K magnitude.
For converting from apparent magnitude to absolute magnitude, the distances given in pm were used.

## 5. CSV Files [WIP]
### praesepe_merged.csv

### targets_abr.csv
Collection of data from our targets:
  Keck names (from Keck_Targets_2018B_Simbad.csv),
  2MASS designation (from pm),
  Gaia DR2 designation (from pm),
  BP-RP color index (calculated using BP and RP columns from pm),
  absolute G magnitudes (calculated using apparent G magnitude and distance from pm),
  absolute K magnitudes (calculated using apparent K magnitude and distance from pm),
  and pm index.
The data in each column was compiled and added to this file using various scripts in module_runner.py.

## 6. Plots [WIP]
