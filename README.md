# YS_eyeballing
Eyeballing code for use with YSD pipeline. If you use this code for your work (or the YSD pipeline), please cite Battley, Pollacco & Armstrong (2020).

Script to generate eyeballing plots for a set of previously detrended light-curves detrended using the YSD pipeline. Also allows for eyeballing notes to be collated on a target by target basis.

This script can be run from either a local machine or external device, depending where the main detrended light-curevs are situated.

Required Input files:
- Period_info_table_[sector]_Cleaned.csv -> CSV file containing the main period info table generated by the YSD pipeline, including the main periods of both the stellar variability and BLS search.
- Young_Star_list.csv -> Copy of the main young star list compiled by Dr Matthew Battley and Dr Ed Gillen 
- [Target]_lc_data.csv -> Light-curve data per object, including time, original lc, background flux, detrended lc and model lc
- [Target]_periodogram_data.csv -> Data from the BLS search per object, to plot the periodogram plot
Note that the path to each of these can be adjusted within the main script

Generates three figures:
1. Comparision between original, detrended and background fluxes
2. Detrended flux folded by the rotation period and 2nd/3rd strongest periods from the BLS search
3. Main eyeballing plot, four panels:
    a) Original lightcurve with overplotted model used to detrend
    b) Detrende lightcurve
    c) BLS periodogram from detrended lc
    d) Detrended lc folded by main BLS period
    
Running the code:

The code should be run from the command line, using the command 'python3 eyeball_from_file.py'

When run, the script will plot the three plots for each object in turn and wait for user input for each target before moving on. To move through the objects one by one first write your eyeballing notes (or leave blank), then click on the main eyeballing plot and press enter to move to the next target.
