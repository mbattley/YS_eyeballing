#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 14:37:45 2020

eyeball_from_file.py
--------------------
This script allows you to replot any of the eyeballing plots from scratch using 
saved data for original, deterended, lowess-fitted lcs and periododram period/powers

Runs without inputs, but requires data files in reachable folders

Basic op

Please use the below codes to aid subsequent analysis, but feel free to also write additional notes

Main eyeballing codes:
EA: Algol-type Eclipsing Binary 
EB: Beta-Lyrae Eclipsing Binary
EW: W Ursae Majoris type Eclipsing Binary
VAR: Other stellar varibility
DIS: Used to put target forward for further discussion - used for good transit-like signals that are less confident than TRA signals
TRA: Clear transit-like signal which needs follow-up
FLARE: Signal is caused by flare
BKG: Signal is caused by background
Scatter: Light-curve evidencing high scatter 
Outlier: Signal caused by one to a few outliers


@author: mbattley
"""

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import csv
#from matplotlib import backend_bases
from astropy.constants import R_sun, R_earth, R_jup
#from astropy.table import Table

# print(matplotlib.get_backend())
matplotlib.use("Qt5Agg")

######################## Set font sizes #######################################
SMALL_SIZE = 8
MEDIUM_SIZE = 12
BIGGER_SIZE = 16

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

###############################################################################
sector = 14
YSL_filename = '/Users/matthewbattley/Documents/Young_Star_lists/For_search/YSL_20221028_S{}_bright.csv'.format(sector)
# save_path = '/Users/matthewbattley/Documents/YS_Surveys/QLP/S{}/S{}/'.format(sector,sector)
save_path = '/Users/matthewbattley/Documents/YS_Surveys/QLP/S{}/Cleaned/'.format(sector)
# info_filename = save_path + 'Period_info_table_S{}.csv'.format(sector)
info_filename = save_path + 'Period_info_table_S{}_Cleaned.csv'.format(sector)

info_table = pd.read_csv(info_filename,skiprows =1)
target_list = np.array(info_table['TIC'])
stellar_info_table = pd.read_csv(YSL_filename)
# test_target_list = [410214986,410214986,410214986]


##### FOR TESTING ######
# filename = '/Users/matthewbattley/Documents/Young_Star_Lists/For_search/YSL_20221028_S14_bright.csv'
# data = Table.read(filename, format='ascii.csv')
# target_list = data['TIC'][0:3]
# start = 0
# end = 3
########################

j = 0
i = 0
start = 310 #For uncleaned: 222
end = 400 #For uncleaned: 250
i = start

for tic in target_list[start:end]:
    ############################ Load Data ####################################
    print(tic)
    period = info_table['Max Period'][i]
    epoch = info_table['Epoch of Max'][i]
    depth = info_table['Depth of max'][i]
    p_rot = info_table['Rotation Period'][i]
    p2 = info_table['2nd highest period'][i]
    p3 = info_table['3rd Highest Period'][i]
    lc_filename = save_path + '{}_lc_data.csv'.format(tic)
    periodogram_filename = save_path + '{}_periodogram_data.csv'.format(tic)
    
    # lc_data = pd.read_csv('data.csv', delimiter='|', names=list(range(7)))
    # OR
    lc_data     = pd.read_csv(lc_filename, header=None,skiprows =1)
    lc_data = lc_data.T
    time_og     = lc_data[0]
    flux_og     = lc_data[1]
    flux_bkg    = lc_data[2]
    time_f      = lc_data[3]
    flux_detr   = lc_data[4]
    lowess_time = lc_data[5]
    lowess_flux = lc_data[6]
    phase = np.mod(time_f-epoch-period/2,period)/period 
    
    # Open periodogram data for lightcurve
    periodogram_data = pd.read_csv(periodogram_filename,header=None).T
    BLS_periods = periodogram_data[0]
    BLS_power = periodogram_data[1]
    
    # Get stellar info
    k = list(stellar_info_table['TIC']).index(tic)
    r_star = stellar_info_table['Rad'][k]
    r_p_m = r_star*np.sqrt(depth)*R_sun # In meters
    r_p_Jup = r_p_m/R_jup           # Radius in Jupiter radii
    r_p_Earth = r_p_m/R_earth       # Radius in Earth radii
    print('R_p = {} R_Jup'.format(r_p_Jup))
    print('R_p = {} R_Earth'.format(r_p_Earth))
    
    ############################ Flux Compar Fig ##############################
    flux_compar_fig = plt.figure()
    mngr = plt.get_current_fig_manager()
    mngr.window.setGeometry(0,0,960, 700)
    ax1_fl = plt.subplot(311)
    ax2_fl = plt.subplot(312, sharex=ax1_fl)
    ax3_fl = plt.subplot(313, sharex=ax1_fl)
    
    ax1_fl.scatter(time_og, flux_og,  s=1, c='k', label='Original Flux')
    ax2_fl.scatter(time_f, flux_detr, s=1, c='g', label='Detrended Flux')
    ax3_fl.scatter(time_og,flux_bkg,  s=1, c='r', label='Background Flux')
    ax3_fl.set_xlabel('Time - 2457000 [BTJD days]')
    #plt.legend()
    plt.subplots_adjust(wspace=0.0)
    flux_compar_fig.tight_layout(h_pad =-1.5)
    ax1_fl.legend(loc='upper right')
    ax2_fl.legend(loc='upper right')
    ax3_fl.legend(loc='upper right')

    
    ############################ Phase Folded Fig #############################
    test_epoch = 1332
    phase_rot = np.mod(time_f-test_epoch-p_rot/2,p_rot)/p_rot
    phase_p2 = np.mod(time_f-epoch-p2/2,p2)/p2
    phase_p3 = np.mod(time_f-epoch-p3/2,p3)/p3
    
    phase_fold_fig  = plt.figure()
    mngr = plt.get_current_fig_manager()
    mngr.window.setGeometry(960,0,960, 700)
    ax1_p = plt.subplot(311)
    ax2_p = plt.subplot(312, sharex=ax1_p)
    ax3_p = plt.subplot(313, sharex=ax1_p)
    
    ax1_p.scatter(phase_rot, flux_detr,  s=1, c='k', label='p_rot = {:0.4} days'.format(p_rot))
    ax2_p.scatter(phase_p2, flux_detr, s=1, c='g', label='p2 = {:0.4} days'.format(p2))
    ax3_p.scatter(phase_p3,flux_detr,  s=1, c='r', label='p3 = {:0.4} days'.format(p3))
    ax3_p.set_xlabel('Phase')
    ax1_p.set_xlim(0,1)
    plt.subplots_adjust(wspace=0.0)
    phase_fold_fig.tight_layout(h_pad =-1.5) #previously h_pad = -1.0
    ax1_p.legend(loc='upper right')
    ax2_p.legend(loc='upper right')
    ax3_p.legend(loc='upper right')

    ########################### Main Eyeballing Fig ###########################
    eyeballing_fig = plt.figure(figsize = (16,7),dpi = 120)
    
    ax1 = plt.subplot(221)
    ax2 = plt.subplot(222, sharex=ax1)
    ax3 = plt.subplot(223)
    ax4 = plt.subplot(224)
    
    # Original DIA with injected transits setup
    ax1.scatter(time_og, flux_og, s=1, c= 'k')
    ax1.set_ylabel('Normalized Flux')
    ax1.set_xlabel('Time- 2457000 [BTJD days]')
    ax1.plot(lowess_time,lowess_flux)
    
    # Detrended figure setup
    ax2.scatter(time_f, flux_detr, c = 'k', s = 1, label = 'TIC-{} residuals after lowess detrending'.format(tic))
    ax2.set_ylabel('Normalized Flux')
    ax2.set_xlabel('Time - 2457000 [BTJD days]')
    ax2.text(0.65,0.90,'R_star = {:0.4} R_sun'.format(r_star),transform=ax2.transAxes)
    
    # Periodogram setup
    ax3.plot(BLS_periods, BLS_power, "k", lw=0.5)
    ax3.set_xlim(np.min(BLS_periods), np.max(BLS_periods))
    ax3.set_xlabel("period [days]")
    ax3.set_ylabel("log likelihood")
    ax3.axvline(period, alpha=0.4, lw=3)
    for n in range(2, 10):
        ax3.axvline(n*period, alpha=0.4, lw=1, linestyle="dashed")
        ax3.axvline(period/n, alpha=0.4, lw=1, linestyle="dashed")
    
    # Folded or zoomed plot setup
    ax4.scatter(phase,flux_detr,c='k',s=1)
    ax4.text(0.75,0.90,'P = {:0.4} days'.format(period),transform=ax4.transAxes)
    ax4.set_xlabel('Phase')
    ax4.set_ylabel('Normalized Flux')
    ax4.set_xlim(0,1)
    
    eyeballing_fig.tight_layout()
    
    mngr = plt.get_current_fig_manager()
    # geom = mngr.window.geometry()
    # x,y,dx,dy = geom.getRect()
    # print('({},{},{},{})'.format(x,y,dx,dy))
    mngr.window.setGeometry(0,1200,1920, 884)
    
    print('Type eyeballing note, then click on plot and press enter')
    
    plt.show(block=False)
    pts = plt.ginput(50, timeout=-1)
    # plt.waitforbuttonpress()
    # eyeballing_fig.canvas.mpl_connect('key_press_event', press)
    note = input('<Enter Eyeballing Notes>')
    print(note)
    
    plt.close('all')
    
    j+=1
    i+=1
    print('Number of lcs eyeballed = {}'.format(j))
    
    with open(save_path+'S{}_eyeballing_notes.csv'.format(sector),'a') as f:
        data_row = [tic,note]
        writer = csv.writer(f, delimiter=',')
        writer.writerow(data_row)

    