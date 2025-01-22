# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 12:04:02 2025

@author: OllieTwigge
"""
# this asks for the directory with a binary and cache folders, sensorlist.txt and deploymentyaml"
# then it asks which dbd file you want to process and spits out all the graphs from that yo
#%% Choose working library - (and setup+definitions)
### Choose working library - (and setup+definitions)
def checkCACfiles(binarydir, cacdir):
    print("Checking CACHE files...")
    missingcacs = []
    # make a list of all the availbe cac files
    caclst = [os.path.splitext(x)[0] for x in os.listdir(cacdir)]
    for file in os.listdir(binarydir):
        if file.endswith(('.ebd', '.dbd', '.tbd', '.sbd')):
            # print(os.path.join(binarydir, file))
            binfilepath = os.path.join(binarydir, file)
            with open(binfilepath, 'rb') as binary_file:
                for line in binary_file:
                    if line.startswith(b'sensor_list_crc'):
                        # print(f"Match found: {line}")
                        cacfile = str(line).split(':')[1].strip()[:8]
                        print(file, cacfile)
                        if cacfile in caclst:
                            continue
                        else:
                            print(file, "MISSING CACHE FILE", cacfile)
                            missingcacs.append(cacfile)
    missingcacs = list(set(missingcacs))
    if not missingcacs:
        print("ALL CACHE FILES AVAILABLE")
    else:
        print('MISSING CAC FILES:', missingcacs)
        
def processGliderlons(lons):
    lon_dd = []
    for val in lons:
        if np.isnan(val):  # Preserve NaN
            lon_dd.append(np.nan)
        else:
            # Handle negative values properly
            sign = -1 if val < 0 else 1
            abs_val = abs(val)
            # Extract degrees and minutes
            # Integer part before the decimal is degrees
            degrees = int(abs_val / 100)
            minutes = abs_val % 100      # Remainder is minutes
            lon_dd.append(sign * (degrees + minutes / 60))
    return lon_dd

def processGliderlats(lats):
    lat_dd = []
    for val in lats:
        if np.isnan(val):  # Preserve NaN
            lat_dd.append(np.nan)
        else:
            # Handle negative values properly
            sign = -1 if val < 0 else 1
            abs_val = abs(val)
            # Extract degrees and minutes
            # Integer part before the decimal is degrees
            degrees = int(abs_val / 100)
            minutes = abs_val % 100      # Remainder is minutes
            lat_dd.append(sign * (degrees + minutes / 60))
    return lat_dd

def average_science(divsds, bin_size = 1):
    depth_var = divsds['sci_rbrctd_depth_00']
    depth_min = depth_var.min().item()
    depth_max = depth_var.max().item()
    print(f"Depth range: {depth_min:.2f} to {depth_max:.2f} meters")
    bin_size = bin_size  # 1-meter bins
    depth_bins = np.arange(np.floor(depth_min), np.ceil(depth_max) + bin_size, bin_size)
    depth_bin_centers = (depth_bins[:-1] + depth_bins[1:]) / 2
    print(f'Number of bins: {len(depth_bins)}')
    
    bsds = divsds.groupby_bins(group= divsds['sci_rbrctd_depth_00'], 
                                       bins = depth_bins, 
                                       labels = depth_bin_centers, 
                                          ).mean()
    bsds = bsds.rename_dims({'sci_rbrctd_depth_00_bins': 'depth_bins'})
    bsds = bsds.rename_vars({'sci_rbrctd_depth_00_bins': 'depth_bins'})
    
    return bsds
def rndup(n, base = 50):
    return n + (base - n) % base

def rnddown(n, base = 50):
    return int(n // base)*base

def autorange(series, increment):
    return rnddown(series.min(), increment),rndup(series.max(), increment)

def autorange_percentiles(series):
    return (np.percentile(series, 5), np.percentile(series, 95))

def autominmax(series):
    return (np.nanmin(series), np.nanmax(series))

#TODO be able to switch between the e/d(big) and s/t(small) bd files 
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import logging
import os
import xarray
import pandas as pd
import pyglider.ncprocess as ncprocess
import pyglider.slocum as slocum
import pyglider.utils as pgutils
from tkinter.filedialog import askopenfilename, askdirectory
import tkinter as tk
from tkinter import filedialog
from datetime import timedelta

def select_directory():
    # Create a visible root window (used as the parent)
    root = tk.Tk()
    root.attributes('-topmost', True)  # Ensure the root window is on top
    root.update()  # Update the root window to apply the attributes

    # Show the directory selection dialog with the root as the parent
    directory = filedialog.askdirectory(title="Select directory with binary, cache, sensorlist and deploymentyaml",
                                        parent=root)

    root.destroy()  # Destroy the root window
    return directory

def select_file():
    # Create a visible root window (used as the parent)
    root = tk.Tk()
    root.attributes('-topmost', True)  # Ensure the root window is on top
    root.withdraw()  # Hide the root window
    
    # Show the file selection dialog with the root as the parent
    file_path = filedialog.askopenfilename(title='Select dive to plot:', #askopenfilename
                            initialdir=rawdir, 
                            filetypes=[('Flight data', '*.dbd.nc'),
                                       ('All NetCDFs', '*.nc'),
                                       ('All Files', '*.*')],
                            parent=root)

    root.destroy()  # Destroy the root window
    return file_path

maindir = select_directory()
# binfile = askopenfilename(title = "Choose the binary file")
os.chdir(maindir)
logging.basicConfig(level='INFO')

binarydir = './raw_dbd_ebd/'
rawdir = './rawnc/'
cacdir = './cac/'
sensorlist = './unit_1104_sensors.txt'
deploymentyaml = './deploymentRealtime.yml'
#l1tsdir = './L0-timeseries/'
#profiledir = './L0-profiles/'
#griddir = './L0-gridfiles/'
scisuffix = 'ebd'  # 'tbd'
glidersuffix = 'dbd'  # 'sbd'
graphs = './graphs/'


checkCACfiles(binarydir, cacdir)

slocum.binary_to_rawnc(
    binarydir, rawdir, cacdir, sensorlist, deploymentyaml,
    incremental=True, scisuffix=scisuffix, glidersuffix=glidersuffix)

if not os.path.exists(graphs):
        os.makedirs(graphs)
        print('CREATED graph folder')
else:
        print('Graphs folder exists')

#%% Plotting individual dives
# Plotting individual dives
# similar to this G:\.shortcut-targets-by-id\1Z7PX__ZLrJ6oL2wHtDtiATh5eAh4Cgg3\Oceanly Science\Instrumentation\Ocean glider\Deployment Tools\example figs and html\betty-2019-263-1-0
# Choosing the file to plot up
divfd = select_file()
divename = os.path.splitext(os.path.splitext(os.path.basename(divfd))[0])[0]

divsd = os.path.splitext(os.path.splitext(os.path.basename(divfd))[0])[0]
divsd = maindir+rawdir+divsd+'.ebd.nc'

divfds = xr.open_dataset(divfd)
divfds = divfds.set_coords('time')
divfds = divfds.swap_dims({'_ind': 'time'})
divfds['time'] = pd.to_datetime(divfds.time.values, unit='s')

try:
    divsds = xr.open_dataset(divsd)
    divsds = divsds.set_coords('time')
    divsds = divsds.swap_dims({'_ind': 'time'})
    divsds['time'] = pd.to_datetime(divsds.time.values, unit='s')
    bsds = average_science(divsds, bin_size= 1)
    skipscience = False
except:
    print(f'SKIPPING {divename} - cannot open or find EBD')
    skipscience = True

fullfn = divfds.attrs['full_filename']
print(f'Dive chosen: {divename}')
divegraphdir = graphs+divename+'/'
if not os.path.exists(graphs+divename):
        os.makedirs(graphs+divename)
        print(f'CREATED dive graph folder {divegraphdir}')
else:
        print('dive graph folder exists')
        
#%% Graphing
#ADD META TABLE TO PRINT
if not skipscience:
    scisens = list(divsds.keys())
    scipoints = len(divsds.time)
else:
    scisens = 'No Science Data'
mdep = f'{np.nanmax(divfds.m_depth.values):.2f}m'
divemet= {
    'First Timestamp': str(pd.to_datetime(divfds.time.values[0])).split('.')[0], 
    'Last Timestamp': str(pd.to_datetime(divfds.time.values[-1])).split('.')[0],
    'Total segment time': str(timedelta(seconds = (divfds.time.values[-1]-divfds.time.values[0])/np.timedelta64(1, 's'))).split('.')[0],
    'Total data points': len(divfds.time),
    'Max Depth': mdep,
    'Number of yos in segment':((np.nanmax(divfds.m_num_half_yos_in_segment.values)+1)/2).astype(int), 
    'Average yo time':	f'{(np.nanmean(divfds.m_avg_yo_time.values)/60):.2f} mins',
    'Average speed': f'{(np.nanmean(divfds.m_avg_speed)):.3f} m/s',
    'Average dive rate': f'{(np.nanmean(divfds.m_avg_dive_rate.values)):.3f} m/s',
    'Average climb rate': f'{(np.nanmean(divfds.m_avg_climb_rate.values)):.3f} m/s',
    'Average downward inflection time':	f'{(np.nanmean(divfds.m_avg_downward_inflection_time.values)).astype(int)}s',
    'Science sensors': scisens,
    'Total science data points': scipoints
    }
plt.ioff() #turns off plotting
#plt.ion()
### Tracks Lat Lon figure
# add start and finish times to the points on the graph

lons = divfds.m_gps_lon.values
lon_gps = processGliderlons(lons)
lats = divfds.m_gps_lat.values
lat_gps = processGliderlats(lats)

lons = divfds.m_lon.values
lon_dr = processGliderlons(lons)
lats = divfds.m_lat.values
lat_dr = processGliderlats(lats)

fig, ax = plt.subplots(figsize=(8, 8))
ax.scatter(lat_dr, lon_dr, label='DR position', s = 6)
ax.scatter(lat_gps, lon_gps, label='GPS position', s = 6)

ax.annotate(xy = (lat_gps[0], lon_gps[0]), xycoords='data',
            text= pd.to_datetime(divfds.m_gps_lon.time.values[0]).strftime('%H:%M'))
gps_ind = pd.Series(lat_gps).last_valid_index()
ax.annotate(xy = (lat_gps[gps_ind], lon_gps[gps_ind]), xycoords='data',
            text= pd.to_datetime(divfds.m_gps_lon.time.values[gps_ind]).strftime('%H:%M'))
# ax.annotate(xy = (lat_dr[0], lon_dr[0]), xycoords='data',
#             text= pd.to_datetime(divfds.m_lon.time.values[0]).strftime('%H:%M'),
#             horizontalalignment='right')
# dr_ind = pd.Series(lat_dr).last_valid_index()
# ax.annotate(xy = (lat_dr[dr_ind], lon_dr[dr_ind]), xycoords='data',
#             text= pd.to_datetime(divfds.m_lon.time.values[-1]).strftime('%H:%M'),
#             horizontalalignment='right')

ax.ticklabel_format(useOffset=False)
#ax.grid(axis = 'both')
plt.legend()
plt.suptitle(f'Tracks for dive:{divename}',fontweight="bold")
plt.title(f'{fullfn}', fontsize = 8)

thistitle = f'a_Tracks_{divename}'
plt.savefig(divegraphdir+thistitle+'.png', format='png')

### Time - Depth - Pitch
fig, ax1 = plt.subplots(figsize=(8, 4))
ax2 = ax1.twinx()
ax1.plot(divfds.time.values, divfds.m_depth.values, label='m_depth',
         lw=1, marker='o', markersize=3)
ax2.plot(divfds.time.values, np.rad2deg(divfds.m_pitch.values), label='m_pitch',
         c='r', lw=1, marker='o', markersize=3)
ax2.plot(divfds.time.values, np.rad2deg(divfds.c_pitch.values), label='c_pitch',
         c='g', lw=1, marker='o', markersize=3)

ax1.invert_yaxis()
ax1.set_ylabel('Depth (m)')
ax1.grid(axis = 'x')
ax2.grid(axis = 'both')

ax2.set_ylabel('Pitch (°)', c='r')
ax2.tick_params(axis='y', color='r')
ax2.yaxis.label.set_color('red')
ax2.spines['right'].set_edgecolor('r')
ax2.tick_params(axis='y', colors='r')
ax2.axhline(26, color='red', linestyle='dotted', linewidth=1.5)
ax2.axhline(-26, color='red', linestyle='dotted', linewidth=1.5)

plt.legend()
plt.suptitle(f'Time-Depth-Pitch: {divename}', fontweight="bold")
plt.title(f'{fullfn}', fontsize = 8)

thistitle = f'b_Time-Depth-Pitch_{divename}'
plt.savefig(divegraphdir+thistitle+'.png', format='png')

### Time - Depth - Battpos
fig, ax1 = plt.subplots(figsize=(8, 4))
ax2 = ax1.twinx()
ax1.plot(divfds.time.values, divfds.m_depth.values, label='m_depth',
         lw=1, marker='o', markersize=3)
ax2.plot(divfds.time.values, divfds.m_battpos.values, label='m_battpos',
         c='r', lw=1, marker='o', markersize=3)
ax2.plot(divfds.time.values, divfds.c_battpos.values, label='c_battpos',
         c='g', lw=1, marker='o', markersize=3)

ax1.invert_yaxis()
ax1.set_ylabel('Depth (m)')
ax1.grid(axis = 'x')
ax2.set_ylabel('Battery Position (in)', c='r')
ax2.tick_params(axis='y', color='r')
ax2.yaxis.label.set_color('red')
ax2.spines['right'].set_edgecolor('r')
ax2.tick_params(axis='y', colors='r')
ax2.grid(axis = 'both')

plt.legend()
plt.suptitle(f'Time-Depth-Battpos: {divename}', fontweight="bold")
plt.title(f'{fullfn}', fontsize = 8)

props = dict(boxstyle='round', facecolor='red', alpha=0.2)
txtstr = f'Max: {np.nanmax(divfds.m_battpos.values):.2f}" \n Min: {np.nanmin(divfds.m_battpos.values):.2f}"'
plt.text(0.85,1.05, txtstr,transform=ax1.transAxes, 
         fontsize=10,bbox = props)
thistitle = f'c_Time-Depth-Battpos_{divename}'
plt.savefig(divegraphdir+thistitle+'.png', format='png')

### Time - Depth - Oil_vol
fig, ax1 = plt.subplots(figsize=(8, 4))
ax2 = ax1.twinx()
ax1.plot(divfds.time.values, divfds.m_depth.values, label='m_depth',
         lw=1, marker='o', markersize=3)
ax1.plot(divfds.time.values, divfds.m_air_pump.values, label='m_air_pump',
         c = 'orange', lw=1, marker='o', markersize=3, alpha = 0.2)
ax2.plot(divfds.time.values, divfds.m_de_oil_vol.values, label='m_ballast_pumped',
         c='r', lw=1, marker='o', markersize=3)
ax2.plot(divfds.time.values, divfds.c_de_oil_vol.values, label='c_ballast_pumped',
         c='g', lw=1, marker='o', markersize=3)
plt.legend()
ax1.invert_yaxis()
ax1.set_ylabel('Depth (m)')
ax1.grid(axis = 'x')
ax2.set_ylabel('1000m (cc)')
ax2.grid(axis = 'both')

plt.legend()
plt.suptitle(f'Time-Depth-Oil_vol for dive: {divename}', fontweight="bold")
plt.title(f'{fullfn}', fontsize = 8)
props = dict(boxstyle='round', facecolor='red', alpha=0.2)
txtstr = f'Max: {np.nanmax(divfds.m_de_oil_vol.values):.2f}cc \n Min: {np.nanmin(divfds.m_de_oil_vol.values):.2f}cc'
plt.text(0.85,1.05, txtstr,transform=ax1.transAxes, 
         fontsize=10,bbox = props)
         
thistitle = f'd_Time-Depth-Oil_vol_{divename}'
plt.savefig(divegraphdir+thistitle+'.png', format='png')

### Time - Depth - Depthrate
fig, ax1 = plt.subplots(figsize=(8, 4))
ax2 = ax1.twinx()
ax1.plot(divfds.time.values, divfds.m_depth.values, label='m_depth',
         lw=1, marker='o', markersize=3)
ax2.plot(divfds.time.values, np.diff(divfds.m_depth.values, prepend= 0), label='m_depthrate',
         c='r', lw=1, marker='o', markersize=3)

ax1.invert_yaxis()
ax1.set_ylabel('Depth (m)')
ax1.grid(axis = 'x')
ax2.set_ylabel('depthrate (ms$^-1$)', c='r')
ax2.tick_params(axis='y', color='r')
ax2.yaxis.label.set_color('red')
ax2.spines['right'].set_edgecolor('r')
ax2.tick_params(axis='y', colors='r')
ax2.grid(axis = 'both')

plt.suptitle(f'Time-Depth-Depthrate for dive: {divename}', fontweight="bold")
plt.title(f'{fullfn}', fontsize = 8)

thistitle = f'e_Time-Depth-Depthrate_{divename}'
plt.savefig(divegraphdir+thistitle+'.png', format='png')

### Time - Depth - Altitude
fig, ax1 = plt.subplots(figsize=(8, 4))

ax1.plot(divfds.time.values, divfds.m_depth.values, label='m_depth',
         lw=1, marker='o', markersize=3)

ax1.plot(divfds.time.values, divfds.m_water_depth.values, label='m_water_depth',
         ls='-',lw=2, marker='o', markersize=3)
ax1.plot(divfds.time.values, divfds.m_altitude.values, label='m_altitude',
         lw=2, marker='o', markersize=3)
ax1.invert_yaxis()
ax1.set_ylabel('Depth (m)')
ax1.grid()
plt.legend()
plt.suptitle(f'Time-Depth-Altitude for dive: {divename}', fontweight="bold")
plt.title(f'{fullfn}', fontsize = 8)

thistitle = f'f_Time-Depth-Altitude_{divename}'
plt.savefig(divegraphdir+thistitle+'.png', format='png')

### Time - Heading
fig, ax1 = plt.subplots(figsize=(8, 4))

ax1.plot(divfds.time.values, np.rad2deg(divfds.m_heading.values), label='m_heading',
         lw=1, marker='o', markersize=3)
ax1.plot(divfds.time.values, np.rad2deg(divfds.c_heading.values), label='c_heading',
         lw=1, marker='o', markersize=3)

ax1.set_ylabel('Degrees')
ax1.set_ylim(0,360)
ax1.grid()

degrees = [0, 45, 90, 135, 180, 225, 270, 315]
compass_labels = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']
plt.yticks(ticks=degrees, labels=compass_labels)
plt.legend()
plt.suptitle(f'Time-Heading for dive: {divename}', fontweight="bold")
plt.title(f'{fullfn}', fontsize = 8)

thistitle = f'g_Time-Heading_{divename}'
plt.savefig(divegraphdir+thistitle+'.png', format='png')

### Time - Roll
fig, ax1 = plt.subplots(figsize=(8, 4))

ax1.plot(divfds.time.values, np.rad2deg(divfds.m_roll.values), label='m_roll',
         lw=1, marker='o', markersize=3)

ax1.set_ylabel('Degrees')
ax1.set_ylim(-15,15)
ax1.grid()

plt.axhline(0, color='red', linestyle='dotted', linewidth=1.5)
plt.legend()
plt.suptitle(f'Time-Roll for dive: {divename}', fontweight="bold")
plt.title(f'{fullfn}', fontsize = 8)

thistitle = f'h_Time-Roll_{divename}'
plt.savefig(divegraphdir+thistitle+'.png', format='png')

### Time - Water_speed
# this one needs work
fig, ax1 = plt.subplots(figsize=(8, 4))

ax1.plot(divfds.time.values, divfds.m_water_vy.values, label='m_water_vy',
         lw=1, marker='o', markersize=3)
ax1.plot(divfds.time.values, divfds.m_water_vy.values, label='m_water_vy',
         lw=1, marker='o', markersize=3)

ax1.invert_yaxis()
ax1.set_ylabel('velocity (ms$^-1$)')
ax1.grid()

plt.legend()
plt.suptitle(f'Time-Water_speed for dive: {divename}', fontweight="bold")
plt.title(f'{fullfn}', fontsize = 8)

thistitle = f'i_Time-Water_speed_{divename}'
plt.savefig(divegraphdir+thistitle+'.png', format='png')

### Time - Depth - Voltage
fig, ax1 = plt.subplots(figsize=(8, 4))
ax2 = ax1.twinx()
ax1.plot(divfds.time.values, divfds.m_depth.values, label='m_depth',
         lw=1, marker='o', markersize=3)
ax2.plot(divfds.time.values, divfds.m_battery.values, label='m_voltage',
         c='r', lw=1, marker='o', markersize=3)

ax1.invert_yaxis()
ax1.set_ylabel('Depth (m)')
ax1.grid(axis = 'x')
ax2.grid(axis = 'y')
ax2.set_ylabel('Voltage (V)', c='r')
ax2.tick_params(axis='y', color='r')
ax2.yaxis.label.set_color('red')
ax2.spines['right'].set_edgecolor('r')
ax2.tick_params(axis='y', colors='r')

plt.suptitle(f'Time-Depth-Voltage for dive: {divename}', fontweight="bold")
plt.title(f'{fullfn}', fontsize = 8)

thistitle = f'j_Time-Depth-Voltage_{divename}'
plt.savefig(divegraphdir+thistitle+'.png', format='png')
### Time - Depth - Current
fig, ax1 = plt.subplots(figsize=(8, 4))
ax2 = ax1.twinx()
ax1.plot(divfds.time.values, divfds.m_depth.values, label='m_depth',
         lw=1, marker='o', markersize=3)
ax2.plot(divfds.time.values, divfds.m_coulomb_amphr.values, label='m_coulomb_amphr',
         c='r', lw=1, marker='o', markersize=3)

ax1.invert_yaxis()
ax1.set_ylabel('Depth (m)')
ax1.grid(axis = 'x')
ax2.grid(axis = 'y')
ax2.set_ylabel('Current (amps)', c='r')
ax2.yaxis.label.set_color('red')
ax2.spines['right'].set_edgecolor('r')
ax2.tick_params(axis='y', colors='r')

plt.suptitle(f'Time-Depth-Current for dive: {divename}', fontweight="bold")
plt.title(f'{fullfn}', fontsize = 8)

thistitle = f'k_Time-Depth-Current_{divename}'
plt.savefig(divegraphdir+thistitle+'.png', format='png')
### Time - Leak Detects
fig, ax1 = plt.subplots(figsize=(8, 4))
# ax2 = ax1.twinx()
ax1.plot(divfds.time.values, divfds.m_leakdetect_voltage_forward.values, 
         label='m_leakdetect_voltage_forward',
         lw=1, marker='o', markersize=3)
ax1.plot(divfds.time.values, divfds.m_leakdetect_voltage.values, 
         label='m_leakdetect_voltage',
         c='r', lw=1, marker='o', markersize=3)

ax1.invert_yaxis()
ax1.set_ylabel('Voltage (V)')
ax1.grid()
# ax2.set_ylabel('Current (amps)', c='r')
# ax2.tick_params(axis='y', color='r')
# ax2.yaxis.label.set_color('red')
# ax2.spines['right'].set_edgecolor('r')  
# ax2.grid(axis = 'y')

plt.legend()
plt.suptitle(f'Time-Leak Detects for dive: {divename}', fontweight="bold")
plt.title(f'{fullfn}', fontsize = 8)

thistitle = f'l_Time-Leak_Detects_{divename}'
plt.savefig(divegraphdir+thistitle+'.png', format='png')
### Time - digifin_leakdetect
fig, ax1 = plt.subplots(figsize=(8, 4))

ax1.plot(divfds.time.values, divfds.m_digifin_leakdetect_reading.values, 
         label='m_digifin_leakdetect_reading',
         lw=1, marker='o', markersize=3)

ax1.set_ylabel('nodim')
ax1.set_ylim(1020, 1025)
ax1.grid()

plt.legend()
plt.suptitle(f'Time-digifin_leakdetect for dive: {divename}', fontweight="bold")
plt.title(f'{fullfn}', fontsize = 8)

thistitle = f'm_Time-digifin_leakdetect_{divename}'
plt.savefig(divegraphdir+thistitle+'.png', format='png')

### Time - COND - TEMP _Press
if not skipscience:
    y_depth = bsds.depth_bins.values
    x_temp = bsds.sci_rbrctd_temperature_00.values
    x_sal = bsds.sci_rbrctd_salinity_00
    x_oxy = bsds.sci_oxy4_oxygen


    fig = plt.figure(figsize = (8,10)) #, {siteproj}
    ax_full = fig.subplots()
                     
    ax_full.invert_yaxis()
    #ax_full.set_ylim(rndup(df.depSM.max(),10),0)

    full = ax_full.plot(x_temp, y_depth, color='black', label='Temp')
    ax_full.spines['bottom'].set_color('black')
    ax_full.tick_params(axis='x', colors='black')
    ax_full.xaxis.label.set_color('black')
    ax_full.set_xlim(autorange(x_temp, 0.1))#-2, -1)
    ax_full.set_xlabel('Temperature [°C]', labelpad = 0.1)
    ax_full.set_ylabel('Depth [m]', labelpad = 0.1)
        
    ax_sal = ax_full.twiny()
    sal = ax_sal.plot(x_sal, y_depth, color='firebrick', label='Sal')
    ax_sal.xaxis.set_ticks_position('bottom')
    ax_sal.xaxis.set_label_position('bottom')
    ax_sal.spines['bottom'].set_position(('axes', -0.055))
    ax_sal.spines['bottom'].set_color('firebrick')
    ax_sal.tick_params(axis='x', colors='firebrick')
    ax_sal.xaxis.label.set_color('firebrick')
    ax_sal.set_xlim(autorange(x_sal,0.02))
    ax_sal.set_xlabel('Salinity [PSU]', labelpad = 0.1)

    ax_oxy = ax_full.twiny()
    oxy = ax_oxy.plot(x_oxy, y_depth, color='blue', label='Oxy')
    ax_oxy.xaxis.set_ticks_position('bottom')
    ax_oxy.xaxis.set_label_position('bottom')
    ax_oxy.spines['bottom'].set_position(('axes', -0.11))
    ax_oxy.spines['bottom'].set_color('blue')
    ax_oxy.tick_params(axis='x', colors='blue')
    ax_oxy.xaxis.label.set_color('blue')
    ax_oxy.set_xlim(autorange(x_oxy, 5))
    ax_oxy.set_xlabel('Oxygen [µmol/L]', labelpad = 0.1)
    
    fig.suptitle(f'TEMP - SAL - OXY for dive: {divename}',fontweight="bold")
    plt.title(f'{fullfn}', fontsize = 8)

    lines = full + sal + oxy# + deoxy
    labels = [l.get_label() for l in lines]
    #ax_full.legend(lines, labels, loc = 'lower right')
    ax_full.grid(linestyle = '--')
    plt.tight_layout()
    
    thistitle = f'n_SCIENCE_{divename}'
    plt.savefig(divegraphdir+thistitle+'.png', format='png')
    
## Optical science
    y_depth = bsds.depth_bins.values
    x_fluoro = bsds.sci_flbbcd_bb_units.values#.interpolate(method = 'polynomial', order = 1)
    x_turb = bsds.sci_flbbcd_chlor_units.values
    x_PAR = bsds.sci_satpar_par
    
    #fig = plt.figure(figsize = (10,10))
    fig = plt.figure(figsize = (8,10))
    ax_fluoro = fig.subplots()
    
    ax_fluoro.invert_yaxis()
    ax_fluoro.set_ylim(rndup(y_depth.max(),10),0)
    
    fluoro = ax_fluoro.plot(x_fluoro, y_depth, color='green', label='Fluorometry')
    # ax_fluoro.xaxis.set_ticks_position('bottom')
    # ax_fluoro.xaxis.set_label_position('bottom')
    #ax_fluoro.spines['bottom'].set_position(('axes', 0.1))
    ax_fluoro.spines['bottom'].set_color('green')
    ax_fluoro.tick_params(axis='x', colors='green')
    ax_fluoro.xaxis.label.set_color('green')
    ax_fluoro.set_xlim(autominmax(x_fluoro))
    ax_fluoro.set_xlabel('Chlorophyll [ug/l]?', labelpad = 0.1)
    
    ax_fluoro.set_ylabel('Depth [m]', labelpad = 0.1)
    
    ax_turb = ax_fluoro.twiny()
    turb = ax_turb.plot(x_turb, y_depth, color='rosybrown', label='Backscatter')
    ax_turb.xaxis.set_ticks_position('bottom')
    ax_turb.xaxis.set_label_position('bottom')
    ax_turb.spines['bottom'].set_position(('axes', -0.055))
    ax_turb.spines['bottom'].set_color('rosybrown')
    ax_turb.tick_params(axis='x', colors='rosybrown')
    ax_turb.xaxis.label.set_color('rosybrown')
    ax_turb.set_xlim(autominmax(x_turb))
    ax_turb.set_xlabel('Backsactter (nodim)', labelpad = 0.1)
    
    ax_PAR = ax_fluoro.twiny()
    PAR = ax_PAR.plot(x_PAR, y_depth, color='magenta', label='PAR')
    ax_PAR.xaxis.set_ticks_position('bottom')
    ax_PAR.xaxis.set_label_position('bottom')
    ax_PAR.spines['bottom'].set_position(('axes', -0.11))
    ax_PAR.spines['bottom'].set_color('magenta')
    ax_PAR.tick_params(axis='x', colors='magenta')
    ax_PAR.xaxis.label.set_color('magenta')
    ax_PAR.set_xlabel(r'PAR [µEm$^{-2}$s$^{-1}$]', labelpad = 0.1)
    
    fig.suptitle(f'CHLOR - BACKSCATTER - PAR for dive: {divename}',fontweight="bold")
    plt.title(f'{fullfn}', fontsize = 8)
    
    lines = fluoro + turb + PAR# + deoxy
    #labels = [l.get_label() for l in lines]
    #ax_full.legend(lines, labels, loc = 'lower right')
    ax_fluoro.grid(linestyle = '--')
    
    #fig.title(f'Ship Station {stnnum}')
    plt.tight_layout()
    
    thistitle = f'o_SCIENCE_{divename}'
    plt.savefig(divegraphdir+thistitle+'.png', format='png')

### This makes the HTML 
import os

def generate_html_for_graphs(folder_path, output_html, 
                             page_title=f'{divename} - {mdep}', table_data=None):
    """
    Generate an HTML file to display all images in a folder with a table at the top.
    
    Args:
        folder_path (str): Path to the folder containing graph images.
        output_html (str): Path to save the generated HTML file.
        page_title (str): Title of the HTML page.
        table_data (dict): Dictionary containing table data (keys and single values).
    """
    # Start HTML content with dynamic title
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>{page_title}</title>
        <style>
            body {{
                font-family: Arial, sans-serif;
                margin: 20px;
                text-align: center;
            }}
            table {{
                width: 50%;
                margin: 0 auto 20px;
                border-collapse: collapse;
                text-align: center;
            }}
            th, td {{
                border: 1px solid #ddd;
                padding: 8px;
                width: 50%; /* Equal column width */
            }}
            th {{
                background-color: #f4f4f4;
                font-weight: bold;
            }}
            img {{
                max-width: 90%;
                height: auto;
                margin-bottom: 20px;
            }}
        </style>
    </head>
    <body>
        <h1>{page_title}</h1>
    """

    # Add table if dictionary is provided
    if table_data:
        html_content += """
        <table>
        <tr><th>Metric</th><th>Value</th></tr>
        """
        for key, value in table_data.items():
            html_content += f"<tr><td>{key}</td><td>{value}</td></tr>"

        html_content += """
        </table>
        """

    # Add each image in the folder to the HTML
    for filename in sorted(os.listdir(folder_path)):
        if filename.endswith(('.png', '.jpg', '.jpeg', '.gif')):  # Include image formats
            image_path = os.path.join(divename, filename)
            html_content += f"""
            <div>
                <!--<h2>{filename}</h2>-->
                <img src="{image_path}" alt="{filename}">
            </div>
            """

    # End HTML content
    html_content += """
    </body>
    </html>
    """

    # Write the HTML content to a file
    with open(output_html, 'w') as f:
        f.write(html_content)

    print(f"HTML file created: {output_html}")

# Example usage
graphs_folder = divegraphdir  # Folder containing graph images
output_file = graphs+divename+"_graphs.html"  # Output HTML file
title = "Graphs and Metrics"  # Custom
generate_html_for_graphs(graphs_folder, output_file, table_data=divemet)

