# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 12:04:02 2025

@author: OllieTwigge
"""
# this asks for the directory with a binary and cache folders, sensorlist.txt and deploymentyaml"
# then it asks which dbd file you want to process and spits out all the graphs from that yo
#%%
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

#%%
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
        
#%%
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
    fig, ax1 = plt.subplots(figsize=(8, 4))
    
    ax1.plot(divsds.time.values, divsds.sci_rbrctd_conductivity_00.values, 
             label='sci_rbrctd_conductivity_00',
             lw=1, marker='o', markersize=3)
    ax1.plot(divsds.time.values, divsds.sci_rbrctd_temperature_00.values, 
             label='sci_rbrctd_temperature_00',
             lw=1, marker='o', markersize=3)
    ax1.plot(divsds.time.values, divsds.sci_rbrctd_pressure_00.values, 
             label='sci_rbrctd_pressure_00',
             lw=1, marker='o', markersize=3)
    ax1.invert_yaxis()
    ax1.set_ylabel('Depth (m)')
    ax1.grid(axis = 'both')
    
    plt.legend()
    plt.suptitle(f'Time - COND - TEMP - PRES for dive: {divename}', fontweight="bold")
    plt.title(f'{fullfn}', fontsize = 8)
    
    thistitle = f'n_SCIENCE_{divename}'
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

