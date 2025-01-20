# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 12:04:02 2025

@author: OllieTwigge
"""
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
mettab = pd.DataFrame()
import glob
for i in glob.glob('./rawnc/*.dbd.nc'):
    divfd = maindir+i
  #need to indent from here  
    #divfd = select_file()
    divename = os.path.splitext(os.path.splitext(os.path.basename(divfd))[0])[0]
    
    divsd = os.path.splitext(os.path.splitext(os.path.basename(divfd))[0])[0]
    divsd = maindir+rawdir+divsd+'.ebd.nc'
    
    divfds = xr.open_dataset(divfd)
    # divsds = xr.open_dataset(divsd)

    divfds = divfds.set_coords('time')
    divfds = divfds.swap_dims({'_ind': 'time'})
    divfds['time'] = pd.to_datetime(divfds.time.values, unit='s')
    
    # divsds = divsds.set_coords('time')
    # divsds = divsds.swap_dims({'_ind': 'time'})
    # divsds['time'] = pd.to_datetime(divsds.time.values, unit='s')
    
    fullfn = divfds.attrs['full_filename']
    print(f'Dive chosen: {divename}')
    divegraphdir = graphs+divename+'/'
    if not os.path.exists(graphs+divename):
            os.makedirs(graphs+divename)
            print(f'CREATED dive graph folder {divegraphdir}')
    else:
            print('dive graph folder exists')
            

    #ADD META TABLE TO PRINT
    mdep = f'{np.nanmax(divfds.m_depth.values):.2f}m'
    divemet= {
        'divename': f'{divename}',
        'First Timestamp': str(pd.to_datetime(divfds.time.values[0])).split('.')[0], 
        'Last Timestamp': str(pd.to_datetime(divfds.time.values[-1])).split('.')[0],
        'Total segment time': str(timedelta(seconds = (divfds.time.values[-1]-divfds.time.values[0])/np.timedelta64(1, 's'))).split('.')[0],
        'Total Data points': len(divfds.time),
        'Max Depth': mdep,
        'Number of yos in segment':((np.nanmax(divfds.m_num_half_yos_in_segment.values)+1)/2).astype(int), 
        'Average yo time':	f'{(np.nanmean(divfds.m_avg_yo_time.values)/60):.2f} mins',
        'Average speed': f'{(np.nanmean(divfds.m_avg_speed)):.3f} m/s',
        'Average dive rate': f'{(np.nanmean(divfds.m_avg_dive_rate.values)):.3f} m/s',
        'Average climb rate': f'{(np.nanmean(divfds.m_avg_climb_rate.values)):.3f} m/s',
        'Average downward inflection time':	f'{(np.nanmean(divfds.m_avg_downward_inflection_time)).astype(int)}s',
        #'Science sensors': list(divsds.keys())
        }
    
    
    thismet = pd.DataFrame.from_dict(divemet, orient = 'index').T
    mettab = pd.concat([mettab, thismet])
    #%%
mettab.to_excel(graphs+'deployment_meta.xlsx')