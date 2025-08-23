# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 09:57:22 2025

@author: OllieTwigge
"""

# -*- coding: utf-8 -*-

# import pyglider.example_data as pexamp 
# pexamp.get_example_data('./') # already having isssues with this

#%%
# -*- coding: utf-8 -*-
import logging
import glob, os
import xarray
import pandas as pd
import pyglider.ncprocess as ncprocess
import pyglider.slocum as slocum
import pyglider.utils as pgutils

#%%
os.chdir(r'C:\temp\glider_dev')
#%%
logging.basicConfig(level='INFO')

binarydir = './raw_dbd_ebd/'
rawdir = './rawnc/'
cacdir = './cac/'
sensorlist = './unit_1104_sensors.txt'
deploymentyaml = './deploymentRealtime.yml'
l1tsdir = './L0-timeseries/'
profiledir = './L0-profiles/'
griddir = './L0-gridfiles/'
scisuffix = 'ebd'#'tbd'
glidersuffix = 'dbd' #'sbd'
#%%
def checkCACfiles(binarydir, cacdir):
    missingcacs = []
    caclst = [os.path.splitext(x)[0] for x in os.listdir(cacdir)] # make a list of all the availbe cac files
    for file in os.listdir(binarydir):
        if file.endswith(('.ebd', '.dbd', '.tbd', '.sbd')):
            #print(os.path.join(binarydir, file))
            binfilepath = os.path.join(binarydir, file)
            with open(binfilepath, 'rb') as binary_file:
                for line in binary_file:
                    if line.startswith(b'sensor_list_crc'):
                        #print(f"Match found: {line}")
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

#%%
checkCACfiles(binarydir, cacdir)
 
#%%
slocum.binary_to_rawnc(
    binarydir, rawdir, cacdir, sensorlist, deploymentyaml,
    incremental=True, scisuffix=scisuffix, glidersuffix=glidersuffix)
#%%
slocum.merge_rawnc(
    rawdir, rawdir, deploymentyaml,
    scisuffix=scisuffix, glidersuffix=glidersuffix)
#%%  
from tkinter.filedialog import askopenfilename
filename = askopenfilename()
print(filename)