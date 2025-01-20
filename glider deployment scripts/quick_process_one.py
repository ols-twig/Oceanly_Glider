# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 15:45:31 2025

@author: OllieTwigge
"""
import pyglider.slocum as slocum
import tkinter as tk
from tkinter import filedialog
import os

def select_file(filetypes=None, title = None):
    # Create a visible root window (used as the parent)
    root = tk.Tk()
    root.attributes('-topmost', True)  # Ensure the root window is on top
    root.withdraw()  # Hide the root window
    
    if filetypes is not None:
        filetypes = filetypes
    else:
        filetypes=[('Flight data', '*.dbd.nc'),
                   ('All NetCDFs', '*.nc'),
                   ('All Files', '*.*')]
        
    # Show the file selection dialog with the root as the parent
    file_path = filedialog.askopenfilename(title=title, #askopenfilename 
                            filetypes = filetypes,
                            parent=root)

    root.destroy()  # Destroy the root window
    return file_path
def select_directory(title = None):
    # Create a visible root window (used as the parent)
    root = tk.Tk()
    root.attributes('-topmost', True)  # Ensure the root window is on top
    root.update()  # Update the root window to apply the attributes
    
    if title is not None:
        title = title
    else:
        title="Select directory with binary, cache, sensorlist and deploymentyaml"
    # Show the directory selection dialog with the root as the parent
    directory = filedialog.askdirectory(title= title,
                                        parent=root)

    root.destroy()  # Destroy the root window
    return directory

#%%
maindir = select_directory()
os.chdir(maindir)

binarydir = './oldraw_dbd_ebd/'#select_directory(title = "Select dbd.file")
rawdir = './oldrawnc/'
cacdir = './cac/'
sensorlist = './unit_1104_sensors.txt'
deploymentyaml = './deploymentRealtime.yml'

slocum.binary_to_rawnc(
    binarydir, rawdir, cacdir, sensorlist, deploymentyaml,
    incremental=True, scisuffix='ebd', glidersuffix='dbd')