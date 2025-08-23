# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 17:31:33 2025

@author: OllieTwigge
"""
from glob import glob
import os

mainpath = r"G:/.shortcut-targets-by-id/1Z7PX__ZLrJ6oL2wHtDtiATh5eAh4Cgg3/Oceanly Science/Expeditions/Solomon Islands 2025/Data/GLIDER_DATA/rawnc"
testpath = mainpath+r"\*.dbd.nc"

fullpath= list(glob(testpath, recursive=False))

#%%
for divfd in fullpath:
    
    divename = os.path.splitext(os.path.splitext(os.path.basename(divfd))[0])[0]
    divsd = os.path.splitext(os.path.splitext(os.path.basename(divfd))[0])[0]
    divsd = maindir+rawdir+divsd+'.ebd.nc'
    
    print(f'\nDive {divename}')
    divfds = xr.open_dataset(divfd)
    divsds = xr.open_dataset(rf"{divsd}")
    
    depth_var = divsds['sci_rbrctd_depth_00']
    depth_min = depth_var.min().item()
    depth_max = depth_var.max().item()
    try:
        print(f'{divfds["m_depth"].min().item()} to {divfds["m_depth"].min().item()}')
        print(f"{depth_min:.2f} to {depth_max:.2f}")
    except:
        print('NOPE')