# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 01:44:19 2025

@author: OllieTwigge
"""

import xarray as xr
import numpy as np
ds = xr.open_dataset(r"C:\temp\solo_25_glider_process\merged_dbd.nc")

np.nanmax(ds.m_depth.values)

#%%
import pandas as pd
import numpy as np

# assume df is your dataframe with m_depth as a column
depth = ds['m_depth'].values

# derivative sign (positive = going deeper, negative = going shallower)
slope = np.sign(np.diff(depth))

# find where slope changes sign (turning points)
turns = np.where(np.diff(slope) != 0)[0]

# number of dives = number of min→max transitions
# number of climbs = number of max→min transitions
dives = np.sum((slope[turns] > 0) & (slope[turns+1] < 0))   # top turns
climbs = np.sum((slope[turns] < 0) & (slope[turns+1] > 0)) # bottom turns

print("Dives:", dives)
print("Climbs:", climbs)
print("Total cycles:", min(dives, climbs))