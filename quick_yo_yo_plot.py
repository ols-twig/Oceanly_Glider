# -*- coding: utf-8 -*-
"""
Created on Fri Feb 13 15:28:14 2026

@author: OllieTwigge
"""

import xarray as xr

gds = xr.open_dataset(r'C:\temp\samoa_glider\rawnc_flat\selkie1104rawdbd.nc')
#%%
import matplotlib.pyplot as plt
import pandas as pd

# --- choose variable ---
p = gds['m_pressure'].values        # dbar (preferred)
# if not available:
# p = gds['m_depth'].values

# time is usually epoch seconds
t = gds['time'].values

# convert to datetime
t = pd.to_datetime(t, unit='s')

# remove NaNs (important)
mask = ~pd.isna(p) & ~pd.isna(t)

plt.figure(figsize=(10,4))
plt.plot(t[mask], p[mask], linewidth=0.5)
plt.gca().invert_yaxis()
plt.xlabel('Time')
plt.ylabel('Pressure (dbar)')
plt.title('Glider Yo-Yo')
plt.tight_layout()
plt.show()
