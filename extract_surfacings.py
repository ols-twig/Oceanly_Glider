# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 00:48:00 2025

@author: OllieTwigge
"""

import xarray as xr

ds = xr.open_dataset(r"C:\temp\solo_25_glider_process\merged_dbd.nc")
surf = ds[['m_appear_to_be_at_surface','m_lat','m_lon']].to_pandas()
#%%

surf = surf[surf['m_appear_to_be_at_surface'] == 1]
surf = surf.dropna()
surf = surf.drop_duplicates()
surf['m_lat'] = surf['m_lat'].round(3)
surf['m_lon'] = surf['m_lon'].round(3)
surf = surf.dropna()
surf = surf.drop_duplicates()

#%%
import pandas as pd

def dms_to_dd(value):
    """Convert float-encoded DMS to decimal degrees."""
    neg = value < 0
    value = abs(value)s

    if value < 10000:  # latitude (DDMM.MMM)
        degrees = int(value // 100)
        minutes = int(value % 100)
        seconds = (value - (degrees * 100 + minutes)) * 60
    else:              # longitude (DDDMM.MMM)
        degrees = int(value // 100)
        minutes = int(value % 100)
        seconds = (value - (degrees * 100 + minutes)) * 60

    dd = degrees + minutes / 60 + seconds / 3600
    if neg:
        dd *= -1
    return dd

# Example usage on your DataFrame
surf["lat_dd"] = surf["m_lat"].apply(dms_to_dd)
surf["lon_dd"] = surf["m_lon"].apply(dms_to_dd)

print(surf[["m_lat", "lat_dd", "m_lon", "lon_dd"]].head())
#%%
surf = surf[['lat_dd','lon_dd']]
surf['lat_dd'] = surf['lat_dd'].round(3)
surf['lon_dd'] = surf['lon_dd'].round(3)
surf = surf.drop_duplicates()

#%%
surf.to_csv(r'G:\.shortcut-targets-by-id\1Z7PX__ZLrJ6oL2wHtDtiATh5eAh4Cgg3\Oceanly Science\Expeditions\Solomon Islands 2025\Data\Glider\glidersurfacings.csv')