# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 00:48:00 2025

@author: OllieTwigge
"""

import xarray as xr
import pandas as pd

ds = xr.open_dataset(r"G:\.shortcut-targets-by-id\1Z7PX__ZLrJ6oL2wHtDtiATh5eAh4Cgg3\Oceanly Science\Expeditions\Samoa 2025\Data\Glider\proc_merged\1\selkie1104rawdbd_1.nc")
surf = ds[['m_appear_to_be_at_surface','m_lat','m_lon','time']].to_pandas()

ds1 = xr.open_dataset(r"G:\.shortcut-targets-by-id\1Z7PX__ZLrJ6oL2wHtDtiATh5eAh4Cgg3\Oceanly Science\Expeditions\Samoa 2025\Data\Glider\proc_merged\2\selkie1104rawdbd_2.nc")
surf1 = ds1[['m_appear_to_be_at_surface','m_lat','m_lon','time']].to_pandas()
#%%
surf = pd.concat([surf,surf1])
#%%
surf['dt'] = pd.to_datetime(surf['time'], unit="s", utc=True)
surf = surf.drop(columns = 'time')
surf["dt"] = surf["dt"].dt.floor("s")
surf = surf.reset_index()
surf = surf.set_index('dt')
#%%

surf = surf[surf['m_appear_to_be_at_surface'] == 1]
surf = surf.dropna()
surf = surf.drop_duplicates(subset=['m_lat','m_lon'])
surf['m_lat'] = surf['m_lat'].round(3)
surf['m_lon'] = surf['m_lon'].round(3)
surf = surf.dropna()
surf = surf.drop_duplicates(subset=['m_lat','m_lon'])


#%%

def dms_to_dd(value):
    """Convert float-encoded DMS to decimal degrees."""
    neg = value < 0
    value = abs(value)

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
surf = surf[['lat_dd','lon_dd','_ind']]
surf['lat_dd'] = surf['lat_dd'].round(3)
surf['lon_dd'] = surf['lon_dd'].round(3)
surf = surf.drop_duplicates(subset=['lat_dd','lon_dd'])

#%% COmparing cause some were way off
oldsurfpath = r"G:\.shortcut-targets-by-id\1Z7PX__ZLrJ6oL2wHtDtiATh5eAh4Cgg3\Oceanly Science\Expeditions\Samoa 2025\Data\Glider\proc_merged\OLDglidersurfacings.csv"
OLDsurf = pd.read_csv(oldsurfpath)
#%%
missing = surf.loc[~surf["_ind"].isin(OLDsurf["_ind"]), "_ind"]
print(missing)
#%% removes those values
surf = surf[~surf["_ind"].isin(missing.values)]

#%%
surf.to_csv(r"G:\.shortcut-targets-by-id\1Z7PX__ZLrJ6oL2wHtDtiATh5eAh4Cgg3\Oceanly Science\Expeditions\Samoa 2025\Data\Glider\proc_merged\glidersurfacings.csv")