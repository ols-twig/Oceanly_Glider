# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 15:30:12 2026

@author: OllieTwigge
"""
import xarray as xr
import numpy as np

gds1 = xr.open_dataset(r"G:\.shortcut-targets-by-id\1Z7PX__ZLrJ6oL2wHtDtiATh5eAh4Cgg3\Oceanly Science\Expeditions\Samoa 2025\Data\Glider\proc_merged\1\L0-timeseries\2025_08_24_0237_selkie_Samoa_2025.nc")
gds2 = xr.open_dataset(r"G:\.shortcut-targets-by-id\1Z7PX__ZLrJ6oL2wHtDtiATh5eAh4Cgg3\Oceanly Science\Expeditions\Samoa 2025\Data\Glider\proc_merged\2\L0-timeseries\2025_09_20_2231_selkie_Samoa_2025.nc")

print(gds1.time.min().values, gds1.time.max().values)
print(gds2.time.min().values, gds2.time.max().values)

overlap = gds1.time.max() >= gds2.time.min()
print("Overlap:", overlap)

#%%
gds_combined = (
    xr.concat([gds1, gds2], dim="time")
    .sortby("time")
    .drop_duplicates(dim="time")
)
#%%

gds_combined.attrs['time_coverage_start'] = str(gds_combined.time.min().values)
gds_combined.attrs['time_coverage_end'] = str(gds_combined.time.max().values)
#%%

outpath = r"G:\.shortcut-targets-by-id\1Z7PX__ZLrJ6oL2wHtDtiATh5eAh4Cgg3\Oceanly Science\Expeditions\Samoa 2025\Data\Glider\L0-timeseries\selkie_20250824T0237_Samoa_2025.nc"
gds_combined.to_netcdf(outpath)

#%%

# Concatenate
gds_combined = xr.concat([gds1, gds2], dim="time")

# Sort by time
gds_combined = gds_combined.sortby("time")

# Drop duplicate timestamps (keep first occurrence)
_, index = np.unique(gds_combined['time'], return_index=True)
gds_combined = gds_combined.isel(time=index)