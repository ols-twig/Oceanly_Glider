# -*- coding: utf-8 -*-
"""
Created on Mon Feb 16 20:02:29 2026

@author: OllieTwigge
"""
#%%
# Using this to make a good start at QA QCing the glider data,
# maybe i make a small subset to do the testing first 
# 



#%% FIRST TO correct the lat lons 
import numpy as np
import xarray as xr

def ddmm_to_decimal(ddmm):
    """
    Convert NMEA-like ddmm.mmmm (lat) or dddmm.mmmm (lon) to decimal degrees.
    Works for negative values too.
    """
    ddmm = np.asarray(ddmm, dtype="float64")
    sign = np.sign(ddmm)
    x = np.abs(ddmm)

    deg = np.floor(x / 100.0)
    minutes = x - deg * 100.0
    dec = deg + minutes / 60.0
    return sign * dec
#%% FIRST TO correct the lat lons 
ds = xr.open_dataset(r"C:\temp\samoa_glider\QAQC\selkie_Samoa_2025.nc")

ds["latitude"]  = (("time",), ddmm_to_decimal(ds["latitude"].values))
ds["longitude"] = (("time",), ddmm_to_decimal(ds["longitude"].values))

# Optional: update attrs to reflect they're now truly decimal degrees
ds["latitude"].attrs["comment"]  = (ds["latitude"].attrs.get("comment","") + " converted from ddmm.mmmm to decimal degrees").strip()
ds["longitude"].attrs["comment"] = (ds["longitude"].attrs.get("comment","") + " converted from ddmm.mmmm to decimal degrees").strip()

#%% #use this to look at where to chooose the extents of the location 
# Will need to zoom in closer and see when the track looks right
import matplotlib.pyplot as plt
plt.plot(ds.latitude.values, ds.longitude.values)
# -172.70, -172.90 (longitude)
# -13.66, -13.45 (latitude)

# bbox (inclusive)
LON_MIN, LON_MAX = -172.90, -172.70
LAT_MIN, LAT_MAX = -13.66, -13.45

in_box = (
    (ds["longitude"] >= LON_MIN) & (ds["longitude"] <= LON_MAX) &
    (ds["latitude"]  >= LAT_MIN) & (ds["latitude"]  <= LAT_MAX)
)

# quick counts
n_total = ds.dims["time"]
n_in = int(in_box.sum().values)
print("time samples total:", n_total)
print("time samples in box:", n_in)
print("time samples NOT in box:", n_total-n_in)

#%%
# Make a small subset (e.g., first 100k in-box points) for speed
ds_in_box = ds.where(in_box, drop=True)

#%%
ds_small.to_netcdf("selkie_Samoa_2025_subset_box.nc")
print("wrote selkie_Samoa_2025_subset_box.nc")
#%% Plotting that track in cartopy

import xarray as xr
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ds = xr.open_dataset(NC)
lat = ds["latitude"].values
lon = ds["longitude"].values

# downsample (~100k points)
n = lat.size
step = max(1, n // 100_000)
lat_ds = lat[::step]
lon_ds = lon[::step]

# QC box
LON_MIN, LON_MAX = -172.90, -172.70
LAT_MIN, LAT_MAX = -13.66, -13.45

# Map extent (pad a bit around your box)
pad_lon = 1.0
pad_lat = 1.0
extent = [LON_MIN - pad_lon, LON_MAX + pad_lon,
          LAT_MIN - pad_lat, LAT_MAX + pad_lat]

fig = plt.figure(figsize=(10, 7))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent(extent, crs=ccrs.PlateCarree())

ax.add_feature(cfeature.LAND, linewidth=0)
ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax.add_feature(cfeature.BORDERS, linewidth=0.4)

gl = ax.gridlines(draw_labels=True, linewidth=0.3)
gl.top_labels = False
gl.right_labels = False

# Track
ax.plot(lon_ds, lat_ds, linewidth=1.0, transform=ccrs.PlateCarree())

# QC box
ax.plot([LON_MIN, LON_MAX, LON_MAX, LON_MIN, LON_MIN],
        [LAT_MIN, LAT_MIN, LAT_MAX, LAT_MAX, LAT_MIN],
        linewidth=2, transform=ccrs.PlateCarree())

ax.set_title("Glider track (downsampled) + QC box")
plt.tight_layout()
plt.show()

