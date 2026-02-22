# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 15:44:58 2026

@author: OllieTwigge
"""
import xarray as xr

gds = xr.open_dataset(r"G:\.shortcut-targets-by-id\1Z7PX__ZLrJ6oL2wHtDtiATh5eAh4Cgg3\Oceanly Science\Expeditions\Samoa 2025\Data\Glider\L0-timeseries\selkie_20250824T0237_Samoa_2025.nc")
dvlds = xr.open_dataset(r"G:\.shortcut-targets-by-id\1Z7PX__ZLrJ6oL2wHtDtiATh5eAh4Cgg3\Oceanly Science\Expeditions\Samoa 2025\Data\Glider\L0-timeseries\ALL_dvl_ei_corr.nc")

#%%
"""
Add absolute depth(time, cell) to a DVL xarray Dataset (dvlds) using glider depth(time) (gds).

Assumes you already have:
    gds   : xarray.Dataset with coords time and data_var "depth" (glider depth, positive down)
    dvlds : xarray.Dataset with coords time, cell (and likely beam) and DVL vars (e.g., echo_intensity)

This script:
  1) Sorts by time
  2) Crops DVL to overlap with glider time (optional but recommended)
  3) Interpolates glider depth onto DVL time
  4) Converts DVL cell index -> range -> vertical offset using beam angle
  5) Adds dvlds["depth"] with dims (time, cell)

NO FILE I/O here; it just modifies dvlds in memory.
"""

import numpy as np
import xarray as xr

# ---------------------------
# USER SETTINGS (EDIT THESE)
# ---------------------------
cell_size_m = 1.0           # meters per cell (bin length)
blanking_distance_m = 0.5   # meters (no data near transducer)
beam_angle_deg = 30.0       # degrees from vertical (Teledyne often 30)
down_looking = True         # set False if upward-looking
dvl_head_offset_m = 0.0     # +ve means DVL is deeper than glider pressure sensor

# How to handle times where DVL extends beyond glider times:
#   "crop"   = drop DVL times outside glider range (recommended / strict)
#   "nearest"= fill glider depth using nearest value (pragmatic)
edge_handling = "crop"      # choose: "crop" or "nearest"


# ---------------------------
# 0) Basic checks
# ---------------------------
if "time" not in gds.coords:
    raise ValueError("gds must have a 'time' coordinate")
if "time" not in dvlds.coords:
    raise ValueError("dvlds must have a 'time' coordinate")
if "depth" not in gds.data_vars and "depth" not in gds.coords:
    raise ValueError("gds must contain a 'depth' variable (glider depth)")

# Make sure 'depth' refers to the data variable (not a coord)
glider_depth = gds["depth"]

# ---------------------------
# 1) Sort by time (interp expects monotonic increasing time)
# ---------------------------
gds = gds.sortby("time")
dvlds = dvlds.sortby("time")

# ---------------------------
# 2) Handle edge overlap (avoid NaNs at start/end)
# ---------------------------
g_time_min = gds["time"].min()
g_time_max = gds["time"].max()

if edge_handling.lower() == "crop":
    dvlds = dvlds.sel(time=slice(g_time_min, g_time_max))
elif edge_handling.lower() == "nearest":
    # no cropping; we will fill with nearest glider depth outside bounds
    pass
else:
    raise ValueError("edge_handling must be 'crop' or 'nearest'")

# ---------------------------
# 3) Interpolate (or reindex) glider depth onto DVL time
# ---------------------------
if edge_handling.lower() == "nearest":
    glider_depth_on_dvl_time = glider_depth.reindex(time=dvlds["time"], method="nearest")
else:
    glider_depth_on_dvl_time = glider_depth.interp(time=dvlds["time"], method="linear")

# Apply fixed sensor offset if needed
glider_depth_on_dvl_time = glider_depth_on_dvl_time + dvl_head_offset_m

# ---------------------------
# 4) Convert cell index -> range (bin center) -> vertical offset
# ---------------------------
# Use bin centers: R = blank + (cell + 0.5) * cell_size
cell = dvlds["cell"].astype("float64")
range_m = blanking_distance_m + (cell + 0.5) * cell_size_m

theta = np.deg2rad(beam_angle_deg)
z_offset_m = range_m * np.cos(theta)  # positive "down" component

# ---------------------------
# 5) Absolute depth(time, cell) and attach to dvlds
# ---------------------------
# Broadcasting creates (time, cell)
if down_looking:
    depth_abs = glider_depth_on_dvl_time + z_offset_m
else:
    depth_abs = glider_depth_on_dvl_time - z_offset_m

# Attach as DataArray (cleanest xarray-native way)
dvlds["depth"] = depth_abs.astype("float32")

dvlds["depth"].attrs.update({
    "long_name": "Depth of DVL measurement cell center",
    "units": "m",
    "positive": "down",
    "comment": (
        "Computed by interpolating glider depth to DVL ping times and adding "
        f"vertical cell offsets. blank={blanking_distance_m} m, cell_size={cell_size_m} m, "
        f"beam_angle={beam_angle_deg} deg, "
        f"{'down-looking' if down_looking else 'up-looking'}, "
        f"dvl_head_offset={dvl_head_offset_m} m, edge_handling={edge_handling}."
    ),
})

# ---------------------------
# 6) Quick sanity prints
# ---------------------------
print("dvlds now has variables:", list(dvlds.data_vars))
print("dvlds['depth']:", dvlds["depth"])
print("Depth min/max (ignoring NaNs):",
      float(dvlds["depth"].min(skipna=True)),
      float(dvlds["depth"].max(skipna=True)))
#%%
dvlds = dvlds.sel(time=slice(gds.time.min(), gds.time.max()))

#%%

dvlds.to_netcdf(r"G:\.shortcut-targets-by-id\1Z7PX__ZLrJ6oL2wHtDtiATh5eAh4Cgg3\Oceanly Science\Expeditions\Samoa 2025\Data\Glider\L0-timeseries\ALL_dvl_ei_corr_depth.nc")