# -*- coding: utf-8 -*-
"""
Created on Thu Feb 12 12:38:58 2026

@author: OllieTwigge
"""

import dolfyn
import os
import glob

import sys
sys.path.append(r'C:\temp\\samoa_glider_dvl\\scripts\\')

import PathfinderTimeSeries

os.chdir(r'C:\temp\\samoa_glider_dvl\\')
filepath = r'./data/'
files = glob.glob(os.path.join(filepath, '*.pd0'))

#%%
from pathlib import Path
import struct
import numpy as np
from datetime import datetime
import xarray as xr

SYNC = b"\x7f\x7f"

def year_from_byte(y: int) -> int:
    return 2000 + y if y < 80 else 1900 + y

def parse_pd0_directory(ens: bytes):
    # PD0 header: 0-1 sync, 2-3 nbytes (little endian), 4 spare, 5 ndat, 6.. offsets (uint16)
    if ens[:2] != SYNC:
        raise ValueError("Not a PD0 ensemble")
    nbytes = struct.unpack_from("<H", ens, 2)[0]
    ndat = ens[5]
    offsets = [struct.unpack_from("<H", ens, 6 + 2*i)[0] for i in range(ndat)]
    return nbytes, ndat, offsets

def parse_time_from_varleader(var: bytes) -> datetime:
    # var starts with ID=0x0080 (2 bytes), then ensemble number (2),
    # then RTC: year, month, day, hour, minute, second, hundredths
    y = year_from_byte(var[4])
    mo, d = var[5], var[6]
    hh, mm, ss = var[7], var[8], var[9]
    hund = var[10]
    return datetime(y, mo, d, hh, mm, ss, hund * 10_000)  # hundredths -> 10ms

def iter_ensembles(blob: bytes):
    i = 0
    while True:
        j = blob.find(SYNC, i)
        if j == -1 or j + 4 > len(blob):
            break
        nbytes = struct.unpack_from("<H", blob, j + 2)[0]
        if j + nbytes > len(blob):
            i = j + 2
            continue
        yield blob[j:j+nbytes]
        i = j + nbytes

def extract_corr_intensity(pd0_path: Path) -> xr.Dataset:
    blob = pd0_path.read_bytes()

    times = []
    corr_list = []
    inten_list = []

    for ens in iter_ensembles(blob):
        _, _, offsets = parse_pd0_directory(ens)

        # offsets correspond to: fixed, var, vel, corr, inten, pg, bt (commonly)
        var = ens[offsets[1]:offsets[2]]
        t = parse_time_from_varleader(var)

        corr_blk = ens[offsets[3]:offsets[4]]
        int_blk  = ens[offsets[4]:offsets[5]]

        corr = np.frombuffer(corr_blk[2:], dtype=np.uint8)   # strip ID
        inten = np.frombuffer(int_blk[2:], dtype=np.uint8)

        beams = 4  # Pathfinder DVL typical
        cells = len(corr) // beams
        corr = corr[:cells*beams].reshape(cells, beams)
        inten = inten[:cells*beams].reshape(cells, beams)

        times.append(np.datetime64(t))
        corr_list.append(corr)
        inten_list.append(inten)

    corr_arr = np.stack(corr_list, axis=0)      # (time, cell, beam)
    inten_arr = np.stack(inten_list, axis=0)

    ds = xr.Dataset(
        data_vars=dict(
            correlation=(("time","cell","beam"), corr_arr),
            echo_intensity=(("time","cell","beam"), inten_arr),
        ),
        coords=dict(
            time=np.array(times),
            cell=np.arange(corr_arr.shape[1]),
            beam=np.arange(corr_arr.shape[2]),
        ),
        attrs=dict(source_file=str(pd0_path.name))
    )
    return ds

def batch_dir(in_dir: str, out_dir: str):
    in_dir = Path(in_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for f in sorted(in_dir.glob("*.pd0")):
        ds = extract_corr_intensity(f)
        ds.to_netcdf(out_dir / (f.stem + "_corr_int.nc"))
        print("Wrote", f.name)

# Example:
# batch_dir(r"C:\data\pd0", r"C:\data\out")
#%% This will convert all of the pd0 to 
# testing
# pd0_file = open(files[28], 'rb').read()
#%% THIS ACTUALLY RUNS ALL OF THE PD0s
batch_dir(r'C:\temp\samoa_glider_dvl\data', './converted/')

#%% USe this to combine them all into one timeseries after 

from pathlib import Path
import xarray as xr
import numpy as np

def combine_corr_int_nc(in_dir: str, out_nc: str, pattern: str = "*_corr_int.nc"):
    in_dir = Path(in_dir)
    files = sorted(in_dir.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files matched {pattern} in {in_dir}")

    ds = xr.open_mfdataset(
        files,
        combine="nested",
        concat_dim="time",
        parallel=False,
        engine="netcdf4",
        coords="minimal",
        data_vars="minimal",
        compat="override",
        join="exact",
    )

    ds = ds.sortby("time")

    # Drop duplicate timestamps (keep first)
    t = ds["time"].values
    keep = np.concatenate(([True], t[1:] != t[:-1]))
    ds = ds.isel(time=keep)

    ds.attrs["source_files"] = ", ".join([f.name for f in files])

    out_nc = Path(out_nc)
    out_nc.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(out_nc)
    ds.close()

    print(f"Wrote combined file: {out_nc}  (from {len(files)} inputs)")

# Example:
# combine_corr_int_nc(r"C:\temp\samoa_glider_dvl\converted", r"C:\temp\samoa_glider_dvl\ALL_corr_int.nc")
#%%

# Example:
combine_corr_int_nc(r"C:\temp\samoa_glider_dvl\converted", r"C:\data\out\ALL_corr_int.nc")


#%% Some plotting and testing

ds = xr.open_dataset(r"C:\temp\samoa_glider_dvl\test\YI092038_corr_int.nc")
#%%
import xarray as xr
import matplotlib.pyplot as plt

ds = xr.open_dataset("test\YI092038_corr_int.nc")

# Select beam 0 (change to 1,2,3 if needed)
da = ds.echo_intensity.sel(beam=0)

plt.figure(figsize=(10,4))
da.T.plot(
    x="time",
    y="cell",
    cmap="viridis",
    robust=True
)
plt.gca().invert_yaxis()   # deeper cells downwards
plt.title("Echo Intensity – Beam 0")
plt.tight_layout()
plt.show()
