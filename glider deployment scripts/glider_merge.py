# -*- coding: utf-8 -*-
"""
Created on Sat Aug 23 18:00:07 2025

@author: OllieTwigge
"""

import xarray as xr
import glob

def safe_flight_merge_rawnc(indir, outdir, deploymentyaml,
                     scisuffix="ebd", glidersuffix="dbd"):
    # Collect input files
    fin = glob.glob(f"{indir}/*.{glidersuffix}.nc")
    print(f"Found {len(fin)} {glidersuffix}.nc files")

    good_files = []
    bad_files = []

    for f in fin:
        try:
            ds = xr.open_dataset(f)
            if all(size > 0 for size in ds.dims.values()):
                good_files.append(f)
            else:
                bad_files.append(f)
        except Exception as e:
            bad_files.append(f)

    print(f"Using {len(good_files)} good files, skipping {len(bad_files)} bad files")

    if len(good_files) == 0:
        raise RuntimeError("No valid files to merge — check raw2nc conversion")

    with xr.open_mfdataset(good_files, decode_times=False, lock=False) as ds:
        outname = f"{outdir}/merged_{glidersuffix}.nc"
        ds.to_netcdf(outname)
        print("Merged file written to:", outname)
        
def safe_science_merge_rawnc(indir, outdir, deploymentyaml,
                     scisuffix="ebd", glidersuffix="dbd"):   
    # Collect input files
    fin = glob.glob(f"{indir}/*.{scisuffix}.nc")
    print(f"Found {len(fin)} {scisuffix}.nc files")

    good_files = []
    bad_files = []

    for f in fin:
        try:
            ds = xr.open_dataset(f)
            if all(size > 0 for size in ds.dims.values()):
                good_files.append(f)
            else:
                bad_files.append(f)
        except Exception as e:
            bad_files.append(f)

    print(f"Using {len(good_files)} good files, skipping {len(bad_files)} bad files")

    if len(good_files) == 0:
        raise RuntimeError("No valid files to merge — check raw2nc conversion")

    with xr.open_mfdataset(good_files, decode_times=False, lock=False) as ds:
        outname = f"{outdir}/merged_{scisuffix}.nc"
        ds.to_netcdf(outname)
        print("Merged file written to:", outname)

    return bad_files
#%%
# instead of slocum.merge_rawnc(...)
bad_flight_files = safe_flight_merge_rawnc(rawdir, rawdir, deploymentyaml,
                             scisuffix=scisuffix, glidersuffix=glidersuffix)
#%%
bad_science_files = safe_science_merge_rawnc(rawdir, rawdir, deploymentyaml,
                             scisuffix=scisuffix, glidersuffix=glidersuffix)
