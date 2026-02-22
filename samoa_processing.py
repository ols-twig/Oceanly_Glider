# -*- coding: utf-8 -*-
"""
Created on Mon Jan 12 15:49:18 2026

@author: OllieTwigge
"""

# -*- coding: utf-8 -*-
import logging
import os
import pyglider.ncprocess as ncprocess
import pyglider.slocum as slocum
import pyglider.utils as pgutils

logging.basicConfig(level='INFO')

binarydir = './raw_dbd_ebd/'
rawdir = './rawnc/'
cacdir = './cac/'
sensorlist = './selkie_sensors.txt'
deploymentyaml = './deploymentRealtime.yml'
l1tsdir = './L0-timeseries/'
profiledir = './L0-profiles/'
griddir = './L0-gridfiles/'
scisuffix = 'ebd'
glidersuffix = 'dbd'
#%%
os.chdir(r'C:\temp\samoa_glider')


#%%
do_direct = False
# only do this for a real run, or something like this
real = False
#%%
if real:
    os.system('rsync -av cproof@sfmc.webbresearch.com:/var/opt/sfmc-dockserver/' +
              'stations/dfo/gliders/ ~/gliderdata/slocum_dockserver/')
    os.system('rsync -av ~/gliderdata/slocum_dockserver/rosie_713/from-glider/* ' +
              binarydir)

    os.system('rm ' + rawdir + 'dfo* ' + rawdir + 'TEMP*.nc ' + l1tsdir + '* ' +
              profiledir + '* ' + griddir + '* ')
#%%
if do_direct:
    # turn *.sdb and *.tbd into timeseries netcdf files
    outname = slocum.binary_to_timeseries(
        binarydir, cacdir, l1tsdir, deploymentyaml, search='*.[s|t]bd',
        profile_filt_time=20, profile_min_time=20)
else:
    # turn *.EBD and *.DBD into *.ebd.nc and *.dbd.nc netcdf files.
    slocum.binary_to_rawnc(
        binarydir, rawdir, cacdir, sensorlist, deploymentyaml,
        incremental=True, scisuffix=scisuffix, glidersuffix=glidersuffix)
#%% Ok when this doesnt work go down to TESTING BAD FILES

    # merge individual neetcdf files into single netcdf files *.ebd.nc and *.dbd.nc
    slocum.merge_rawnc(
        rawdir, rawdir, deploymentyaml,
        scisuffix=scisuffix, glidersuffix=glidersuffix)
#%%
    # Make level-1 timeseries netcdf file from th raw files...
    outname = slocum.raw_to_timeseries(
        rawdir, l1tsdir, deploymentyaml,
        profile_filt_time=10, profile_min_time=30)
#%%
if True:
    # make profile netcdf files for ioos gdac...
    ncprocess.extract_timeseries_profiles(outname, profiledir, deploymentyaml)

# make grid of dataset....

outname2 = ncprocess.make_gridfiles(outname, griddir, deploymentyaml)

pgutils.example_gridplot(outname2, './gridplot2.png', ylim=[150, 0],
                         toplot=['potential_temperature', 'salinity',
                                 'oxygen_concentration', 'chlorophyll', 'cdom'])


#%% TESTING BAD FILES This was good at identifying the bad files 
import glob
import xarray as xr

rawdir = "./rawnc"
files = sorted(glob.glob(f"{rawdir}/**/*.nc"))

bad = []
for f in files:
    try:
        ds = xr.open_dataset(f, decode_times=False)
        # common: time dimension empty
        time_len = ds.sizes.get("time", None)
        if time_len == 0 or any(v == 0 for v in ds.sizes.values()):
            bad.append((f, dict(ds.sizes)))
        ds.close()
    except Exception as e:
        bad.append((f, f"FAILED TO OPEN: {e}"))

print("Bad files:")
for item in bad:
    print(item)
#%% This then renames those bad files to singletons
for f, reason in bad:
    newname = f + ".singleton"

    # avoid double-renaming if script is re-run
    if not f.endswith(".singleton") and not os.path.exists(newname):
        os.rename(f, newname)
        print(f"Renamed: {f} -> {newname}")
    else:
        print(f"Skipped (already renamed): {f}")
#%% This is combining everything after cleaning out the bad files 
rawdir= './rawnc_flat/'
slocum.merge_rawnc(
    rawdir, rawdir, deploymentyaml,
    scisuffix=scisuffix, glidersuffix=glidersuffix)
#%% I halved them and it seemed to work way better - no issues
rawdir= './half/'
slocum.merge_rawnc(
    rawdir, rawdir, deploymentyaml,
    scisuffix=scisuffix, glidersuffix=glidersuffix)
#%%
rawdir= './rawnc_flat/'
outname = slocum.raw_to_timeseries(
    rawdir, l1tsdir, deploymentyaml,
    profile_filt_time=100, profile_min_time=300)
#%%
rawdir= './half/'
outname = slocum.raw_to_timeseries(
    rawdir, l1tsdir, deploymentyaml,
    profile_filt_time=100, profile_min_time=300)
#%% Just plotting!
outname2 = ncprocess.make_gridfiles(outname, griddir, deploymentyaml)
pgutils.example_gridplot(outname2, './gridplot2.png', ylim=[200, 0],
                         toplot=['potential_temperature', 'salinity',
                                 'oxygen_concentration', 'chlorophyll', 'cdom'])

#%%
#FROM HERE BELOW IS ANOTHER METHOD I TRIED SPLITING EACH YOYO into files
#Could come in handy if wanting to look at single profiles later
#
outdir = './rawnc/combd'
#%% Using this to test if it worked 
# merge individual neetcdf files into single netcdf files *.ebd.nc and *.dbd.nc
rawdir = './rawnc/0140'
slocum.merge_rawnc(
    rawdir, rawdir, deploymentyaml,
    scisuffix=scisuffix, glidersuffix=glidersuffix)

#%% this runs though each file 
from glob import glob
testpath = './rawnc/*'
fullpath= list(glob(testpath, recursive=True))

for i in fullpath[:-2]:
    print(i)
    slocum.merge_rawnc(
        i, outdir, deploymentyaml,
        scisuffix=scisuffix, glidersuffix=glidersuffix)

#%% Moving them back to a single file of only good ones
import glob
import os
import shutil

rawdir = "./rawnc"
outdir = "./rawnc_flat"

os.makedirs(outdir, exist_ok=True)

files = sorted(glob.glob(f"{rawdir}/**/*.nc", recursive=True))

print(f"Found {len(files)} files")

for f in files:
    fname = os.path.basename(f)
    dest = os.path.join(outdir, fname)

    # handle filename collisions
    if os.path.exists(dest):
        base, ext = os.path.splitext(fname)
        i = 2
        while True:
            new_name = f"{base}__{i}{ext}"
            dest = os.path.join(outdir, new_name)
            if not os.path.exists(dest):
                break
            i += 1

    shutil.copy2(f, dest)
    print(f"Copied: {f} -> {dest}")



