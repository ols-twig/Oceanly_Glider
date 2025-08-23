# -*- coding: utf-8 -*-

# import pyglider.example_data as pexamp 
# pexamp.get_example_data('./') # already having isssues with this

#%%
# -*- coding: utf-8 -*-
import logging
import os
import xarray
import pandas as pd
import pyglider.ncprocess as ncprocess
import pyglider.slocum as slocum
import pyglider.utils as pgutils
#%%
os.chdir(r'G:\.shortcut-targets-by-id\1Z7PX__ZLrJ6oL2wHtDtiATh5eAh4Cgg3\Oceanly Science\
         Expeditions\Solomon Islands 2025\GLIDER_DATA')
#C:\GitHub\Tonga_Oceanography\Glider\predep_test')
         #C:\GitHub\Tonga_Oceanography\Glider\testing\example-slocum')
#%%
logging.basicConfig(level='INFO') 

binarydir = './raw_sbd_tbd/'
rawdir = './rawnc_sbd_tbd/'
cacdir = './cac/'
sensorlist = './unit_1104_sensors.txt'
deploymentyaml = './deploymentRealtime.yml'
l1tsdir = './L0-timeseries/'
profiledir = './L0-profiles/'
griddir = './L0-gridfiles/'
scisuffix = 'tbd'#'ebd'#
glidersuffix = 'sbd'#'dbd' #
#%%
do_direct = False # changed from true
# only do this for a real run, or something like this
real = False
if real:
    os.system('rsync -av cproof@sfmc.webbresearch.com:/var/opt/sfmc-dockserver/' +
              'stations/dfo/gliders/ ~/gliderdata/slocum_dockserver/')
    os.system('rsync -av ~/gliderdata/slocum_dockserver/rosie_713/from-glider/* ' +
              binarydir)

    os.system('rm ' + rawdir + 'dfo* ' + rawdir + 'TEMP*.nc ' + l1tsdir + '* ' +
              profiledir + '* ' + griddir + '* ')

if do_direct:
    # turn *.sdb and *.tbd into timeseries netcdf files
    outname = slocum.binary_to_timeseries( #changed from binary_to_timeseries
        binarydir, cacdir, l1tsdir, deploymentyaml, search='*.[s|t]bd',
        profile_filt_time=20, profile_min_time=20)
else:
    # turn *.EBD and *.DBD into *.ebd.nc and *.dbd.nc netcdf files.
    slocum.binary_to_rawnc(
        binarydir, rawdir, cacdir, sensorlist, deploymentyaml,
        incremental=True, scisuffix=scisuffix, glidersuffix=glidersuffix)
    
   # merge individual netcdf files into single netcdf files *.ebd.nc and *.dbd.nc
   #%%
    slocum.merge_rawnc(
        rawdir, rawdir, deploymentyaml,
        scisuffix=scisuffix, glidersuffix=glidersuffix)
#%%
    # Make level-1 timeseries netcdf file from th raw files...
    outname = slocum.raw_to_timeseries(
        rawdir, l1tsdir, deploymentyaml,
        profile_filt_time=100, profile_min_time=300)
#%%    
if True:
    # make profile netcdf files for ioos gdac...
    ncprocess.extract_timeseries_profiles(outname, profiledir, deploymentyaml)

#%%
# make grid of dataset....
outname2 = ncprocess.make_gridfiles(outname, griddir, deploymentyaml)
pgutils.example_gridplot(outname2, './gridplot2.png', ylim=[410, 0], # !!changed potential_density to density in function 
                         toplot=['potential_temperature', 'salinity',
                                 'oxygen_concentration', 'chlorophyll', 'cdom'])
#%% Plotting flight data

import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
plt.rcParams["figure.figsize"] = (10,5)
plt.rcParams['figure.dpi'] = 200

ollie=xr.open_dataset('./rawnc/unit_1104rawdbd.nc')
ollie2=xr.open_dataset('./rawnc/unit_1104rawebd.nc')
#./L0-timeseries/unit_1104_Tonga_2024.nc')
# '/rawnc/unit_1104rawdbd.nc')

fig,axes = plt.subplots(4,sharex=True)
plt.rcParams['axes.grid']=True
#fig.tight_layout()#pad = 1.5)

# depth
g = axes[0].scatter(pd.to_datetime(ollie.time.values, unit='s'),
                ollie.m_depth.values, s = 4)    
#axes[0].set_ylim(20,-2)
axes[0].invert_yaxis()
axes[0].set_ylabel('Depth (m)')
axes[0].axhline(y=0, c= 'k', linestyle = '--')

# Angles - Pitch
g = axes[1].scatter(pd.to_datetime(ollie.time.values, unit='s'),
                ollie.m_pitch.values*10, c = 'r', s = 4)
axes[1].set_ylabel('Angle (Pitch - Deg)')
axes[1].axhline(y=0, c= 'k', linestyle = '--')
axes[1].set_ylim(20,-20)

# Dive and Climb Rates
axes[2].scatter(pd.to_datetime(ollie.time.values, unit='s'),
                ollie.m_avg_dive_rate.values, c = 'g', s=4)

axes[2].scatter(pd.to_datetime(ollie.time.values, unit='s'),
                ollie.m_avg_climb_rate.values, s=4, c= 'orange')
axes[2].set_ylabel('Dive and Climb rates')
axes[2].legend(['dive','climb'], loc='upper left')
axes[2].axhline(y=0, c= 'k', linestyle = '--')

axes[3].scatter(pd.to_datetime(ollie2.time.values, unit='s'),
                ollie2.sci_rbrctd_temperature_00.values, s =4, c = 'k')
axes[3].scatter(pd.to_datetime(ollie2.time.values, unit='s'),
                ollie2.sci_rbrctd_salinity_00.values, s =4, c = 'blue')

axes[3].set_ylabel('Sci_temp')

print('Maximum angle =', np.nanmax(ollie.m_pitch.values*10))
print('Minimum angle =', np.nanmin(ollie.m_pitch.values*10))
#%% add temp and salinity


#%% Looking into the L0/L1 files 
#sort out potentiial density nomenalcture


ollie= xarray.open_dataset('./L0-gridfiles//dfo-rosie713-20190615_grid.nc')
list(ollie.keys())

#%%
import xarray as xr
def example_gridplot(filename, outname,
                     toplot=['potential_temperature', 'salinity',
                             'oxygen_concentration'],
                     pdenlevels=np.arange(10, 30, 0.5), dpi=200, ylim=None):
    """
    Smoketest for plotting the netcdf files.
    """

    import matplotlib.pyplot as plt

    ntoplot = len(toplot)
    with xr.open_dataset(filename) as ds:
        fig, axs = plt.subplots(nrows=ntoplot, constrained_layout=True,
                                figsize=(7, 3*ntoplot),
                                sharex=True, sharey=True)
        for ax, vname in zip(axs, toplot):
            ds[vname].plot.pcolormesh(ax=ax)
            (ds['density']-1000).plot.contour(ax=ax, levels=pdenlevels) #changed from potential_denisty
            if ylim:
                ax.set_ylim(ylim)
        fig.savefig(outname, dpi=dpi)
#%%
import cmocean.cm as cmo
scidat = xr.open_dataset('./rawnc/unit_1104rawebd.nc')
#flydat
fig, axs = plt.subplots(1,1)
axs.scatter(pd.to_datetime(scidat['time'].values, unit='s'), scidat['sci_rbrctd_depth_00'],
            c = scidat['sci_rbrctd_temperature_00'], s = 4, cmap = cmo.thermal)
axs.invert_yaxis()



#%%
ollie.close()
ollie2.close()

