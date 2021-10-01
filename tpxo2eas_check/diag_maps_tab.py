#
# Script for HARMONIC ANALYSIS AND POST_PROC
#
# imports
import matplotlib.pyplot as plt
import matplotlib as mpl # Palettes
import numpy as np
import netCDF4 as NC
import os
import sys
import warnings
#
warnings.filterwarnings("ignore") # Avoid warnings
#
from scipy.optimize import curve_fit
from scipy import stats
import collections
import pandas as pd
import csv
import math
import datetime
from datetime import date, timedelta
from datetime import datetime
from operator import itemgetter 
import plotly
from plotly import graph_objects as go # for bar plot
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LogNorm
from operator import itemgetter # to order lists
from statsmodels.distributions.empirical_distribution import ECDF # empirical distribution functions
#
# by AC Goglio (CMCC)
# annachiara.goglio@cmcc.it
#
# Written: 07/2021
# Modified: 08/07/2021
#
# Script to list and plot the correction to the bug in the detiding procedure 
#
#################################################################
# The user should modify the following lines to set his run
#################################################################
# General run parameters:
#---------------------
# Work dir path:
workdir= '/work/oda/ag15419/tmp/detiding_check/points_list' 
# input files:
tpxo_bug='/data/opa/mfs/Med_static/MFS_TPXO_V0/2019/tpxoextr_20190203.nc'
tpxo_debug='/data/opa/mfs-dev/Med_static/MFS_TPXO_V1/2019/tpxoextr_20190203.nc'
# diff obtained with the CDOs cdo eq in1 in2 out
tpxo_cdodiff=workdir+'/diffs_mask.nc'

# Bathymetry and mes mask
model_bathy='/work/oda/ag15419/PHYSW24_DATA/TIDES/DATA0/bathy_meter.nc'
model_meshmask='/work/oda/ag15419/PHYSW24_DATA/TIDES/DATA0/mesh_mask.nc'

########################################################
# DO NOT CHANGE THE CODE BELOW THIS LINE!!!
########################################################
# Open the files:   
bug = NC.Dataset(tpxo_bug,'r')
debug = NC.Dataset(tpxo_debug,'r')
cdodiff = NC.Dataset(tpxo_cdodiff,'r')

field_bug=bug.variables['tide_z'][0,:,:]
field_debug=debug.variables['tide_z'][0,:,:]
field_cdodiff=cdodiff.variables['tide_z'][0,:,:]

lons = debug.variables['nav_lon'][:]
lats = debug.variables['nav_lat'][:]

# Compute the diffs
pydiff =  field_debug - field_bug

# Split AtlanticBox from Med Sea
def which_domain(longitude,latitude):
    # Check if the grid point is inside the whole domain
    if longitude < 36.30 and longitude > -17.2917 and latitude < 45.97917 and latitude > 30.18:
       if longitude > 0.000:
          domain='Med'
       elif longitude < -6.000:
          domain='AtlBox'
       else:
          if latitude < 42.000:
             domain='Med'
          else:
             if longitude < -2.000:
                domain='AtlBox'
             else:
                domain='OUT'
    else:
       domain='OUT'
    return domain


lons_bug=[]
lats_bug=[]

print ('py diff:')
py_idx=1
for yidx in range (0,len(field_debug[:,0])):
   for xidx in range (0,len(field_debug[0,:])):
       if pydiff[yidx,xidx] != 0:
          loc=which_domain(lons[yidx,xidx],lats[yidx,xidx])
          print (py_idx,yidx,xidx,lats[yidx,xidx],lons[yidx,xidx],loc,'\\\\')
          lons_bug.append(lons[yidx,xidx])
          lats_bug.append(lats[yidx,xidx])
          py_idx=py_idx+1

print ('cdo diff:')
cdo_idx=1
for yidx in range (0,len(field_debug[:,0])):
   for xidx in range (0,len(field_debug[0,:])):
       if field_cdodiff[yidx,xidx] == 0:
          loc=which_domain(lons[yidx,xidx],lats[yidx,xidx])
          print (cdo_idx,yidx,xidx,lats[yidx,xidx],lons[yidx,xidx],loc,'\\\\',field_debug[yidx,xidx],field_bug[yidx,xidx])
          cdo_idx=cdo_idx+1

# Plot the map

# Open bathymetry and mes_mask
model3 = NC.Dataset(model_bathy,'r')
vals_bathy=model3.variables['Bathymetry'][:]
model4 = NC.Dataset(model_meshmask,'r')
landmask_mesh=model4.variables['tmask'][0,0,:,:]

lon_0 = lons.mean()
llcrnrlon = lons.min()
llcrnrlon = -10.0
urcrnrlon = lons.max()
lat_0 = lats.mean()
llcrnrlat = lats.min()
urcrnrlat = lats.max()


# Fig
plt.figure(figsize=(20,10))
plt.rc('font', size=13) #  size=12
# Plot Title
plt.title ('Grid points affected by the procedure bug')
m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,resolution='c',projection='merc',lat_0=lat_0,lon_0=lon_0)
xi, yi = m(lons, lats)
# Plot the frame to the map
plt.rcParams["axes.linewidth"]  = 1.25
m.drawparallels(np.arange(30., 46., 5), labels=[1,0,0,0], fontsize=10) 
m.drawmeridians(np.arange(-20., 40., 10), labels=[0,0,0,1], fontsize=10)

# Plot the bathy
cmap = mpl.cm.Blues(np.linspace(0,1,20))
cmap = mpl.colors.ListedColormap(cmap[5:,:-1])
cmap =  cmap.reversed()
cs = m.pcolor(xi,yi,-np.squeeze(vals_bathy),cmap=cmap,vmax=-5000,vmin=0) 
contourf = plt.contourf(xi,yi,landmask_mesh,[0.000,0.999],colors='gray')

# Plot the legend and its label
cbar = m.colorbar(cs, location='right', pad="10%")
cbar.set_label('Bathimetry [m]')

# Add points
for p2plot_idx in range(0,len(lons_bug)) :
    lon_ok=lons_bug[p2plot_idx]
    lat_ok=lats_bug[p2plot_idx]
    xp, yp = m(lon_ok,lat_ok)
    plt.scatter(xp, yp, s=100, alpha=1,marker='o', c='red')

# Save and close 
plt.savefig(workdir+'/bug_points.png')
plt.clf()

