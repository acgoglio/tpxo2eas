# imports
import matplotlib.pyplot as plt
import numpy as np
import os
import netCDF4 as NC
from scipy.optimize import curve_fit
from scipy import stats
import collections
import pandas as pd
import csv
import math
import datetime
from datetime import datetime
from operator import itemgetter
import plotly
#
workdir='/work/oda/ag15419/tmp/new_tpxoextr/' 
inpath='/data/oda/ag15419/TPXO9_DATA/new_tpxoextr/'  # '/data/opa/mfs/Med_static/MFS_TPXO_V1/' #'/data/opa/mfs/Med_static/MFS_TPXO_V0/'

point_idx=281
point_idy=823
###point_lat=39.0625 #39.08161
###point_lon=17.16666603088379 #17.13704

field_name='tide_z'

infile_prename='tpxoextr_'
outfile_prename='mask_tpxo2check'
#
for extr_year in range(2022,2023): 
  for extr_month in ('01','02','03','04','05','06','07','08','09','10','11','12'):
     outfile=open(workdir+outfile_prename+str(extr_year)+str(extr_month)+'.txt',"w")
     print ('Outfile: ',outfile)

     first_day=1
     for extr_day in ('01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31'):
   
        infiles_pre=infile_prename+str(extr_year)+str(extr_month)+str(extr_day)
        print ('Input file ',infiles_pre)

        # Open input files 
        for file in os.listdir(inpath+'/'+str(extr_year)+'/'): 
           if file.startswith(infiles_pre):
              fT1_mod = NC.Dataset(inpath+'/'+str(extr_year)+'/'+file,'r') #+'/'+str(extr_year)+'/'+file,'r')
              print ('I am working on ',str(extr_year),str(extr_month),str(extr_day))
              print ('Input files template: ',file)
   
              # Read coordinates
              ###if first_day != 0:
              if first_day == 1:
                 latitudes=fT1_mod.variables['nav_lat']
                 longitudes=fT1_mod.variables['nav_lon']
                 ###for lat_idx in range(0,len(latitudes[:,0])):
                 ###    for lon_idx in range(0,len(latitudes[0,:])):
                 ###        print ('Prova ',lat_idx,lon_idx,float(latitudes[lat_idx,lon_idx]),float(longitudes[lat_idx,lon_idx]),float(point_lat),float(point_lon))
                 ###        if float(latitudes[lat_idx,lon_idx]) == float(point_lat) and float(longitudes[lat_idx,lon_idx]) == float(point_lon):
                 ###           point_idx=lat_idx
                 ###           point_idy=lon_idx
                 ###           print ('Point coordinates: ',latitudes[point_idx,point_idy],longitudes[point_idx,point_idy])  
                 print ('Point coordinates: ',latitudes[point_idx,point_idy],longitudes[point_idx,point_idy])

              # Read and write the field
              ssh=fT1_mod.variables[field_name]
              for idx_time in range(0,len(ssh[:,point_idx,point_idy])):
                 print (ssh[idx_time,point_idx,point_idy],file=outfile)
        first_day=0
