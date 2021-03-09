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
workdir='/work/oda/ag15419/tmp/OTPS/ts_check_all/'
inpath='/data/opa/mfs-dev/Med_static/MFS_TPXO_V0/'
outfile=open(workdir+".txt","w")
point_idx=138
point_idy=354
#
for extr_year in range(2019,2020): #2024
  for extr_month in ('01','02','03','04','05','06','07','08','09','10','11','12'):

     outfile=open(workdir+'tpxo2check_'+str(extr_year)+str(extr_month)+'.txt',"w")
     print ('Outfile: ',outfile)

     first_day=1
     for extr_day in ('01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31'):
   
        infiles_pre='tpxoextr_'+str(extr_year)+str(extr_month)+str(extr_day)
   
        # Open input files 
        for file in os.listdir(inpath+'/'+str(extr_year)+'/'):
           if file.startswith(infiles_pre):
              fT1_mod = NC.Dataset(inpath+'/'+str(extr_year)+'/'+file,'r')
              print ('I am working on ',str(extr_year),str(extr_month),str(extr_day))
              print ('Input files template: ',file)
   
              # Read coordinates
              if first_day != 0:
                 latitudes=fT1_mod.variables['nav_lat']
                 longitudes=fT1_mod.variables['nav_lon']
                 print ('Point coordinates: ',latitudes[point_idx,point_idy],longitudes[point_idx,point_idy])

              # Read and write the field
              ssh=fT1_mod.variables['tide_z']
              for idx_time in range(0,len(ssh[:,point_idx,point_idy])):
                 print (ssh[idx_time,point_idx,point_idy],file=outfile)
        first_day=0
