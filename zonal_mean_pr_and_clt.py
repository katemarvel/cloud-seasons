import glob
import sys
import cdms2 as cdms
import numpy as np
import MV2 as MV
import difflib
import pickle
#import ENSO_years_piControl as en

global crunchy
import socket
if socket.gethostname().find("crunchy")>=0:
    crunchy = True
else:
    crunchy = False

import peakfinder as pf
import cdtime,cdutil,genutil
from eofs.cdms import Eof
from eofs.multivariate.cdms import MultivariateEof
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as patches
### Set classic Netcdf (ver 3)
cdms.setNetcdfShuffleFlag(0)
cdms.setNetcdfDeflateFlag(0)
cdms.setNetcdfDeflateLevelFlag(0)


import scipy.stats as stats
from scipy.interpolate import interp1d

if crunchy:
    sys.path.append("/work/marvel1/python-utils")
else:
    sys.path.append("~/Google Drive/python-utils")
from Plotting import *
import CMIP5_tools as cmip5

def historical_rcp85_zonal_mean(x):
    start = '1900-1-1'
    stop = '2100-1-1'
    fgrid = cdms.open("~/precip.mon.mean.nc")
    the_grid = fgrid["precip"].getGrid()
    data = x.regrid(the_grid,regridTool='regrid2')(time=(start,stop))
    return cdutil.averager(data,axis='x')

def historical_rcp85_mma(variable):
    prefix = "/work/cmip5/historical-rcp85/atm/mo/"
    direc = prefix+variable+"/"
    mma = cmip5.multimodel_average(direc,variable,func=historical_rcp85_zonal_mean)
    return mma


    
def OnePct_zonal_mean(x):
  
    fgrid = cdms.open("~/precip.mon.mean.nc")
    the_grid = fgrid["precip"].getGrid()
    data = x.regrid(the_grid,regridTool='regrid2')[:140*12]
    if len(data.getTime())!= 140*12:
        return
    else:
        return cdutil.averager(data,axis='x')

def OnePct_mma(variable):
    prefix = "/work/cmip5/1pctCO2/atm/mo/"
    direc = prefix+variable+"/"
    mma = cmip5.multimodel_average(direc,variable,func=OnePct_zonal_mean,verbose=True)
    fw = cdms.open(variable+".1pctCO2.zonalmean.nc","w")
    mma.id=variable
    fw.write(mma)
    fw.close()
    return mma


