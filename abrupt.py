#!/usr/local/anaconda2/envs/latest/bin/python
import glob
import sys
import cdms2 as cdms
import numpy as np
import MV2 as MV
import difflib
import scipy.stats as stats
global crunchy
import socket
import string
if socket.gethostname().find("crunchy")>=0:
    crunchy = True
else:
    crunchy = False

#import peakfinder as pf
import cdtime,cdutil,genutil
from eofs.cdms import Eof
from eofs.multivariate.cdms import MultivariateEof
import matplotlib.pyplot as plt
import matplotlib.cm as cm
### Set classic Netcdf (ver 3)
cdms.setNetcdfShuffleFlag(0)
cdms.setNetcdfDeflateFlag(0)
cdms.setNetcdfDeflateLevelFlag(0)

from Plotting import *
import CMIP5_tools as cmip5
import DA_tools 
#import peakfinder
#import obs_trends_seasonal_high_low as ot

def regrid_cloud(C):
 
    axes = C.getAxisIds()

 
 
    plev_ax = C.getAxisIds().index("plev")
    Call = MV.sum(C,axis=plev_ax)#(time=('1979-1-1','2005-12-31'))
   # cdutil.setTimeBoundsMonthly(Call)
   
 
    if crunchy:
        fobs = cdms.open("/work/marvel1/CLOUD_SEASONS/cloud-seasons/CLOUD_OBS/clt_ISCCP_corrected_198301-200912.nc")
    else:
        fobs = cdms.open("/Users/kmarvel/Google Drive/CLOUD_SEASONS/cloud-seasons/CLOUD_OBS/clt_ISCCP_corrected_198301-200912.nc")
    the_grid = fobs["clt"].getGrid()
    Cnew = MV.zeros(Call.shape[:2]+the_grid.shape)
    for i in range(7):
        Cold = Call[:,i,:,:]
        Cnew[:,i,:,:] = Cold.regrid(the_grid,regridTool='regrid2')


    Cnew.setAxis(0,C.getTime())
    Cnew.setAxis(1,C.getLevel())
    Cnew.setAxis(2,the_grid.getLatitude())
    Cnew.setAxis(3,the_grid.getLongitude())
    fobs.close()
    Cnew.id = C.id
    return Cnew
if __name__ == "__main__":
    abrupt = cmip5.get_datafiles("abrupt4xCO2","clisccp")
    abrupt_r1 = np.array(abrupt)[np.where(np.array([x.find("r1i")>=0 for x  in abrupt]))]

    for fname in abrupt_r1:
        writename = 'Regridded112017/abrupt4xCO2/'+fname.split("/")[-1].replace("xml","nc")
        fwrite = cdms.open(writename,"w")
        f = cdms.open(fname)
        C = f("clisccp")
        Cnew = regrid_cloud(C)
        f.close()
        fwrite.write(Cnew)
        fwrite.close()

    piC = cmip5.get_datafiles("piControl","clisccp")

    for fname in piC:
        writename = 'Regridded112017/piControl/'+fname.split("/")[-1].replace("xml","nc")
        fwrite = cdms.open(writename,"w")
        f = cdms.open(fname)
        C = f("clisccp")
        Cnew = regrid_cloud(C)

        fwrite.write(Cnew)
        f.close()
        fwrite.close()

    historical = cmip5.get_datafiles("historical","clisccp")
    historical_no_ipsl = np.array(historical)[np.where(np.array([x.find("IPSL")<0 for x in historical]))]
    for fname in historical_no_ipsl:
        writename = 'Regridded112017/historical/'+fname.split("/")[-1].replace("xml","nc")
        fwrite = cdms.open(writename,"w")
        f = cdms.open(fname)
        C = f("clisccp")
        Cnew = regrid_cloud(C)

        fwrite.write(Cnew)
        fwrite.close()
        f.close()
    OnePct = cmip5.get_datafiles("1pctCO2","clisccp")
    OnePct_no_ipsl = np.array(OnePct)[np.where(np.array([x.find("IPSL")<0 for x in OnePct]))]
    for fname in OnePct_no_ipsl:
        writename = 'Regridded112017/1pctCO2/'+fname.split("/")[-1].replace("xml","nc")
        fwrite = cdms.open(writename,"w")
        f = cdms.open(fname)
        C = f("clisccp")
        Cnew = regrid_cloud(C)
        f.close()
        fwrite.write(Cnew)
        fwrite.close()

