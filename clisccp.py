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
import peakfinder
import obs_trends_seasonal_high_low as ot

def low_cloud(clisccp,thresh=440,nyears=200):
 
    axes = clisccp.getAxisIds()

 
 
    plev_ax = clisccp.getAxisIds().index("plev")
    low = MV.sum(clisccp(level=(1000*100,thresh*100)),axis=plev_ax)#(time=('1979-1-1','2005-12-31'))
    cdutil.setTimeBoundsMonthly(low)
    annual_cycle_removed= cdutil.ANNUALCYCLE.departures(low)[0:nyears*12]
 
    if crunchy:
        fobs = cdms.open("/work/marvel1/CLOUD_SEASONS/cloud-seasons/CLOUD_OBS/clt_ISCCP_corrected_198301-200912.nc")
    else:
        fobs = cdms.open("/Users/kmarvel/Google Drive/CLOUD_SEASONS/cloud-seasons/CLOUD_OBS/clt_ISCCP_corrected_198301-200912.nc")
    the_grid = fobs["clt"].getGrid()
    annual_cycle_removed_regrid = annual_cycle_removed.regrid(the_grid,regridTool='regrid2')
    fobs.close()
    return annual_cycle_removed_regrid





def high_cloud(clisccp,thresh=440,nyears=200):
 
    axes = clisccp.getAxisIds()

 
 
    plev_ax = clisccp.getAxisIds().index("plev")
    high = MV.sum(clisccp(level=(thresh*100,0)),axis=plev_ax)#(time=('1979-1-1','2005-12-31'))
    cdutil.setTimeBoundsMonthly(high)
    annual_cycle_removed= cdutil.ANNUALCYCLE.departures(high)[0:nyears*12]
 
    if crunchy:
        fobs = cdms.open("/work/marvel1/CLOUD_SEASONS/cloud-seasons/CLOUD_OBS/clt_ISCCP_corrected_198301-200912.nc")
    else:
        fobs = cdms.open("/Users/kmarvel/Google Drive/CLOUD_SEASONS/cloud-seasons/CLOUD_OBS/clt_ISCCP_corrected_198301-200912.nc")
    the_grid = fobs["clt"].getGrid()
    annual_cycle_removed_regrid = annual_cycle_removed.regrid(the_grid,regridTool='regrid2')
    fobs.close()
    return annual_cycle_removed_regrid

def write_1pctCO2_cloud():
    direc = "/kate/cl_regrid_isccp/1pctCO2/XML/"
    search_string = "*r1*xml"
    HIGH = cmip5.get_ensemble(direc,"pgrid_cl",search_string=search_string,func=high_cloud)
    HIGH.id="high_clt"
    fh = cdms.open("cmip5.1pctCO2.ensemble.high_clt.nc","w")
    fh.write(HIGH)
    fh.close()
    
    LOW = cmip5.get_ensemble(direc,"pgrid_cl",search_string=search_string,func=low_cloud)
    LOW.id="low_mid_clt"
    fl = cdms.open("cmip5.1pctCO2.ensemble.lowmid_clt.nc","w")
    fl.write(LOW)
    fl.close()
    return HIGH,LOW

def write_piControl_cloud(high=False):
    direc = "/kate/cl_regrid_isccp/piControl/XML/"
    search_string = "*r1*xml"

    if high:
        HIGH = cmip5.get_ensemble(direc,"pgrid_cl",search_string=search_string,func=high_cloud)
        HIGH.id="high_clt"
        fh = cdms.open("cmip5.piControl.ensemble.high_clt.nc","w")
        fh.write(HIGH)
        fh.close()
    
    LOW = cmip5.get_ensemble(direc,"pgrid_cl",search_string=search_string,func=low_cloud)
    LOW.id="low_mid_clt"
    fl = cdms.open("cmip5.piControl.ensemble.lowmid_clt.nc","w")
    fl.write(LOW)
    fl.close()
    #return HIGH,LOW




### Historical
def low_cloud_historical(clisccp,thresh=440,start='1983-7-1',stop='2005-12-31'):
 
    axes = clisccp.getAxisIds()

 
 
    plev_ax = clisccp.getAxisIds().index("plev")
    low = MV.sum(clisccp(level=(1000*100,thresh*100)),axis=plev_ax)(time=(start,stop))
    cdutil.setTimeBoundsMonthly(low)
    annual_cycle_removed= cdutil.ANNUALCYCLE.departures(low)
 
    if crunchy:
        fobs = cdms.open("/work/marvel1/CLOUD_SEASONS/cloud-seasons/CLOUD_OBS/clt_ISCCP_corrected_198301-200912.nc")
    else:
        fobs = cdms.open("/Users/kmarvel/Google Drive/CLOUD_SEASONS/cloud-seasons/CLOUD_OBS/clt_ISCCP_corrected_198301-200912.nc")
    the_grid = fobs["clt"].getGrid()
    annual_cycle_removed_regrid = annual_cycle_removed.regrid(the_grid,regridTool='regrid2')
    fobs.close()
    return annual_cycle_removed_regrid

def high_cloud_historical(clisccp,thresh=440,start='1983-7-1',stop='2005-12-31'):
 
    axes = clisccp.getAxisIds()

 
 
    plev_ax = clisccp.getAxisIds().index("plev")
    high = MV.sum(clisccp(level=(thresh*100,0)),axis=plev_ax)(time=(start,stop))
    cdutil.setTimeBoundsMonthly(high)
    annual_cycle_removed= cdutil.ANNUALCYCLE.departures(high)
 
    if crunchy:
        fobs = cdms.open("/work/marvel1/CLOUD_SEASONS/cloud-seasons/CLOUD_OBS/clt_ISCCP_corrected_198301-200912.nc")
    else:
        fobs = cdms.open("/Users/kmarvel/Google Drive/CLOUD_SEASONS/cloud-seasons/CLOUD_OBS/clt_ISCCP_corrected_198301-200912.nc")
    the_grid = fobs["clt"].getGrid()
    annual_cycle_removed_regrid = annual_cycle_removed.regrid(the_grid,regridTool='regrid2')
    fobs.close()
    return annual_cycle_removed_regrid

def write_historical_cloud():
    direc = "/kate/cl_regrid_isccp/historical-rcp85/XML/"
    search_string = "*r1*xml"
    HIGH = cmip5.get_ensemble(direc,"pgrid_cl",search_string=search_string,func=high_cloud_historical,stop='2099-12-31')
    HIGH.id="high_clt"
    fh = cdms.open("cmip5.historical-rcp85.ensemble.high_clt.nc","w")
    fh.write(HIGH)
    fh.close()
    
    LOW = cmip5.get_ensemble(direc,"pgrid_cl",search_string=search_string,func=low_cloud_historical,stop='2099-12-31')
    LOW.id="low_mid_clt"
    fl = cdms.open("cmip5.historical-rcp85.ensemble.lowmid_clt.nc","w")
    fl.write(LOW)
    fl.close()
   # return HIGH,LOW

def write_Nat_cloud():
    direc = "/kate/cl_regrid_isccp/historicalNat/XML/"
    search_string = "*r1*xml"
    HIGH = cmip5.get_ensemble(direc,"pgrid_cl",search_string=search_string,func=high_cloud_historical,stop='2005-12-31')
    HIGH.id="high_clt"
    fh = cdms.open("cmip5.historicalNat.ensemble.high_clt.nc","w")
    fh.write(HIGH)
    fh.close()
    
    LOW = cmip5.get_ensemble(direc,"pgrid_cl",search_string=search_string,func=low_cloud_historical,stop='2005-12-31')
    LOW.id="low_mid_clt"
    fl = cdms.open("cmip5.historicalNat.ensemble.lowmid_clt.nc","w")
    fl.write(LOW)
    fl.close()
    #return HIGH,LOW

if __name__ == "__main__":
    #write_piControl_cloud(high=False)
    write_Nat_cloud()
    write_historical_cloud()
    
