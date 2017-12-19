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

 
 
    tau_ax = C.getAxisIds().index("tau")
    Call = MV.sum(C,axis=tau_ax)#(time=('1979-1-1','2005-12-31'))
   # cdutil.setTimeBoundsMonthly(Call)
   
 
    if crunchy:
        fobs = cdms.open("/work/marvel1/CLOUD_SEASONS/cloud-seasons/CLOUD_OBS/clt_ISCCP_corrected_198301-200912.nc")
    else:
        fobs = cdms.open("/Users/kmarvel/Google Drive/CLOUD_SEASONS/cloud-seasons/CLOUD_OBS/clt_ISCCP_corrected_198301-200912.nc")
    obs_clt = fobs("clt")
    the_grid = obs_clt.getGrid()
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
def write_data():
    abrupt_files = cmip5.get_datafiles("abrupt4xCO2","clisccp")
    abrupt_r1 = np.array(abrupt_files)[np.where(np.array([x.find("r1i")>=0 for x  in abrupt_files]))]

    for fname in abrupt_r1:
        writename = '/kate/Regridded112017/abrupt4xCO2/'+fname.split("/")[-1].replace("xml","nc")
        fwrite = cdms.open(writename,"w")
        f = cdms.open(fname)
        C = f("clisccp")
        Cnew = regrid_cloud(C)
        
        fwrite.write(Cnew)
        f.close()
        fwrite.close()

    piC = cmip5.get_datafiles("piControl","clisccp")

    for fname in piC:
        writename = '/kate/Regridded112017/piControl/'+fname.split("/")[-1].replace("xml","nc")
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
        writename = '/kate/Regridded112017/historical/'+fname.split("/")[-1].replace("xml","nc")
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
        writename = '/kate/Regridded112017/1pctCO2/'+fname.split("/")[-1].replace("xml","nc")
        fwrite = cdms.open(writename,"w")
        f = cdms.open(fname)
        C = f("clisccp")
        Cnew = regrid_cloud(C)
        f.close()
        fwrite.write(Cnew)
        fwrite.close()

def low_cloud(C):
    plev_ax =C.getAxisIds().index("plev")
    low = MV.sum(C(level=(1000*100,680*100)),axis=plev_ax)
    low.id = "low_clt"
    return low
def mid_cloud(C):
    plev_ax =C.getAxisIds().index("plev")
    mid = MV.sum(C(level=(680*100,440*100)),axis=plev_ax)
    mid.id = "mid_clt"
    return mid
def high_cloud(C):
    plev_ax =C.getAxisIds().index("plev")
    high = MV.sum(C(level=(440*100,0)),axis=plev_ax)
    high.id="high_clt"
    return high

def abrupt_changes(C):
    start = cmip5.start_time(C)
    start20 = start.add(20,cdtime.Years)
    hundred = start20.add(100,cdtime.Years)
    stop = hundred.add(20,cdtime.Years)

    first20 = MV.average(C(time=(start,start20)),axis=0)
    last20 = MV.average(C(time=(hundred,stop)),axis=0)
    return last20-first20


def hist_ensembles():
    hist=sorted(glob.glob("/kate/Regridded112017/historical/*"))
    hist_no_had=[hist[0]]+hist[2:]
    start,stop=cmip5.get_common_timeax(hist_no_had)
    Nmod = len(hist_no_had)
    i=0
    fname = hist_no_had[i]
    f = cdms.open(fname)
    C=f("clisccp",time=(start,stop))
    f.close()
    low_0 = low_cloud(C)
    mid_0 = mid_cloud(C)
    high_0 = high_cloud(C)
    ### Initialise arrays
    LOW = MV.zeros((Nmod,)+low_0.shape)
    MID = MV.zeros((Nmod,)+mid_0.shape)
    HIGH = MV.zeros((Nmod,)+high_0.shape)
    ###
    LOW[i]=low_0
    MID[i]=mid_0
    HIGH[i]=high_0
    for i in range(len(hist_no_had))[1:]:
        fname = hist_no_had[i]
        f = cdms.open(fname)
        C=f("clisccp",time=(start,stop))
        f.close()
        low_i = low_cloud(C)
        mid_i = mid_cloud(C)
        high_i = high_cloud(C)
    
        LOW[i]=low_i
        MID[i]=mid_i
        HIGH[i]=high_i
    modax = cmip5.make_model_axis(hist_no_had)
    LOW.id = "low_clt"
    LOW.setAxis(0,modax)
    LOW.setAxis(1,low_0.getTime())
    LOW.setAxis(2,low_0.getLatitude())
    LOW.setAxis(3,low_0.getLongitude())

    MID.id = "mid_clt"
    MID.setAxis(0,modax)
    MID.setAxis(1,mid_0.getTime())
    MID.setAxis(2,mid_0.getLatitude())
    MID.setAxis(3,mid_0.getLongitude())

    HIGH.id = "high_clt"
    HIGH.setAxis(0,modax)
    HIGH.setAxis(1,high_0.getTime())
    HIGH.setAxis(2,high_0.getLatitude())
    HIGH.setAxis(3,high_0.getLongitude())

    fwrite = cdms.open("/kate/Regridded112017/MMA/cmip5.MMA.historical.r1i1p1.mo.atm.cfMon.clisccp.ver-1.latestX.nc","w")
    fwrite.write(LOW)
    fwrite.write(MID)
    fwrite.write(HIGH)
    fwrite.close()
    return LOW,MID,HIGH
    
    
    
def abrupt_allfiles():
    ab=sorted(glob.glob("/kate/Regridded112017/abrupt4xCO2/*"))
    ab_no_ipsl=ab[:2]+ab[4:]
    Nmod = len(ab_no_ipsl)
    AC = MV.zeros((Nmod,)+(7,72,144))
    i=0
    for fname in ab_no_ipsl:
        f=cdms.open(fname)
        C=f("clisccp")
        test = abrupt_changes(C)
        AC[i]=test
        f.close()
        i+=1
    AC.setAxis(0,cmip5.make_model_axis(ab_no_ipsl))
    for axi in range(4)[1:]:
        AC.setAxis(axi,C.getAxis(axi))
    return AC
    
    
    
def write_zonal(experiment):
    thefiles = sorted(glob.glob("/kate/Regridded112017/"+experiment+"/cmip5.*"))
    writedir = "/kate/Regridded112017/ZONAL/"+experiment+"/"
    for fil in thefiles:
        if fil.find("IPSL")<0:
            f = cdms.open(fil)
            C = f("clisccp")
            Cz = cdutil.averager(C,axis='x')
            high = high_cloud(Cz)
            mid = mid_cloud(Cz)
            low = low_cloud(Cz)
            fw = cdms.open(writedir+fil.split("/")[-1],"w")
            fw.write(high)
            fw.write(mid)
            fw.write(low)
            fw.close()
            f.close()
        
