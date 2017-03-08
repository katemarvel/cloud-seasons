#!/usr/local/anaconda2/envs/latest/bin/python
import glob
import sys
import cdms2 as cdms
import numpy as np
import MV2 as MV

#Working remotely?
global crunchy
import socket
if socket.gethostname().find("crunchy")>=0:
    crunchy = True
else:
    crunchy = False

import cdutil,genutil

### Import plotting tools
import matplotlib.pyplot as plt
import matplotlib.cm as cm

### Set classic Netcdf (ver 3)
cdms.setNetcdfShuffleFlag(0)
cdms.setNetcdfDeflateFlag(0)
cdms.setNetcdfDeflateLevelFlag(0)

### Import modules I wrote that (git repo python-utils)
if crunchy:
    sys.path.append("/work/marvel1/python-utils")
else:
    sys.path.append("~/Google Drive/python-utils")
from Plotting import *
import CMIP5_tools as cmip5

def HadleyCell(VA):
    if len(VA.getTime())< 140*12:
        return
    plev = VA.getLevel()[:]
    plev_units = VA.getLevel().units
    if plev_units == "Pa":
        fac = 100
    else:
        fac = 1
    i500 = np.where(plev/100 == 500)[0] #mid troposphere = 500hPa
    fgrid = cdms.open("~/precip.mon.mean.nc")
    the_grid = fgrid["precip"].getGrid()
    data= VA.regrid(the_grid,regridTool='regrid2')[:140*12]
    VAz = cdutil.averager(data,axis='x')
    dp = np.diff(plev)
    dp = np.append(dp,0-plev[-1]) #Top of atmosphere = pressure level 0
    mass_streamf = MV.array(np.cumsum(VAz*dp[np.newaxis,:,np.newaxis],axis=1))[:,i500,:]
    mass_streamf.setAxis(0,VAz.getTime())
    mass_streamf.setAxis(1,VAz.getLatitude())
    mass_streamf.id = "mass_streamfunction"
    return mass_streamf

    
def OnePercentHadleyMMA():
    prefix = "/work/cmip5/1pctCO2/atm/mo/"
    direc = prefix+"va/"
    mma = cmip5.get_ensemble(direc,"va",func=HadleyCell,verbose=True)
    fw = cdms.open("ZonalMeanData/mass_streamfunction.1pctCO2.zonalmean.nc","w")
    mma.id=variable
    fw.write(mma)
    fw.close()
    return mma

    
def HadleyCell_historical(VA):
    VA = VA(time=('1860-1-1','2005-12-31'))
    plev = VA.getLevel()[:]
    plev_units = VA.getLevel().units
    if plev_units == "Pa":
        fac = 100
    else:
        fac = 1
    i500 = np.where(plev/100 == 500)[0] #mid troposphere = 500hPa
    fgrid = cdms.open("~/precip.mon.mean.nc")
    the_grid = fgrid["precip"].getGrid()
    data= VA.regrid(the_grid,regridTool='regrid2')
    VAz = cdutil.averager(data,axis='x')
    dp = np.diff(plev)
    dp = np.append(dp,0-plev[-1]) #Top of atmosphere = pressure level 0
    mass_streamf = MV.array(np.cumsum(VAz*dp[np.newaxis,:,np.newaxis],axis=1))[:,i500,:]
    mass_streamf.setAxis(0,VAz.getTime())
    mass_streamf.setAxis(1,VAz.getLatitude())
    mass_streamf.id = "mass_streamfunction"
    return mass_streamf

    
def historicalHadleyMMA():
    prefix = "/work/cmip5/1pctCO2/atm/mo/"
    direc = prefix+"va/"
    mma = cmip5.get_ensemble(direc,"va",func=HadleyCell_historical,verbose=True)
    fw = cdms.open("ZonalMeanData/mass_streamfunction.historical.zonalmean.nc","w")
    mma.id=variable
    fw.write(mma)
    fw.close()
    return mma    
    


def historical_rcp85_zonal_mean(x):
    """ Function to regrid and calculate zonal mean of spliced historical-rcp85 runs """
    start = '1900-1-1'
    stop = '2100-1-1'
    fgrid = cdms.open("~/precip.mon.mean.nc")
    the_grid = fgrid["precip"].getGrid()
    data = x.regrid(the_grid,regridTool='regrid2')(time=(start,stop))
    return cdutil.averager(data,axis='x')


def historical_zonal_mean(x):
    """ Function to regrid and calculate zonal mean of historical runs """
    start = '1860-1-1'
    stop = '2006-1-1'
    
    fgrid = cdms.open("~/precip.mon.mean.nc")
    the_grid = fgrid["precip"].getGrid()
    data = x.regrid(the_grid,regridTool='regrid2')(time=(start,stop))
    return cdutil.averager(data,axis='x')
def rcp85_zonal_mean(x):
    """ Function to regrid and calculate zonal mean of rcp85 runs """
    start = '2006-1-1'
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

def historical_rcp85_ensemble(variable):
    prefix = "/work/cmip5/historical-rcp85/atm/mo/"
    direc = prefix+variable+"/"
    mma = cmip5.get_ensemble(direc,variable,func=historical_rcp85_zonal_mean)
    fw = cdms.open("ZonalMeanData/"+variable+".historical_rcp85_ensemble.zonalmean.nc","w")
    mma.id=variable
    fw.write(mma)
    fw.close()
    return mma

def historical_ensemble(variable):
    prefix = "/work/cmip5/historical/atm/mo/"
    direc = prefix+variable+"/"
    mma = cmip5.get_ensemble(direc,variable,func=historical_zonal_mean)
    fw = cdms.open("ZonalMeanData/"+variable+".historical.zonalmean.nc","w")
    mma.id=variable
    fw.write(mma)
    fw.close()
    return mma


def rcp85_ensemble(variable):
    prefix = "/work/cmip5/rcp85/atm/mo/"
    direc = prefix+variable+"/"
    mma = cmip5.get_ensemble(direc,variable,func=rcp85_zonal_mean)
    fw = cdms.open("ZonalMeanData/"+variable+".rcp85.zonalmean.nc","w")
    mma.id=variable
    fw.write(mma)
    fw.close()
    return mma


def make_splice_dictionary(historical,rcp85):
    hmodels = eval(historical.getAxis(0).models)
    rcp85models = eval(rcp85.getAxis(0).models)
    get_id = lambda x: x.split(".")[1]+"."+x.split(".")[3]
    h_id = map(get_id,hmodels)
    r_id = map(get_id,rmodels)
    direc = {}
    for i in range(len(h_id)):
        cor=np.where(np.array(r_id) == h_id[i])[0]
        if len(cor) == 1:
            direc[hmodels[i]] = rcp85models[cor]
    return direc
def check_splice_times(direc):
    timediff = []
    for k in sorted(direc.keys()):
        fh = cdms.open(k)
        fr = cdms.open(direc[k])
        hunits = fh["clt"].getTime().units
        rstart = fr["clt"].getTime().asComponentTime()[0]
        hstop = fh["clt"].getTime().asComponentTime()[-1]
        timediff += [rstart.torel(hunits).value-hstop.torel(hunits).value]
        fh.close()
        fr.close()
    return timediff

def splice(historical,rcp85):
    hmodels = eval(historical.getAxis(0).models)
    rcp85models = eval(rcp85.getAxis(0).models)
    hi = []
    ri = []
    get_id = lambda x: x.split(".")[1]+"."+x.split(".")[3]
    h_id = map(get_id,hmodels)
    r_id = map(get_id,rmodels)
    direc = {}
    for i in range(len(h_id)):
        cor=np.where(np.array(r_id) == h_id[i])[0]
        if len(cor) == 1:
            hi += [i]
            ri += [int(cor)]
    hm_good = np.array(hmodels)[hi]
    rm_good = np.array(rcp85models)[ri]
    L = len(hm_good)
    splice_models = [hm_good[i]+" + "+ rm_good[i]]
    spliceax = cdms.createAxis(np.arange(L))
    spliceax.id = "model"
    spliceax.models = str(splice_models)
    historical_trunc = MV.array(historical.asma()[hi])
    historical_trunc.setAxis(1,historical.getAxis(1))
    historical_trunc.setAxis(2,historical.getAxis(2))
    rcp85_trunc = MV.array(rcp85.asma()[ri])
    rcp85_trunc.setAxis(1,rcp85.getAxis(1))
    rcp85_trunc.setAxis(2,rcp85.getAxis(2))
    spliced = MV.concatenate([historical_trunc,rcp85_trunc],axis=1)
    spliced.id = "clt"
    fw = cdms.open("ZonalMeanData/clt.historical-rcp85_ensemble.zonalmean.nc","w")
    fw.write(spliced)
    fw.close()
    

def historicalMisc_ensemble(variable,forcing="AA"):
    prefix = "/work/cmip5/historicalMisc/atm/mo/"
    path = prefix+variable+"/"
    hm = cmip5.HistoricalMisc()
    search_strings = getattr(hm,forcing)
    ens = cmip5.get_ensemble(path,variable,search_string = search_strings[0],func=historical_zonal_mean)
    models = eval(ens.getAxis(0).models)
    for i in range(len(hm.AA))[1:]:
        print hm.AA[i]
        try:
            ens_new = cmip5.get_ensemble(path,variable,search_string = search_strings[i],func=historical_zonal_mean)
        except:
            continue
        if ens_new.shape[1:] == ens.shape[1:]:
            models += eval(ens_new.getAxis(0).models)
            ens = MV.concatenate([ens,ens_new],axis=0,axisattributes={"models":str(models)})
        else:
            print "PROBLEM: SHAPE IS "+str(ens_new.shape)
    ens.id=variable
    fw = cdms.open("ZonalMeanData/"+variable+"."+forcing+"_ensemble.zonalmean.nc","w")
    fw.write(ens)
    fw.close()
    return ens


def clisccp_zonal_mean(data):
    start = '1980-1-1'
    stop = '2006-1-1'
    Tot = MV.sum(data,axis=1)
    x= MV.sum(Tot,axis=1)
    fgrid = cdms.open("~/precip.mon.mean.nc")
    the_grid = fgrid["precip"].getGrid()
    data = x.regrid(the_grid,regridTool='regrid2')(time=(start,stop))
    return cdutil.averager(data,axis='x')



def clisccp_mma():
    variable="clisccp"
    prefix = "/work/cmip5/historical/atm/mo/"
    direc = prefix+variable+"/"
    mma = cmip5.multimodel_average(direc,variable,func=clisccp_zonal_mean)
    fw = cdms.open("ZonalMeanData/"+variable+".clisccp.zonalmean.nc","w")
    mma.id=variable
    fw.write(mma)
    fw.close()
    return mma

def clisccp_ensemble():
    variable="clisccp"
    prefix = "/work/cmip5/historical/atm/mo/"
    
    direc = prefix+variable+"/"
    mma = cmip5.get_ensemble(direc,variable,func=clisccp_zonal_mean)
    fw = cdms.open("ZonalMeanData/"+variable+".clisccp_ensemble.zonalmean.nc","w")
    mma.id=variable
    fw.write(mma)
    fw.close()
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
    fw = cdms.open("ZonalMeanData/"+variable+".1pctCO2.zonalmean.nc","w")
    mma.id=variable
    fw.write(mma)
    fw.close()
    return mma

    
def piC_zonal_mean(x):
  
    fgrid = cdms.open("~/precip.mon.mean.nc")
    the_grid = fgrid["precip"].getGrid()
    data = x.regrid(the_grid,regridTool='regrid2')[:250*12]
    if len(data.getTime())!= 250*12:
        return
    else:
        return cdutil.averager(data,axis='x')

def piControl_ensemble(variable):
    prefix = "/work/cmip5/piControl/atm/mo/"
    direc = prefix+variable+"/"
    mma = cmip5.get_ensemble(direc,variable,func=piC_zonal_mean)
    fw = cdms.open("ZonalMeanData/"+variable+".piControl.zonalmean.nc","w")
    mma.id=variable
    fw.write(mma)
    fw.close()
    return mma

def amip_zonal_mean(x):
    start = '1979-1-1'
    stop = '2008-12-31'
    fgrid = cdms.open("~/precip.mon.mean.nc")
    the_grid = fgrid["precip"].getGrid()
    data = x.regrid(the_grid,regridTool='regrid2')(time=(start,stop))
    return cdutil.averager(data,axis='x')

def amip_mma(variable):
    prefix = "/work/cmip5/amip/atm/mo/"
    direc = prefix+variable+"/"
    mma = cmip5.multimodel_average(direc,variable,func=amip_zonal_mean)
    fw = cdms.open("ZonalMeanData/"+variable+".amip_mma.zonalmean.nc","w")
    mma.id=variable
    fw.write(mma)
    fw.close()
    return mma

def amip_ensemble(variable):
    prefix = "/work/cmip5/amip/atm/mo/"
    direc = prefix+variable+"/"
    mma = cmip5.get_ensemble(direc,variable,func=amip_zonal_mean)
    fw = cdms.open("ZonalMeanData/"+variable+".amip_ensemble.zonalmean.nc","w")
    mma.id=variable
    fw.write(mma)
    fw.close()
    return mma

if __name__ == "__main__":
    historicalHadleyMMA()
    OnePercentHadleyMMA()
    # import Index3 as i3
    # for variable in ["clwvi"]:
    #   historical_ensemble(variable)
      
    #   rcp85_ensemble(variable)
      
    #   #for forcing in ["AA","Oz","Vl","Ant","LU","Sl"]:
    #    # historicalMisc_ensemble(variable,forcing)
    # i3.write_model_data("historical_ensemble",write_pr=True)
    # i3.write_model_data("rcp85_ensemble",write_pr=True)

