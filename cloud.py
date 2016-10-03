import glob
import sys
import cdms2 as cdms
import numpy as np
import MV2 as MV
import difflib
import pickle
import ENSO_years_piControl as en

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
### Set classic Netcdf (ver 3)
cdms.setNetcdfShuffleFlag(0)
cdms.setNetcdfDeflateFlag(0)
cdms.setNetcdfDeflateLevelFlag(0)

import sc_utilities as sc
import scipy.stats as stats
from scipy.interpolate import interp1d

class CLOUDS():
    def __init__(self):
        self.datasets = ["ISCCP","ISCCP_raw","PATMOSX","PATMOSX_raw"]
        f = cdms.open("OBS/clt_ISCCP_corrected_198301-200912.nc")
        fp = cdms.open("OBS/clt_PATMOSX_corrected_198301-200912.nc")
        
        f_old = cdms.open("OBS/clt_ISCCP_198307-200806.nc")
        fp_old = cdms.open("OBS/clt_PATMOSX_198200-200912.nc")
    
        self.ISCCP = f("clt",time=('1983-7-1','2009-12-31'))
        self.ISCCP = MV.masked_where(np.isnan(self.ISCCP),self.ISCCP)
        cdutil.setTimeBoundsMonthly(self.ISCCP)

        self.PATMOSX = fp("clt",time=('1983-7-1','2009-12-31'))
        self.PATMOSX = MV.masked_where(np.isnan(self.PATMOSX),self.PATMOSX)
        cdutil.setTimeBoundsMonthly(self.PATMOSX)

        self.ISCCP_raw = f_old("clt",time=('1983-7-1','2008-6-31'))
        self.ISCCP_raw = MV.masked_where(np.isnan(self.ISCCP_raw),self.ISCCP_raw)
        cdutil.setTimeBoundsMonthly(self.ISCCP_raw)

        self.PATMOSX_raw = fp_old("clt",time=('1982-1-1','2009-12-31'))
        self.PATMOSX_raw = MV.masked_where(np.isnan(self.PATMOSX_raw),self.PATMOSX_raw)
        cdutil.setTimeBoundsMonthly(self.PATMOSX_raw)
    def get_colors(self,label):
        d={}
        d["ISCCP"]=cm.Blues(.9)
        d["ISCCP_raw"]=cm.Blues(.5)
        d["PATMOSX"]=cm.Reds(.9)
        d["PATMOSX_raw"]=cm.Reds(.5)
        return d[label]
    def plot_seasonal_climatologies(self,season,**kwargs):
        raw=kwargs.pop("raw",False)
        if not raw:
            datasets = ["ISCCP","PATMOSX"]
        else:
             datasets = ["ISCCP_raw","PATMOSX_raw"]
        for dataset in datasets:
            clim=getattr(cdutil,season).climatology(getattr(self,dataset),criteriaarg=(1,None))[0]
            plt.plot(clim.getLatitude()[:],clim.asma(),lw=3,color=self.get_colors(dataset),label=dataset)

    def plot_seasonal_trends(self,season,**kwargs):
        raw=kwargs.pop("raw",False)
        if not raw:
            datasets = ["ISCCP","PATMOSX"]
        else:
             datasets = ["ISCCP_raw","PATMOSX_raw"]
        for dataset in datasets:
            data = getattr(cdutil,season).departures(getattr(self,dataset),criteriaarg=(1,None))
            trends = genutil.statistics.linearregression(data,axis=0,nointercept=1)*120
            lat = data.getLatitude()[:]
            plt.plot(lat,trends.asma(),"-",color=self.get_colors(dataset))#,mec=self.get_colors(dataset))
        plt.axhline(0,ls=":",c="k")

    def get_extrema(self,season,dataset):
        clim=getattr(cdutil,season).climatology(getattr(self,dataset),criteriaarg=(1,None))(latitude=(-58,58))
        smoothed = pf.spatially_smooth(clim)
        xmax,xmin,ymax,ymin = pf.find_all_peaks(smoothed[0],return_maxmin=True)
        return xmax,xmin
    def get_spatially_smoothed(self,season,dataset,sigma=5):
        return pf.spatially_smooth(getattr(cdutil,season)(getattr(self,dataset),criteriaarg=(1,None))(latitude=(-58,58)),sigma=sigma)
    def get_old_indices(self,season,dataset,sigma=5):
        data=self.get_spatially_smoothed(season,dataset,sigma)
        data.season = season
        return pf.thermo_and_dynamic(data,own_bounds=True)
    def regimes(self,season,dataset):
        climatology = getattr(cdutil,season).climatology(getattr(self,dataset),criteriaarg=(1,None))(latitude=(-58,58))
        xclim,yclim = pf.thermo_and_dynamic(pf.spatially_smooth(climatology))
        Sdry_clim=float(xclim[0,1])
        Ndry_clim=float(xclim[0,3])
        SH_clim = float(cdutil.averager(climatology[0](latitude=(-60,Sdry_clim))))
        tropics_clim = float(cdutil.averager(climatology[0](latitude=(Sdry_clim,Ndry_clim))))
        NH_clim = float(cdutil.averager(climatology[0](latitude=(Ndry_clim,60))))
                                            
        x,y=self.get_old_indices(season,dataset)
        Sdry = x[:,1]
        Ndry = x[:,3]
        data = getattr(cdutil,season)(getattr(self,dataset))
        nt = data.shape[0]
        latitudes = MV.zeros((nt,2))
        cloud_fractions = MV.zeros((nt,3))
        for i in range(nt):
            data_timestep=data[i]
            SH_dry=float(x[i,1])
            NH_dry=float(x[i,3])
            if (np.isnan(SH_dry) or np.isnan(NH_dry)):
                latitudes[i]=[np.nan,np.nan]
                cloud_fractions[i]=[np.nan,np.nan,np.nan]
            else:
            

                latitudes[i]=[SH_dry-Sdry_clim,NH_dry-Ndry_clim]
                SH = cdutil.averager(data_timestep(latitude=(-60,SH_dry)))
                tropics = cdutil.averager(data_timestep(latitude=(SH_dry,NH_dry)))
                NH = cdutil.averager(data_timestep(latitude=(NH_dry,60)))
                cloud_fractions[i]=[SH-SH_clim,tropics-tropics_clim,NH-NH_clim]
        latitudes = MV.masked_where(np.isnan(latitudes),latitudes)
        cloud_fractions = MV.masked_where(np.isnan(cloud_fractions),cloud_fractions)
        latitudes.setAxis(0,data.getTime())
        cloud_fractions.setAxis(0,data.getTime())
        return latitudes,cloud_fractions
    def all_regimes(self,dataset):
        months = ["JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"]
        startdata = getattr(self,dataset)
        nt = startdata.shape[0]
        startmonth = startdata.getTime().asComponentTime()[0].month -1
        newmonths = np.roll(months,-startmonth)
        alllats = MV.zeros((nt,2))
        allcf=MV.zeros((nt,3))
        for i in range(12):
            lats,cf=self.regimes(newmonths[i],dataset)
            alllats[i::12]=lats
            allcf[i::12]=cf
        alllats.setAxis(0,startdata.getTime())
        allcf.setAxis(0,startdata.getTime())
        return alllats,allcf
        
        
        
            
            
        
        
        
        




