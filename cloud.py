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
from Plotting import *
import CMIP5_tools as cmip5

class CLOUDS():
    def __init__(self):
        """ Read in cloud and precipitation obs data"""
        self.datasets = ["ISCCP","ISCCP_raw","PATMOSX","PATMOSX_raw"]
        f = cdms.open("OBS/clt_ISCCP_corrected_198301-200912.nc")
        fp = cdms.open("OBS/clt_PATMOSX_corrected_198301-200912.nc")
        
        f_old = cdms.open("OBS/clt_ISCCP_198307-200806.nc")
        fp_old = cdms.open("OBS/clt_PATMOSX_198200-200912.nc")

        fgpcp = cdms.open("OBS/GPCP.precip.mon.mean.nc")
        fcmap = cdms.open("OBS/CMAP.std.precip.mon.mean.nc")
        
    
        self.ISCCP = f("clt",time=('1984-1-1','2009-12-31'))
        self.ISCCP = MV.masked_where(np.isnan(self.ISCCP),self.ISCCP)
        cdutil.setTimeBoundsMonthly(self.ISCCP)

        self.PATMOSX = fp("clt",time=('1984-1-1','2009-12-31'))
        self.PATMOSX = MV.masked_where(np.isnan(self.PATMOSX),self.PATMOSX)
        cdutil.setTimeBoundsMonthly(self.PATMOSX)

        self.ISCCP_raw = f_old("clt",time=('1984-1-1','2008-6-31'))
        self.ISCCP_raw = MV.masked_where(np.isnan(self.ISCCP_raw),self.ISCCP_raw)
        cdutil.setTimeBoundsMonthly(self.ISCCP_raw)

        self.PATMOSX_raw = fp_old("clt",time=('1982-1-1','2009-12-31'))
        self.PATMOSX_raw = MV.masked_where(np.isnan(self.PATMOSX_raw),self.PATMOSX_raw)
        cdutil.setTimeBoundsMonthly(self.PATMOSX_raw)

        self.GPCP = cdutil.averager(fgpcp("precip",time=('1979-1-1','2014-12-31'),latitude=(-90,90)),axis='x')
        cdutil.setTimeBoundsMonthly(self.GPCP)
        self.CMAP = cdutil.averager(fcmap("precip",time=('1979-1-1','2014-12-31'),latitude=(-90,90)),axis='x')
        self.CMAP.setAxis(0,self.GPCP.getTime())
        cdutil.setTimeBoundsMonthly(self.CMAP)
    def trough_trends(self,dataset,smooth = None):
        start='1984-1-1'
        stop = '2008-6-31'
        seasons = ["DJF","MAM","JJA","SON"]
        nhtrends = {}
        shtrends = {}

        data = getattr(self,dataset)(time=(start,stop))
        tax_units = data.getTime().units.split()[0]
        if tax_units == "days":
            fac = 3650 #per day -> per decade
        elif tax_units == "months":
            fac = 120 #per month -> per decade
        elif tax_units == "hours":
            fac = 3650*24 #per month -> per decade
        else:
            print "units not recognized"
            fac=1
        for season in seasons:
            seasondata = getattr(cdutil,season)(data,criteriaarg=(1,None))
            
            if smooth is not None:
                seasondata = pf.spatially_smooth(seasondata,sigma=smooth)
            NH = NH_trough(seasondata)
            SH = SH_trough(seasondata)
            slope,error = genutil.statistics.linearregression(NH,axis=0,nointercept=1,error=1)
           
            nhtrends[season]=[fac*float(slope),fac*float(error),np.ma.average(NH)]

            slope,error = genutil.statistics.linearregression(SH,axis=0,nointercept=1,error=1)
           
            shtrends[season]=[fac*float(slope),fac*float(error),np.ma.average(SH)]
        return nhtrends,shtrends
            
    def get_all_trough_trends(self,smooth=None):
        NH = {}
        SH = {}
        datasets = ["ISCCP","ISCCP_raw","PATMOSX","PATMOSX_raw","GPCP","CMAP"]
        for dataset in datasets:
            nhtrends,shtrends = self.trough_trends(dataset,smooth=smooth)
            NH[dataset] = nhtrends
            SH[dataset] = shtrends
        return NH,SH
    def test_smooth(self):
        self.NH_nosmooth,self.SH_nosmooth = self.get_all_trough_trends()
        self.NH5,self.SH5 = self.get_all_trough_trends(smooth=5)
        self.NH2,self.SH2 = self.get_all_trough_trends(smooth=2)
    
        
        
    def trough_correlations(self,season,smooth = None):
        NH = np.zeros((6,6))
        SH = np.zeros((6,6))
        start='1984-1-1'
        stop = '2008-6-31'
        datasets = ["ISCCP","ISCCP_raw","PATMOSX","PATMOSX_raw","GPCP","CMAP"]
        for i in range(6):
            dataset1 = getattr(self,datasets[i])(time=(start,stop))
            if smooth is not None:
                dataset1 = pf.spatially_smooth(dataset1,sigma=smooth)
            j = 5
            
            while j>i:
                dataset2 = getattr(self,datasets[j])(time=(start,stop))
                if smooth is not None:
                    dataset2 = pf.spatially_smooth(dataset2,sigma=smooth)
                
                SH[i,j]=float(genutil.statistics.correlation(SH_trough(getattr(cdutil,season)(dataset1,criteriaarg=(1,None))),SH_trough(getattr(cdutil,season)(dataset2,criteriaarg=(1,None)))))
                NH[i,j]=float(genutil.statistics.correlation(NH_trough(getattr(cdutil,season)(dataset1,criteriaarg=(1,None))),NH_trough(getattr(cdutil,season)(dataset2,criteriaarg=(1,None)))))
                j -= 1
        NH = MV.masked_where(NH==0,NH)
        SH = MV.masked_where(SH==0,SH)
        plt.subplot(121)
        plt.pcolor(NH,vmin=-1,vmax=1)
        plt.xticks(np.arange(6)+.5,datasets,rotation=90)
        plt.yticks(np.arange(6)+.5,datasets)
        for i in range(6):
            for j in range(6):
                if not NH.mask[i,j]:
                    plt.text(j+.5,i+.5,str(np.round(NH[i,j],2)))
            
        plt.subplot(122)
        plt.pcolor(SH,vmin=-1,vmax=1)
        plt.xticks(np.arange(6)+.5,datasets,rotation=90)
        plt.yticks(np.arange(6)+.5,datasets)
        for i in range(6):
            for j in range(6):
                    if not SH.mask[i,j]:    
                        plt.text(j+.5,i+.5,str(np.round(SH[i,j],2)))
                
        
        return SH,NH
    def amplitudes(self,dataset):
        """ Calculate the amplitude of the annual cycle for each year """
        if ((dataset =="GPCP") or (dataset == "CMAP")):
            start = '1979-1-1'
            stop = '2014-12-31'
        else:
            start = '1984-1-1'
            if dataset == "ISCCP_raw":
                stop = '2007-12-31'
            else:
                stop = '2009-12-31'
        X = getattr(self,dataset)(time=(start,stop))
        R,P = sc.fast_annual_cycle(X)
        return MV.masked_where(np.isnan(R),R)
    def phases(self,dataset):
        """ Calculate phases of the annual cycle for each year """
        start = '1984-1-1'
        if dataset == "ISCCP_raw":
            stop = '2007-12-31'
        else:
            stop = '2009-12-31'
        X = getattr(self,dataset)(time=(start,stop))
        R,P = sc.fast_annual_cycle(X)
        return MV.masked_where(np.isnan(P),P)
    def get_colors(self,label):
        d={}
        d["ISCCP"]=cm.Blues(.9)
        d["ISCCP_raw"]=cm.Blues(.5)
        d["PATMOSX"]=cm.Reds(.9)
        d["PATMOSX_raw"]=cm.Reds(.5)
        d["GPCP"] = cm.PuOr(.1)
        d["CMAP"]=cm.PuOr(.9)
        return d[label]
    def plot_seasonal_climatologies(self,season,**kwargs):
        raw=kwargs.pop("raw",False)
        precip = kwargs.pop("precip",False)
        print precip
        if not raw:
            if precip:
                datasets = ["GPCP","CMAP"]
            else:
                datasets = ["ISCCP","PATMOSX"]
            
        else:
             datasets = ["ISCCP_raw","PATMOSX_raw"]
        for dataset in datasets:
            clim=getattr(cdutil,season).climatology(getattr(self,dataset))[0]
            plt.plot(clim.getLatitude()[:],clim.asma(),lw=3,color=self.get_colors(dataset),label=dataset)

    def plot_seasonal_trends(self,season,**kwargs):
        raw=kwargs.pop("raw",False)
        if not raw:
            datasets = ["ISCCP","PATMOSX"]
        else:
             datasets = ["ISCCP_raw","PATMOSX_raw"]
        for dataset in datasets:
            data = getattr(cdutil,season).departures(getattr(self,dataset))
            trends = genutil.statistics.linearregression(data,axis=0,nointercept=1)*120
            lat = data.getLatitude()[:]
            plt.plot(lat,trends.asma(),"-",color=self.get_colors(dataset))#,mec=self.get_colors(dataset))
        plt.axhline(0,ls=":",c="k")

    def get_amplitude_extrema(self,dataset,smooth=None,raw=False,precip=False,normalize_by_annual = False):
        functions = [SH_peak_stormtrack,SH_eqflank_stormtrack,SH_ITCZ_extent,thermal_equator,NH_ITCZ_extent,NH_eqflank_stormtrack,NH_peak_stormtrack]
        d = {}
        L = len(functions)
        for i in range(L):
            
            func = functions[i]
            
            R = getattr(self,"amplitudes")(dataset)
            if normalize_by_annual:
                R = cmip5.cdms_clone(R/cdutil.YEAR(getattr(self,dataset)),R)
            if smooth is not None:
                R = pf.spatially_smooth(R,sigma=smooth)
            test_i = func(R)
            d[func.__name__]=test_i
        return d
         
    def amplitude_extrema(self,anom=False,smooth=None,raw=False,precip=False,normalize_by_annual = False):
        functions = [SH_peak_stormtrack,SH_eqflank_stormtrack,SH_ITCZ_extent,thermal_equator,NH_ITCZ_extent,NH_eqflank_stormtrack,NH_peak_stormtrack]
        if raw:
            isccp="ISCCP_raw"
            patmos = "PATMOSX_raw"

        elif precip:
            isccp = "GPCP"
            patmos = "CMAP"
        else:
            isccp = "ISCCP"
            patmos = "PATMOSX"
        D = {}
        D[isccp] = {}
        D[patmos] = {}
        L = len(functions)
        for i in range(L):
            plt.subplot(4,2,i+1)
            func = functions[i]
            
            
            R = getattr(self,"amplitudes")(isccp)
            if normalize_by_annual:
                R = cmip5.cdms_clone(R/cdutil.YEAR(getattr(self,isccp)),R)
            if smooth is not None:
                R = pf.spatially_smooth(R,sigma=smooth)
            test_i = func(R)
            D[isccp][func.__name__]=test_i
            if anom:
                test_i = test_i.anom()
            time_plot(test_i,color=self.get_colors(isccp),marker="o",mec=self.get_colors(isccp))
        
            R = getattr(self,"amplitudes")(patmos)
            if normalize_by_annual:
                R = cmip5.cdms_clone(R/cdutil.YEAR(getattr(self,patmos)),R)
            if smooth is not None:
                R = pf.spatially_smooth(R,sigma=smooth)
            test_p = func(R)
            D[patmos][func.__name__]=test_p
            if anom:
                test_p = test_p.anom()
            time_plot(test_p,color=self.get_colors(patmos),marker="o",mec=self.get_colors(patmos))
            plt.title(func.__name__)
           # try:
            #    print "CORRELATION FOR "+func.__name__+ "IS :"+str(genutil.statistics.correlation(test_p,test_i))
            #except:
             #   print "CORRELATION FOR "+func.__name__+ "doesn't exist because records are different lengths"
            if anom:
                plt.ylim(-2.5,2.5)
        return D
    def get_mma_for_fingerprint(self,z,piControl=False):
        functions = [SH_peak_stormtrack,SH_eqflank_stormtrack,SH_ITCZ_extent,thermal_equator,NH_ITCZ_extent,NH_eqflank_stormtrack,NH_peak_stormtrack]
        if piControl:
            X = z.piControl_raw("amp")
        else:
            X = z.amp
        nmod,nt,nlat = X.shape
        nfunc = len(functions)
        test = MV.zeros((nmod,nt,nfunc))+1.e20
        for i in range(nfunc):
            func = functions[i]
            for modeli in range(20):
                try:
                    test[modeli,:,i] = func(X[modeli])
                except:
                    continue
        test = MV.masked_where(test>1.e10,test)
        test.setAxis(0,X.getAxis(0))
        test.setAxis(1,X.getTime())
        funcax = cdms.createAxis(np.arange(nfunc))
        funcax.points = str([x.__name__ for x in functions])
        test.setAxis(-1,funcax)
        return test
    
        

        
    def plot_trends_amplitude_extrema(self,raw=False,smooth=None,precip=False,normalize_by_annual=False):
        labels = []
        functions = [SH_peak_stormtrack,SH_eqflank_stormtrack,SH_ITCZ_extent,thermal_equator,NH_ITCZ_extent,NH_eqflank_stormtrack,NH_peak_stormtrack]
        L = len(functions)
        for i in range(L):
           # plt.subplot(4,2,i+1)
            func = functions[i]
            if raw:
                isccp="ISCCP_raw"
                patmos = "PATMOSX_raw"
            elif precip == True:
                isccp = "GPCP"
                patmos = "CMAP"
            else:
                isccp = "ISCCP"
                patmos = "PATMOSX"
            R = getattr(self,"amplitudes")(isccp)
            if normalize_by_annual:
                R = cmip5.cdms_clone(R/cdutil.YEAR(getattr(self,isccp)),R)
            if smooth is not None:
                R = pf.spatially_smooth(R,sigma=smooth)
            test_i = func(R)
            trend_i,error_i,Pt1,Pt2,Pf1,Pf2 = genutil.statistics.linearregression(test_i,axis=0,nointercept=1,probability=1,error=1)
            plt.errorbar([i],trend_i*120,yerr=error_i*120,color=self.get_colors(isccp),marker="o",mec=self.get_colors(isccp),markersize=10,lw=3,label=isccp)
             #units are months so turn into per decade
           
            
        
            R = getattr(self,"amplitudes")(patmos)
            if normalize_by_annual:
                R = cmip5.cdms_clone(R/cdutil.YEAR(getattr(self,patmos)),R)
            if smooth is not None:
                R = pf.spatially_smooth(R,sigma=smooth)
            test_p = func(R)
            trend_p,error_p,Pt1,Pt2,Pf1,Pf2 = genutil.statistics.linearregression(test_p,axis=0,nointercept=1,probability=1,error=1)
            plt.errorbar([i+.4],trend_p*120,yerr=error_p*120,color=self.get_colors(patmos),marker="o",mec=self.get_colors(patmos),markersize=10,lw=3,label=patmos)
             #units are months so turn into per de
            labels += [func.__name__]
        plt.xticks(np.arange(L)+.2,labels,rotation=90)
            
    def plot_amplitude_and_climatology(self,dataset):
        start = '1984-1-1'
        if dataset == "ISCCP_raw":
            stop = '2007-12-31'
        else:
            stop = '2009-12-31'
        X = getattr(self,dataset)(time=(start,stop))
        lat_plot(MV.average(X,axis=0),color=self.get_colors(dataset),ls="--",lw=4,label = "Climatology ("+dataset+")")
        R = MV.average(self.amplitudes(dataset),axis=0)
        lat_plot(R,color=self.get_colors(dataset),ls="-",lw=4,label = "Annual cycle amplitude ("+dataset+")")
        x,y = pf.find_all_peaks(R(latitude=(-40,40)))
        plt.plot(x,y,"o",color=self.get_colors(dataset),mec=self.get_colors(dataset))
        for thing in x:
            plt.axvline(thing,color=self.get_colors(dataset),ls=":")
            
    def plot_all_amplitudes_and_climatologies(self):
        plt.subplot(211)
        self.plot_amplitude_and_climatology("ISCCP")
        self.plot_amplitude_and_climatology("ISCCP_raw")
        plt.xticks(np.arange(-90,90+15,15))
        latitude_label_ticks(plt.gca())
        plt.xlim(-90,90)
        plt.ylabel("Total Cloud Fraction")
        plt.legend(loc=0,ncol=2,numpoints=1,fontsize=10)
        plt.title("ISCCP")

              
        plt.subplot(212)
        self.plot_amplitude_and_climatology("PATMOSX")
        self.plot_amplitude_and_climatology("PATMOSX_raw")
        plt.xticks(np.arange(-90,90+15,15))
        latitude_label_ticks(plt.gca())
        plt.xlim(-90,90)
        plt.ylabel("Total Cloud Fraction")
        plt.legend(loc=0,ncol=2,numpoints=1,fontsize=10)
        plt.title("PATMOS-x")
    def plot_seasonal_variations(self,dataset,**kwargs):
        if (dataset == "ISCCP") or (dataset == "PATMOSX"):
            variable="Total Cloud Fraction"
        else:
            variable = "Precip (mm/day)"
        linecolor = kwargs.pop("c","k")
        d = self.get_amplitude_extrema(dataset,**kwargs)
        keys = ['NH_peak_stormtrack','NH_ITCZ_extent','SH_ITCZ_extent','SH_peak_stormtrack']
        for k in keys:
            plt.axvline(np.ma.average(d[k]),c=linecolor)
        
        data = getattr(self,dataset)
        months = ["JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"]
        
        colors = [cm.hsv(i/12.) for i in range(12)]
        i=0
        for season in months:
            
            lat_plot(getattr(cdutil,season).climatology(data)[0],lw=3,color=colors[i])
            i+=1
        latitude_label_ticks(plt.gca())
        plt.ylabel(variable)
        plt.title(dataset)
        a = plt.axes([.65, .6, .2, .2])
        patches,texts=a.pie(np.ones(12),startangle=90,colors=colors,labels=  months)
       
            
              
        
        
        
        
import seasonal_cycle_utils as sc
def get_amplitude(dataset):
      start = '1984-1-1'
      stop = '2009-12-31'
      X = dataset(time=(start,stop))
      R,P = sc.fast_annual_cycle(X)
      return MV.masked_where(np.isnan(R),R)
def plot_amplitude(dataset,cmap=cm.RdYlBu):
      R = get_amplitude(dataset)
      nt,nlat = R.shape
      for i in range(nt):
            lat_plot(R[i],c=cmap(i/float(nt)))

def SH_ITCZ_extent(R,value=False):
      nt,nlat = R.shape
      lat_bounds = (-20,0)
      test = []
      climatology = MV.average(R,axis=0)
      Xclim = climatology(latitude=lat_bounds)
      clim_max,clim_min,clim_ymax,clim_ymins=pf.find_all_peaks(Xclim,return_maxmin=True)
      if len(clim_max)>1:
            idx = np.argmax(clim_ymax)
            clim_max = clim_max[idx]
      for i in range(nt):
            X = R[i](latitude=lat_bounds)
            maxes,mins,ymax,ymins=pf.find_all_peaks(X,return_maxmin=True)
            if len(maxes) == 1:
                  if value:
                        test += [ymax[0]]
                  else:
                        test+=[ maxes[0]]
                        
            elif len(maxes) >1:
                  if value:
                        test+=[ymaxes[np.argmin(np.abs(maxes-clim_max))]]
                  else:
                        test+=[maxes[np.argmin(np.abs(maxes-clim_max))]]
            else:
                  test+=[1.e20]
      test = MV.masked_where(np.array(test)>1.e10,test)
      test.setAxis(0,R.getTime())
      return test
                  
def NH_ITCZ_extent(R,value=False):
      nt,nlat = R.shape
      lat_bounds = (0,20)
      test = []
      climatology = MV.average(R,axis=0)
      Xclim = climatology(latitude=lat_bounds)
      clim_max,clim_min,clim_ymax,clim_ymins=pf.find_all_peaks(Xclim,return_maxmin=True)
      if len(clim_max)>1:
            idx = np.argmax(clim_ymax)
            clim_max = clim_max[idx]
      for i in range(nt):
            X = R[i](latitude=lat_bounds)
            maxes,mins,ymax,ymins=pf.find_all_peaks(X,return_maxmin=True)
            if len(maxes) == 1:
                  if value:
                        test += [ymax[0]]
                  else:
                        test+=[ maxes[0]]
                        
            elif len(maxes) >1:
                  if value:
                        test+=[ymaxes[np.argmin(np.abs(maxes-clim_max))]]
                  else:
                        test+=[maxes[np.argmin(np.abs(maxes-clim_max))]]
            else:
                  test+=[1.e20]
      test = MV.masked_where(np.array(test)>1.e10,test)
      test.setAxis(0,R.getTime())
      return test

def NH_peak_stormtrack(R,value=False):
      nt,nlat = R.shape
      lat_bounds = (30,50)
      test = []
      climatology = MV.average(R,axis=0)
      Xclim = climatology(latitude=lat_bounds)
      clim_max,clim_min,clim_ymax,clim_ymins=pf.find_all_peaks(Xclim,return_maxmin=True)
      if len(clim_max)>1:
            idx = np.argmax(clim_ymax)
            clim_max = clim_max[idx]
      for i in range(nt):
            X = R[i](latitude=lat_bounds)
            maxes,mins,ymax,ymins=pf.find_all_peaks(X,return_maxmin=True)
            if len(maxes) == 1:
                  if value:
                        test += [ymax[0]]
                  else:
                        test+=[ maxes[0]]
                        
            elif len(maxes) >1:
                  if value:
                        test+=[ymaxes[np.argmin(np.abs(maxes-clim_max))]]
                  else:
                        test+=[maxes[np.argmin(np.abs(maxes-clim_max))]]
            else:
                  test+=[1.e20]
      test = MV.masked_where(np.array(test)>1.e10,test)
      test.setAxis(0,R.getTime())
      return test
def SH_peak_stormtrack(R,value=False):
      nt,nlat = R.shape
      lat_bounds = (-50,-30)
      test = []
      climatology = MV.average(R,axis=0)
      Xclim = climatology(latitude=lat_bounds)
      clim_max,clim_min,clim_ymax,clim_ymins=pf.find_all_peaks(Xclim,return_maxmin=True)
      if len(clim_max)>1:
            idx = np.argmax(clim_ymax)
            clim_max = clim_max[idx]
      for i in range(nt):
            X = R[i](latitude=lat_bounds)
            maxes,mins,ymax,ymins=pf.find_all_peaks(X,return_maxmin=True)
            if len(maxes) == 1:
                  if value:
                        test += [ymax[0]]
                  else:
                        test+=[ maxes[0]]
                        
            elif len(maxes) >1:
                  if value:
                        test+=[ymaxes[np.argmin(np.abs(maxes-clim_max))]]
                  else:
                        test+=[maxes[np.argmin(np.abs(maxes-clim_max))]]
            else:
                  test+=[1.e20]
      test = MV.masked_where(np.array(test)>1.e10,test)
      test.setAxis(0,R.getTime())
      return test
def thermal_equator(R,value=False):
      nt,nlat = R.shape
      lat_bounds = (-10,10)
      test = []
      climatology = MV.average(R,axis=0)
      Xclim = climatology(latitude=lat_bounds)
      clim_max,clim_min,clim_ymax,clim_ymins=pf.find_all_peaks(Xclim,return_maxmin=True)
      if len(clim_min)>1:
            idx = np.argmin(clim_ymins)
            clim_min = clim_min[idx]
      
      for i in range(nt):
            X = R[i](latitude=lat_bounds)
            maxes,mins,ymax,ymins=pf.find_all_peaks(X,return_maxmin=True)
            if len(mins) == 1:
                  if value:
                        test += [ymins[0]]
                  else:
                        test+=[ mins[0]]
                        
            elif len(mins) >1:
                  if value:
                        test+=[ymins[np.argmin(np.abs(mins-clim_min))]]
                  else:
                        test+=[mins[np.argmin(np.abs(mins-clim_min))]]
            else:
                  test+=[1.e20]
      test = MV.masked_where(np.array(test)>1.e10,test)
      test.setAxis(0,R.getTime())
      return test             
def NH_eqflank_stormtrack(R,value=False):
      nt,nlat = R.shape
      lat_bounds = (10,40)
      test = []
      climatology = MV.average(R,axis=0)
      Xclim = climatology(latitude=lat_bounds)
      clim_max,clim_min,clim_ymax,clim_ymins=pf.find_all_peaks(Xclim,return_maxmin=True)
      if len(clim_min)>1:
            idx = np.argmin(clim_ymins)
            clim_min = clim_min[idx]
      
      for i in range(nt):
            X = R[i](latitude=lat_bounds)
            maxes,mins,ymax,ymins=pf.find_all_peaks(X,return_maxmin=True)
            if len(mins) == 1:
                  if value:
                        test += [ymins[0]]
                  else:
                        test+=[ mins[0]]
                        
            elif len(mins) >1:
                  if value:
                        test+=[ymins[np.argmin(np.abs(mins-clim_min))]]
                  else:
                        test+=[mins[np.argmin(np.abs(mins-clim_min))]]
            else:
                  test+=[1.e20]
      test = MV.masked_where(np.array(test)>1.e10,test)
      test.setAxis(0,R.getTime())
      return test          
def SH_eqflank_stormtrack(R,value=False):
      nt,nlat = R.shape
      lat_bounds = (-40,-10)
      test = []
      climatology = MV.average(R,axis=0)
      Xclim = climatology(latitude=lat_bounds)
      clim_max,clim_min,clim_ymax,clim_ymins=pf.find_all_peaks(Xclim,return_maxmin=True)
      if len(clim_min)>1:
            idx = np.argmin(clim_ymins)
            clim_min = clim_min[idx]
      
      for i in range(nt):
            X = R[i](latitude=lat_bounds)
            maxes,mins,ymax,ymins=pf.find_all_peaks(X,return_maxmin=True)
            if len(mins) == 1:
                  if value:
                        test += [ymins[0]]
                  else:
                        test+=[ mins[0]]
                        
            elif len(mins) >1:
                  if value:
                        test+=[ymins[np.argmin(np.abs(mins-clim_min))]]
                  else:
                        test+=[mins[np.argmin(np.abs(mins-clim_min))]]
            else:
                  test+=[1.e20]
      test = MV.masked_where(np.array(test)>1.e10,test)
      test.setAxis(0,R.getTime())
      return test          
def get_locations(dataset):
      #Climatological locations of max and min
      R = get_amplitude(dataset)
      Rsmooth=pf.spatially_smooth(R,sigma=2)
      Rsmooth = MV.masked_where(np.abs(Rsmooth)>1.e10,Rsmooth)
      
      nt,nlat = R.shape
      PEAKS = MV.zeros((nt,3))+1.e20
      TROUGHS = MV.zeros((nt,2))+1.e20
      PEAKS[:,0]=SH_ITCZ_extent(R)
      PEAKS[:,1]=NH_ITCZ_extent(R)
      PEAKS[:,2]=NH_peak_stormtrack(R)
      TROUGHS[:,0]=SH_eqflank_stormtrack(R)
      TROUGHS[:,1]=NH_eqflank_stormtrack(R)
                  
      PEAKS = MV.masked_where(PEAKS>1.e10,PEAKS)
      TROUGHS = MV.masked_where(TROUGHS>1.e10,TROUGHS)
      PEAKS.setAxis(0,R.getTime())
      TROUGHS.setAxis(0,R.getTime())
      return PEAKS,TROUGHS
            
def SH_trough(R,value=False):
      nt,nlat = R.shape
      lat_bounds = (-40,-10)
      test = []
      climatology = MV.average(R,axis=0)
      Xclim = climatology(latitude=lat_bounds)
      clim_max,clim_min,clim_ymax,clim_ymins=pf.find_all_peaks(Xclim,return_maxmin=True)
      if len(clim_min)>1:
            clim_min = np.min(clim_min)
      
      for i in range(nt):
            X = R[i](latitude=lat_bounds)
            if len(np.where(X.mask)[0]) == len(X):
                  test+=[1.e20]
                  continue
            maxes,mins,ymax,ymins=pf.find_all_peaks(X,return_maxmin=True)
            if len(mins) == 1:
                  if value:
                        test += [ymins[0]]
                  else:
                        test+=[ mins[0]]
            elif len(mins) >1:
                  if value:
                        test += [ymins[np.argmin(np.abs(mins-clim_min))]]
                  else:
                        test+=[mins[np.argmin(np.abs(mins-clim_min))]]
            else:
                  test+=[1.e20]
      test = MV.masked_where(np.array(test)>1.e10,test)
      test.setAxis(0,R.getTime())
      return test                      
      
def central_max(R,value=False):
      nt,nlat = R.shape
      lat_bounds = (-20,20)
      test = []
      climatology = MV.average(R,axis=0)
      Xclim = climatology(latitude=lat_bounds)
      clim_max,clim_min,clim_ymax,clim_ymins=pf.find_all_peaks(Xclim,return_maxmin=True)
      if len(clim_min)>1:
            clim_min = np.min(clim_min)
      
      for i in range(nt):
            X = R[i](latitude=lat_bounds)
            if len(np.where(X.mask)[0]) == len(X):
                  test+=[1.e20]
                  continue
            maxes,mins,ymax,ymins=pf.find_all_peaks(X,return_maxmin=True)
            if len(maxes) == 1:
                  if value:
                        test += [ymax[0]]
                  else:
                        test+=[ maxes[0]]
            elif len(mins) >1:
                  if value:
                        test += [ymax[np.argmin(np.abs(maxes-clim_max))]]
                  else:
                        test+=[maxes[np.argmin(np.abs(maxes-clim_max))]]
            else:
                  test+=[1.e20]
      test = MV.masked_where(np.array(test)>1.e10,test)
      test.setAxis(0,R.getTime())
      return test         
def NH_trough(R,value=False):
      nt,nlat = R.shape
      lat_bounds = (10,40)
      test = []
      climatology = MV.average(R,axis=0)
      Xclim = climatology(latitude=lat_bounds)
      clim_max,clim_min,clim_ymax,clim_ymins=pf.find_all_peaks(Xclim,return_maxmin=True)
      if len(clim_min)>1:
            clim_min = np.min(clim_min)
      
      for i in range(nt):
            X = R[i](latitude=lat_bounds)
            if len(np.where(X.mask)[0]) == len(X):
                  test+=[1.e20]
                  continue
            maxes,mins,ymax,ymins=pf.find_all_peaks(X,return_maxmin=True)
            if len(mins) == 1:
                  if value:
                        test += [ymins[0]]
                  else:
                        test+=[ mins[0]]
            elif len(mins) >1:
                  if value:
                        test += [ymins[np.argmin(np.abs(mins-clim_min))]]
                  else:
                        test+=[mins[np.argmin(np.abs(mins-clim_min))]]
            else:
                  test+=[1.e20]
      test = MV.masked_where(np.array(test)>1.e10,test)
      test.setAxis(0,R.getTime())
      return test                      
                  
def plot_trough_trends(NH,SH):
    datasets = ["ISCCP","ISCCP_raw","PATMOSX","PATMOSX_raw","GPCP","CMAP"]
    seasons =["DJF","MAM","JJA","SON"]
    ax1 = plt.subplot(111)
    plt.xlim(-.5,4.5)
    plt.ylim(-60,60)
    seasoncolors = {}
    seasoncolors["JJA"]=cm.Reds(.8)
    seasoncolors["DJF"]=cm.Blues(.8)
    seasoncolors["MAM"]=cm.Greens(.8)
    seasoncolors["SON"]=cm.Oranges(.5)
    diffs = np.linspace(0,.6,6)
    for dataset in datasets:
        spacing = diffs[datasets.index(dataset)]
        for season in seasons:
            i = seasons.index(season)
            x = spacing+i
            trend,err,mean = NH[dataset][season]
            y = mean
            dx = 0
            dy = 10*trend
            width=0.2
            ax1.add_patch(patches.Arrow(x,y,dx,dy,width,color=seasoncolors[season]))
            
            ax1.text(x,0,dataset,rotation="vertical",verticalalignment="center",fontsize=10)
            trend,err,mean = SH[dataset][season]
            y = mean
            dx = 0
            dy = 10*trend
            width=0.2
            ax1.add_patch(patches.Arrow(x,y,dx,dy,width,color=seasoncolors[season]))

def latplot_check(c,dataset,smooth=None):
    counter =1
    for season in ["DJF","MAM","JJA","SON"]:
        plt.subplot(2,2,counter)
        data=getattr(cdutil,season)(getattr(c,dataset),criteriaarg=(1,None))
        if smooth is not None:
            data = pf.spatially_smooth(data,sigma=smooth)
        nt,nlat = data.shape
        for i in range(nt):
            lat_plot(data[i],c=cm.RdYlBu(i/float(nt)),lw=1)
        plt.title(season)
        counter +=1
           
   
sys.path.append("/Users/kmarvel/Google Drive/SEASONAL/")
from Final_code import ZonalData
import DA_tools as DA
def noise_histograms(control=None):
    if control is None:
        z = ZonalData(restrict_poles=(-60,60))
        c = CLOUDS()
        test = c.get_mma_for_fingerprint(z)
        clim = MV.average(test(time=('1979-1-1','2014-12-31')),axis=1)
        
        ocontrol = c.get_mma_for_fingerprint(z,piControl=True)
        control_anoms = cmip5.cdms_clone(ocontrol.asma()-clim.asma()[:,np.newaxis,:],ocontrol)
        control = DA.concatenate_this(control_anoms)
    N = MV.zeros((91,7))
    count = 0
    for i in np.arange(10,101,1):
        for j in range(7):
            N[count,j] = np.ma.std(DA.get_slopes(control[:,j],i))
        count +=1
    return N
        
    
    
        
        




