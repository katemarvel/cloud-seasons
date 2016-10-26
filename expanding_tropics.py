import glob
import sys
import cdms2 as cdms
import numpy as np
import MV2 as MV
import difflib
import pickle
import datetime

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

class TropicTracker():
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
        self.datasets = ["ISCCP","ISCCP_raw","PATMOSX","PATMOSX_raw","GPCP","CMAP"]
        self.functions = [SH_trough,center_max,NH_trough]
        self.seasons = ["DJF","MAM","JJA","SON"]

        
    def get_season_data(self,dataset,season,smooth=None):
        data = getattr(self,dataset)
        seasondata = getattr(cdutil,season)(data,criteriaarg=(1,None))
        if smooth is not None:
            seasondata = pf.spatially_smooth(seasondata, sigma=smooth)
        return seasondata

    def write_indices_to_file(self,smooth=None):
        direc = "ExpandingTropicIndices/"
        
        for func in self.functions:
            latwriteprefix = direc+"Smooth."+str(smooth)+"."+func.__name__+"_lat"
            valwriteprefix = direc+"Smooth."+str(smooth)+"."+func.__name__+"_value"
            for season in self.seasons:
                latstub = latwriteprefix+"."+season
                valstub = valwriteprefix+"."+season
                
                for dataset in self.datasets:
                    latfile = cdms.open(latstub+"."+dataset+".nc","w")
                    data = self.get_season_data(dataset,season,smooth=smooth)
                    towrite = func(data)
                    towrite.id = "lat"
                    towrite.info = "Written by KM on "+str(datetime.datetime.today())
                    latfile.write(towrite)
                    latfile.close()

                    valfile = cdms.open(valstub+"."+dataset+".nc","w")
                    data = self.get_season_data(dataset,season,smooth=smooth)
                    towrite = func(data,value=True)
                    if dataset in ["GPCP","CMAP"]:
                        towrite.id = "pr"
                    else:
                        towrite.id="clt"
                    towrite.info = "Written by KM on "+str(datetime.datetime.today())
                    valfile.write(towrite)
                    valfile.close()

    def read_indices_from_file(self,smooth=None):
        direc = "ExpandingTropicIndices/"
        
        for func in self.functions:
            latwriteprefix = direc+"Smooth."+str(smooth)+"."+func.__name__+"_lat"
            valwriteprefix = direc+"Smooth."+str(smooth)+"."+func.__name__+"_value"

            setattr(self,func.__name__+"_lat",{})
            setattr(self,func.__name__+"_value",{})
            
            for season in self.seasons:
                getattr(self,func.__name__+"_lat")[season]={}
                getattr(self,func.__name__+"_value")[season]={}
                latstub = latwriteprefix+"."+season
                valstub = valwriteprefix+"."+season
                
                for dataset in self.datasets:
                    latfile = cdms.open(latstub+"."+dataset+".nc","r")
                    getattr(self,func.__name__+"_lat")[season][dataset]=latfile("lat")
                    latfile.close()
                         
                    valfile = cdms.open(valstub+"."+dataset+".nc","r")
                    
                    if dataset in ["GPCP","CMAP"]:
                        variable = "pr"
                    else:
                        variable ="clt"
                    getattr(self,func.__name__+"_value")[season][dataset] = valfile(variable)

                    valfile.close()
                    
                    
    def indices(self,smooth=None):
        for func in self.functions:
            setattr(self,func.__name__+"_lat",{})
            setattr(self,func.__name__+"_value",{})
            for season in self.seasons:
               
                getattr(self,func.__name__+"_lat")[season]={}
                getattr(self,func.__name__+"_value")[season]={}
                for dataset in self.datasets:
          
                    data = self.get_season_data(dataset,season,smooth=smooth)
                    getattr(self,func.__name__+"_lat")[season][dataset] = func(data)
                    getattr(self,func.__name__+"_value")[season][dataset] = func(data,value=True)

    
            
    def plot_indices(self,season,anom=False,value=False):
        if not hasattr(self,"SH_trough_lat"):
            self.indices()
        plt.subplot(313)
        for dataset in self.datasets:
            if value:
                to_plot = self.SH_trough_value[season][dataset]
            else:
                to_plot = self.SH_trough_lat[season][dataset]
            if anom:
                to_plot = to_plot.anom()
            time_plot(to_plot,color=self.get_colors(dataset),mec=self.get_colors(dataset),lw=3,marker = "o")
        plt.legend(loc=0,fontsize=8,ncol=3)
        plt.title("Southern Hemisphere")
        ax2 = plt.subplot(312)
        for dataset in self.datasets:
            if value:
                to_plot = self.center_max_value[season][dataset]
            else:
                to_plot = self.center_max_lat[season][dataset]
            if anom:
                to_plot = to_plot.anom()
            time_plot(to_plot,color=self.get_colors(dataset),mec=self.get_colors(dataset),lw=3,marker="o")
        plt.title("Central Peak")
        ax3 = plt.subplot(311)
        for dataset in self.datasets:
            if value:
                to_plot = self.NH_trough_value[season][dataset]
            else:
                to_plot = self.NH_trough_lat[season][dataset]
            if anom:
                to_plot = to_plot.anom()
            time_plot(to_plot,color=self.get_colors(dataset),lw=3,mec=self.get_colors(dataset),label=dataset,marker="o")
        plt.title("Northern Hemisphere")

    def plot_trends(self,smooth):
        self.read_indices_from_file(smooth=smooth)
        diffs = np.linspace(0,.6,6)
        ax1 = plt.subplot(111)
        plt.xlim(-.2,3.8)
        plt.ylim(-50,50)
        seasoncolors = {}
        seasoncolors["JJA"]=cm.Reds(.8)
        seasoncolors["DJF"]=cm.Blues(.8)
        seasoncolors["MAM"]=cm.Greens(.8)
        seasoncolors["SON"]=cm.Oranges(.5)
        for dataset in self.datasets:
            spacing = diffs[self.datasets.index(dataset)]
            data = getattr(self,dataset)
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
            for season in self.seasons:
                i = self.seasons.index(season)
                x = spacing+i
                data = self.SH_trough_lat[season][dataset]
                mean = np.ma.average(data)
                trend = float(genutil.statistics.linearregression(data,axis=0,nointercept=1))*fac
                
                y = mean
                dx = 0.
                dy = 10*trend
                width=0.2
                
                #ax1.add_patch(patches.Arrow(x,y,dx,dy,width=width,color=seasoncolors[season]))
                ax1.add_patch(patches.Arrow(x,y,dx,dy,width=width,color=self.get_colors(dataset)))
                if season == "DJF":
                    ax1.text(x,0,dataset,rotation="vertical",verticalalignment="center",fontsize=10)

                data = self.center_max_lat[season][dataset]
                mean = np.ma.average(data)
                trend = float(genutil.statistics.linearregression(data,axis=0,nointercept=1))*fac
                y = mean
                dx = 0.
                dy = 10*trend
                width=0.2
                #ax1.add_patch(patches.Arrow(x,y,dx,dy,width,color=seasoncolors[season]))
                ax1.add_patch(patches.Arrow(x,y,dx,dy,width=width,color=self.get_colors(dataset)))
                              
                data = self.NH_trough_lat[season][dataset]
                mean = np.ma.average(data)
                trend = float(genutil.statistics.linearregression(data,axis=0,nointercept=1))*fac
                y = mean
                dx = 0.
                dy = 10*trend
                width=0.2
                #ax1.add_patch(patches.Arrow(x,y,dx,dy,width,color=seasoncolors[season]))
                ax1.add_patch(patches.Arrow(x,y,dx,dy,width=width,color=self.get_colors(dataset)))
        plt.xticks(np.arange(4)+.3,self.seasons)
        [plt.axvline(.8+i,color=cm.Greys(.5),lw=4) for i in range(4)]

        

        
        
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


       
  
def SH_trough(R,value=False):
      nt,nlat = R.shape
      lat_bounds = (-50,-10)
      test = []
      climatology = MV.average(R,axis=0)
      Xclim = climatology(latitude=lat_bounds)
      clim_max,clim_min,clim_ymax,clim_ymins=pf.find_all_peaks(Xclim,return_maxmin=True)
      if len(clim_min)>1:
            clim_min = np.min(clim_min)
      elif len(clim_min) == 0:
            clim_min = -23.43 #tropic of Capricorn
      
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
      
def center_max(R,value=False):
      nt,nlat = R.shape
      lat_bounds = (-20,20)
      test = []
      climatology = MV.average(R,axis=0)
      Xclim = climatology(latitude=lat_bounds)
      clim_max,clim_min,clim_ymax,clim_ymins=pf.find_all_peaks(Xclim,return_maxmin=True)
      if len(clim_max)>1:
            clim_max = np.max(clim_max)
      elif len(clim_max) ==0:
            clim_max =0. #equator
      
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
            elif len(maxes) >1:
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
      lat_bounds = (10,50)
      test = []
      climatology = MV.average(R,axis=0)
      Xclim = climatology(latitude=lat_bounds)
      clim_max,clim_min,clim_ymax,clim_ymins=pf.find_all_peaks(Xclim,return_maxmin=True)
      if len(clim_min)>1:
            clim_min = np.min(clim_min)
      elif len(clim_min) == 0:
            clim_min = 23.43 #tropic of Cancer
      
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
           
def write_model_data(smooth=5):
    direc = "ExpandingTropicIndices/historical-rcp85/"
    functions = [SH_trough,center_max,NH_trough]
    fclt = cdms.open("clt_hist85.mma.nc")
    cltdata = fclt("clt")
    cdutil.setTimeBoundsMonthly(cltdata)
    fclt.close()
    fpr = cdms.open("pr_hist85.mma.nc")
    prdata = fpr("pr")
    cdutil.setTimeBoundsMonthly(prdata)
    fpr.close()
    
    for func in functions:
        latwriteprefix = direc+"Smooth."+str(smooth)+"."+func.__name__+"_lat"
        valwriteprefix = direc+"Smooth."+str(smooth)+"."+func.__name__+"_value"
        for season in ["DJF","MAM","JJA","SON"]:
            latstub = latwriteprefix+"."+season
            valstub = valwriteprefix+"."+season
            
            seasonal_data_pr = getattr(cdutil,season)(prdata,criteriaarg=(1,None))
            nmod,nt,nlat = seasonal_data_pr.shape
            indices_pr = MV.zeros((nmod,nt))+1.e20
            valindices_pr = MV.zeros((nmod,nt))+1.e20
            
            for model_i in range(nmod):
                if (len(np.where(seasonal_data_pr[model_i][1].mask)[0]) == nt):
                    continue
                model_data_pr = seasonal_data_pr[model_i]
                if smooth is not None:
                    model_data_pr = pf.spatially_smooth(model_data_pr,sigma=smooth)

                indices_pr[model_i] = func(model_data_pr)
                valindices_pr[model_i] = func(model_data_pr,value=True)
            indices_pr = MV.masked_where(indices_pr>1.e10,indices_pr)
            valindices_pr = MV.masked_where(valindices_pr>1.e10,valindices_pr)
            
            indices_pr.id = "lat"
            indices_pr.setAxisList(seasonal_data_pr.getAxisList()[:-1])
            valindices_pr.id = "pr"
            valindices_pr.setAxisList(seasonal_data_pr.getAxisList()[:-1])
            
            fwrite_pr = cdms.open(latstub+"historical-rcp85.pr.nc","w")
            fwrite_pr.write(indices_pr)
            fwrite_pr.close()

            vfwrite_pr = cdms.open(valstub+"historical-rcp85.pr.nc","w")
            vfwrite_pr.write(valindices_pr)
            vfwrite_pr.close()
            ###### CLT ######
            seasonal_data_clt = getattr(cdutil,season)(cltdata,criteriaarg=(1,None))
            nmod,nt,nlat = seasonal_data_clt.shape
            indices_clt = MV.zeros((nmod,nt))+1.e20
            valindices_clt = MV.zeros((nmod,nt))+1.e20
            
            for model_i in range(nmod):
                if (len(np.where(seasonal_data_clt[model_i][1].mask)[0]) == nt):
                    continue
                model_data_clt = seasonal_data_clt[model_i]
                if smooth is not None:
                    model_data_clt = pf.spatially_smooth(model_data_clt,sigma=smooth)

                indices_clt[model_i] = func(model_data_clt)
                valindices_clt[model_i] = func(model_data_clt,value=True)
            indices_clt = MV.masked_where(indices_clt>1.e10,indices_clt)
            valindices_clt = MV.masked_where(valindices_clt>1.e10,valindices_clt)
            
            indices_clt.id = "lat"
            indices_clt.setAxisList(seasonal_data_clt.getAxisList()[:-1])
            valindices_clt.id = "clt"
            valindices_clt.setAxisList(seasonal_data_clt.getAxisList()[:-1])
            
            fwrite_clt = cdms.open(latstub+"historical-rcp85.clt.nc","w")
            fwrite_clt.write(indices_clt)
            fwrite_clt.close()

            vfwrite_clt = cdms.open(valstub+"historical-rcp85.clt.nc","w")
            vfwrite_clt.write(valindices_clt)
            vfwrite_clt.close()
            

def check_model(model,season,typ):
    """ Show the time series and locations of troughs/peak for a given model"""
    
    funcs =  ["SH_trough","center_max","NH_trough"]
    fz = cdms.open(typ+"_hist85.mma.nc")
    zonal = fz(typ)
    models = np.array(eval(zonal.getAxis(0).models))
    i = np.where(models == model)[0][0]
    modzonal = zonal[i]
    cdutil.setTimeBoundsMonthly(modzonal)
    seasonal_data_unsmooth=getattr(cdutil,season)(modzonal,criteriaarg=(1,None))
    nt,nlat = seasonal_data_unsmooth.shape
    [lat_plot(seasonal_data_unsmooth[j],color=cm.RdYlBu(j/float(nt))) for j in range(nt)]
    plt.figure()
    seasonal_data = modzonal_smooth = pf.spatially_smooth(seasonal_data_unsmooth,sigma=5)
    
    
    
    ax1 = plt.subplot(221)
    [lat_plot(seasonal_data[j],color=cm.RdYlBu(j/float(nt))) for j in range(nt)]
    fz.close()
  
    for L in range(3):
        axn = plt.subplot(2,2,L+2)
        func = funcs[L]
        flat = cdms.open("ExpandingTropicIndices/historical-rcp85/Smooth.5."+func+"_lat."+season+"historical-rcp85."+typ+".nc")
        lats = flat("lat")
        print lats.shape
        fval = cdms.open("ExpandingTropicIndices/historical-rcp85/Smooth.5."+func+"_value."+season+"historical-rcp85."+typ+".nc")
        vals = fval(typ)
        x = lats[i].asma()
        y = vals[i].asma()
        t = cmip5.get_plottable_time(lats[i])
        axn.plot(t,x)
        for j in range(nt):
            ax1.plot([x[j]],[y[j]],"o",color=cm.RdYlBu(j/float(nt)))
        fval.close()
        flat.close()
    
    
                        
def get_mma_array(season,typ,xy="lat"):
    
    funcs =  ["SH_trough","center_max","NH_trough"]
    func = funcs[0]
    f = cdms.open("ExpandingTropicIndices/historical-rcp85/Smooth.5."+func+"_lat."+season+"historical-rcp85."+typ+".nc")
    shape = f["lat"].shape
    f.close()
    S = MV.zeros(shape+(3,))+1.e20
    for i in range(3):
        func = funcs[i]
        f = cdms.open("ExpandingTropicIndices/historical-rcp85/Smooth.5."+func+"_lat."+season+"historical-rcp85."+typ+".nc")
        data = f("lat")
        S[:,:,i]=data
        modax = data.getAxis(0)
        tax = data.getTime()
        f.close()
    S = MV.masked_where(S>1.e10,S)
    S.setAxis(0,modax)
    S.setAxis(1,tax)
    return S
def get_seasoncolor(season):
    seasoncolors = {}
    seasoncolors["JJA"]=cm.Reds(.8)
    seasoncolors["DJF"]=cm.Blues(.8)
    seasoncolors["MAM"]=cm.Greens(.8)
    seasoncolors["SON"]=cm.Oranges(.5)
    return seasoncolors[season]
def plot_time_series(season,typ,anom=True,include_individual=True):
   
    S = get_mma_array(season,typ)
    if anom:
        S = cmip5.time_anomalies(S,start='1984-1-1',stop='2008-12-31')
    nmod = S.shape[0]
    mma = MV.average(S,axis=0)[1:-1]
    for i in range(3):
        plt.subplot(1,3,i+1)
        if include_individual:
            [time_plot(S[m,:,i],color="k",alpha=.5) for m in range(nmod)]
        time_plot(mma[:,i],color=get_seasoncolor(season),lw=3)
    
    
def plot_eof(season,typ,**kwargs):
    S = get_mma_array(season,typ)
    Sclim = MV.average(MV.average(S(time=('1984-1-1','2008-12-31')),axis=1),axis=0)
    Sa = cmip5.time_anomalies(S,start='1984-1-1',stop='2008-12-31')
    mma = MV.average(Sa,axis=0)[1:-1]
    solver = Eof(mma)
    plt.subplot(211)
    fac = cmip5.get_orientation(solver)
    eof1 = fac*solver.eofs()[0]
    plt.plot(Sclim.asma(),eof1.asma(),"o-",**kwargs)
    plt.axhline(0,c="k",ls=":")

    plt.subplot(212)
    time_plot(fac*solver.pcs()[:,0],**kwargs)

    print solver.varianceFraction()[0]
    
def barplot_trends(typ,start='1900-1-1',stop='2100-1-1'):
    if typ == "pr":
        negsup = "Driest"
        possup="Wettest"
    else:
        negsup="Clearest"
        possup = "Cloudiest"
    
    ax1 = plt.subplot(313)
    ax1.set_title("SH "+negsup+" Latitude")
    ax2 = plt.subplot(312)
    ax2.set_title(possup+" Latitude")
    ax3 = plt.subplot(311)
    ax3.set_title("NH "+negsup+" Latitude")
    axes = [ax1,ax2,ax3]
    counter = 0.
    
    for season in ["DJF","MAM","JJA","SON"]:
        
        S = get_mma_array(season,typ)
        xbar = np.arange(S.shape[0])+counter
        for i in range(3):
            seasonlats = S[:,:,i](time=(start,stop))
            trends = genutil.statistics.linearregression(seasonlats,axis=1,nointercept=1)*3650
            axes[i].bar(xbar,trends.asma(),color=get_seasoncolor(season),ec=get_seasoncolor(season),width=.2,label=season)
        counter += .2
    ax1.set_xticks(np.arange(S.shape[0])+.4)
    ax1.set_xticklabels(eval(S.getAxis(0).models),rotation=90)
    ax2.set_xticks([])
    ax3.set_xticks([])
    ax1.legend(loc=0,ncol=4)
        
    
        



