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
import DA_tools as da

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
        direc = "ExpandingTropicIndices/OBS/"
        
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
        direc = "OBS/"
        
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
        diffs = np.linspace(0,.6,4)
        ax1 = plt.subplot(111)
        plt.xlim(-.2,3.8)
        plt.ylim(-50,50)
        seasoncolors = {}
        seasoncolors["JJA"]=cm.Reds(.8)
        seasoncolors["DJF"]=cm.Blues(.8)
        seasoncolors["MAM"]=cm.Greens(.8)
        seasoncolors["SON"]=cm.Oranges(.5)
        for dataset in self.datasets[:-2]:
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
      if R.getTime() is not None:
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
      if R.getTime() is not None:
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
      if R.getTime() is not None:
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
           
def write_model_data(experiment, smooth=5,verbose=False,write_pr = False):
    direc = "ExpandingTropicIndices/"+experiment+"/"
    if experiment not in os.listdir("ExpandingTropicIndices/"):
        os.system("mkdir "+direc)
    functions = [SH_trough,center_max,NH_trough]
    fclt = cdms.open("ZonalMeanData/clt."+experiment+".zonalmean.nc")
    cltdata = fclt("clt")
    cdutil.setTimeBoundsMonthly(cltdata)
    fclt.close()
    if write_pr:
        fpr = cdms.open("ZonalMeanData/pr."+experiment+".zonalmean.nc")
    
        prdata = fpr("pr")
        cdutil.setTimeBoundsMonthly(prdata)
        fpr.close()
    
    for func in functions:
        if verbose:
            print func.__name__
            
        latwriteprefix = direc+"Smooth."+str(smooth)+"."+func.__name__+"_lat"
        valwriteprefix = direc+"Smooth."+str(smooth)+"."+func.__name__+"_value"
        for season in ["DJF","MAM","JJA","SON"]:
            if verbose:
                print season
            latstub = latwriteprefix+"."+season
            valstub = valwriteprefix+"."+season

            if write_pr:
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

                fwrite_pr = cdms.open(latstub+"."+experiment+".pr.nc","w")
                fwrite_pr.write(indices_pr)
                fwrite_pr.close()

                vfwrite_pr = cdms.open(valstub+"."+experiment+".pr.nc","w")
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
            
            fwrite_clt = cdms.open(latstub+"."+experiment+".clt.nc","w")
            fwrite_clt.write(indices_clt)
            fwrite_clt.close()

            vfwrite_clt = cdms.open(valstub+"."+experiment+".clt.nc","w")
           
            vfwrite_clt.write(valindices_clt)
            vfwrite_clt.close()
            
import matplotlib as mpl
from mpl_toolkits.axes_grid.inset_locator import inset_axes
def check_model(model,season,typ,experiment="historical-rcp85"):
    """ Show the time series and locations of troughs/peak for a given model"""
    
    funcs =  ["SH_trough","center_max","NH_trough"]
    fz = cdms.open("ZonalMeanData/"+typ+"."+experiment+".zonalmean.nc")
    
    zonal = fz(typ)
    if model == "MMA":
        modzonal = MV.average(zonal,axis=0)
    else:
        models = np.array(eval(zonal.getAxis(0).models))
        i = np.where(models == model)[0][0]
        print i
        modzonal = zonal[i]
    cdutil.setTimeBoundsMonthly(modzonal)
    seasonal_data_unsmooth=getattr(cdutil,season)(modzonal,criteriaarg=(1,None))
    if experiment == "OBS":
        seasonal_data_unsmooth = 100.*seasonal_data_unsmooth
    nt,nlat = seasonal_data_unsmooth.shape
    [lat_plot(seasonal_data_unsmooth[j],color=cm.RdYlBu(j/float(nt))) for j in range(nt)]
    plt.figure()
    seasonal_data = pf.spatially_smooth(seasonal_data_unsmooth,sigma=5)
    
    
    
    #ax1 = plt.subplot(221)
    ax1 = plt.subplot(111)
    [lat_plot(seasonal_data[j],color=cm.RdYlBu(j/float(nt)),alpha=.5) for j in range(nt)]
    

    fz.close()
    plt.figure()
    for L in range(3):

        axn = plt.subplot(2,2,L+2)
        func = funcs[L]
        flat = cdms.open("ExpandingTropicIndices/"+experiment+"/Smooth.5."+func+"_lat."+season+"."+experiment+"."+typ+".nc")
        lats = flat("lat")
        print lats.shape
        fval = cdms.open("ExpandingTropicIndices/"+experiment+"/Smooth.5."+func+"_value."+season+"."+experiment+"."+typ+".nc")
        vals = fval(typ)
        if model == "MMA":
            x = MV.average(lats,axis=0).asma()
            y = MV.average(vals,axis=0).asma()
        else:
            x = lats[i].asma()
            y = vals[i].asma()
        if experiment == "OBS":
            y = y*100.
        t = cmip5.get_plottable_time(lats[0])
        axn.plot(t,x,c="k")
        axn.set_title(func)
        latitude_label_ticks(axn,axis='y')
        for j in range(nt):
            ax1.plot([x[j]],[y[j]],"o",color=cm.RdYlBu(j/float(nt)),mec=cm.RdYlBu(j/float(nt)))
        fval.close()
        flat.close()
    ax1.set_xlim(-50,50)
    #ax1.set_ylim(30,70)
    latitude_label_ticks(ax1,axis='x')
    iax = inset_axes(ax1,width="90%",height="10%",loc=9)
    cmap = cm.RdYlBu
    
    t = get_plottable_time(seasonal_data)
    norm = mpl.colors.Normalize(vmin=min(t),vmax=max(t))
    cb1 = mpl.colorbar.ColorbarBase(iax,cmap=cmap,norm=norm,orientation="horizontal")
    return ax1
    
                        
def get_mma_array(season,typ,xy="lat",experiment="historical-rcp85"):
    if xy == "lat":
        variab = "lat"
    else:
        variab = typ
    funcs =  ["SH_trough","center_max","NH_trough"]
    func = funcs[0]
    f = cdms.open("ExpandingTropicIndices/"+experiment+"/Smooth.5."+func+"_"+xy+"."+season+"."+experiment+"."+typ+".nc")
    shape = f[variab].shape
    f.close()
    S = MV.zeros(shape+(3,))+1.e20
    for i in range(3):
        func = funcs[i]
        f = cdms.open("ExpandingTropicIndices/"+experiment+"/Smooth.5."+func+"_"+xy+"."+season+"."+experiment+"."+typ+".nc")
        data = f(variab)
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
def get_lighter_seasoncolor(season):
    seasoncolors2={}
    seasoncolors2["JJA"]=cm.Reds(.6)
    seasoncolors2["DJF"]=cm.Blues(.6)
    seasoncolors2["MAM"]=cm.Greens(.6)
    seasoncolors2["SON"]=cm.Oranges(.3)
    return seasoncolors2[season]
def plot_time_series(season,typ,anom=True,include_individual=True,experiment="historical-rcp85",xy="lat"):
    
    S = get_mma_array(season,typ,experiment=experiment,xy=xy)
    if anom:
        if experiment is "historical-rcp85":
            S = cmip5.time_anomalies(S,start='1984-1-1',stop='2008-12-31')
        else:
            S = cmip5.time_anomalies(S)
    nmod = S.shape[0]
    mma = MV.average(S,axis=0)[1:-1]
    axes=[]
    for i in range(3):
        axes+=[plt.subplot(1,3,i+1)]
        if include_individual:
            seasoncolors2={}
            seasoncolors2["JJA"]=cm.Reds(.6)
            seasoncolors2["DJF"]=cm.Blues(.6)
            seasoncolors2["MAM"]=cm.Greens(.6)
            seasoncolors2["SON"]=cm.Oranges(.3)
            [time_plot(S[m,:,i],color=seasoncolors2[season],alpha=.3) for m in range(nmod)]
        time_plot(mma[:,i],color=get_seasoncolor(season),lw=3)
    if typ == "pr":
        negsup = "Driest"
        possup="Wettest"
    else:
        negsup="Clearest"
        possup = "Cloudiest"
    
    ax1,ax2,ax3=axes
    ax1.set_title("SH "+negsup+" Latitude")

    ax2.set_title(possup+" Latitude")

    ax3.set_title("NH "+negsup+" Latitude")

        
def get_eof(season,typ,experiment,xy="lat"):
    
    S = get_mma_array(season,typ,experiment=experiment,xy=xy)
    if experiment == "historical-rcp85":
        start = "1981-1-1"
        stop = "2009-12-31"
    else:
        start = cmip5.start_time(S)
        stop = cmip5.stop_time(S)
    
    

    Sa = cmip5.time_anomalies(S,start=start,stop=stop)
    mma = MV.average(Sa,axis=0)[1:-1]
    solver = Eof(mma)
    return solver

def projections(season,typ,xy="lat",project_on = "1pctCO2"):
    S = get_mma_array(season,typ,experiment="OBS")
    S = cmip5.time_anomalies(S)
    P = MV.zeros(S.shape[:-1])
    solver = get_eof(season,typ,project_on,xy=xy)
    fac = cmip5.get_orientation(solver)
    Sf = cmip5.cdms_clone(MV.filled(S),S)
    for i in range(2):
        proj = solver.projectField(Sf[i])[:,0]*fac
        proj = MV.masked_where(np.abs(proj)>1.e10,proj)
        P[i] = proj
    P.setAxisList(S.getAxisList()[:-1])
    return P

def noise(season,typ,xy='lat',anom=True):
    S = get_mma_array(season,typ,experiment="piControl",xy=xy)
    if anom:
        S = cmip5.time_anomalies(S)
    tax = S.getTime()
    latax = S.getAxis(-1)
    modax = S.getAxis(0)
    good = np.where(~S[:,10,0].mask)[0]
    newmodels = np.array(eval(modax.models))[good]
    Snew = MV.array(S.asma()[good])
    Snew.setAxis(1,tax)
    Snew.setAxis(2,latax)
    newmodax = cdms.createAxis(np.arange(Snew.shape[0])*1.)
    newmodax.models=str(newmodels.tolist())
    newmodax.id= "model"
    Snew.setAxis(0,newmodax)
    return Snew

def noise_individual_latitudes(season,typ,xy='lat',axes=None,for_giss=False):
    Snew = noise(season,typ,xy=xy)
    Snew_conc = da.concatenate_this(Snew)
    OBS = get_mma_array(season,typ,experiment="OBS",xy=xy)
    OBS = cmip5.time_anomalies(OBS)
    start = cmip5.start_time(OBS)
    stop = cmip5.stop_time(OBS)
    H = get_mma_array(season,typ,experiment="historical-rcp85_ensemble",xy=xy)(time=(start,stop))
    H = cmip5.time_anomalies(H)
    mod_trends = genutil.statistics.linearregression(H,axis=1,nointercept=1)*3650
    if for_giss:
        giss= np.where(np.array([x.find("GISS")>=0 for x in eval(H.getAxis(0).models)]))[0]
        mod_trends_giss = mod_trends.asma()[giss]
    
    
    obs_trends = genutil.statistics.linearregression(OBS,axis=1,nointercept=1)*120 #obs time axis in months; multiply by 120 to get trends per decade
    if axes is None:
        ownaxes = False
        axes = []
    else:
        ownaxes = True
    for i in range(3):
        if not ownaxes:
            ax = plt.subplot(1,3,i+1)
            axes +=[ax]
        else:
            ax = axes[i]
        test =da.get_slopes(Snew_conc[:,i],26)
        noisenorm = np.std(test)
        #ax.hist(test/noisenorm,normed=True,color=get_seasoncolor(season),ec="w")
        #ax.hist(mod_trends[:,i]/noisenorm,normed=True,color="w",ec=get_seasoncolor(season))
        da.fit_normals_to_data(test/noisenorm,color=get_seasoncolor(season),lw=4,label="piControl",ax=ax)
        da.fit_normals_to_data(mod_trends[:,i]/noisenorm,a=4,color=get_lighter_seasoncolor(season),lw=4,label="Hist-RCP85",ax=ax)
        
        #ax.hist(mod_trends[:,i],normed=True,color = get_lighter_seasoncolor(season),ec="w")
        isccp,patmos = obs_trends[:,i]
        ax.axvline(isccp/noisenorm,lw=3,ls="-",color=cm.Greys(.9),label="ISCCP")
        if for_giss:
                for gissthing in mod_trends_giss[:,i]:
                        ax.axvline(gissthing/noisenorm,lw=3,ls=":",color=cm.Greys(.4))
                ax.axvline(mod_trends_giss[:,i][0]/noisenorm,lw=3,ls=":",color=cm.Greys(.4),label="GISS models")
        #ax.axvline(patmos/noisenorm,lw=3,ls="--",color=cm.Greys(.9),label="PATMOSx")
        
        maxsignal=max([np.abs(isccp/noisenorm),np.abs(patmos/noisenorm)])
        if maxsignal >= stats.norm.interval(.9)[1]:
            ax.set_axis_bgcolor('beige')
        ax.set_xlim(-4,4)
        if i == 0:
            ax.legend(ncol=4,fontsize=6)
        
        
def noise_individual_latitudes_ALL_SEASONS(typ,xy='lat',axes=None,for_giss=False):
    seasons = ["DJF","MAM","JJA","SON"]
    for s in range(len(seasons)):
        season = seasons[s]
        axes = [plt.subplot(4,3,i+1+3*s) for i in range(3)]
        noise_individual_latitudes(season,typ,xy=xy,axes=axes,for_giss=for_giss)
    axes = plt.gcf().axes
    [ax.set_xlabel("S/N") for ax in axes[-3:]]
    [ax.legend(loc=0,fontsize=6) for ax in axes]
    axes[0].set_title("SH trough")
    axes[2].set_title("NH trough")
    axes[1].set_title("Central Peak")
    axes[0].set_ylabel("DJF")
    axes[3].set_ylabel("MAM")
    axes[6].set_ylabel("JJA")
    axes[9].set_ylabel("SON")
        
        
def signal_to_noise_histograms(season,typ,xy='lat',project_on="1pctCO2" ):
    proj = noise_projection(season,typ,xy=xy,project_on=project_on)
    test =da.get_slopes(proj,26)
    noise = np.std(test)
    
    obs = projections(season,typ,xy=xy,project_on=project_on)
    plt.hist(test/noise,normed=True,color=get_seasoncolor(season),ec="w")
    isccp,patmos = genutil.statistics.linearregression(obs,axis=1,nointercept=1)*120
    plt.axvline(isccp/noise,lw=3,ls="-",color=cm.Greys(.9),label="ISCCP")
    plt.axvline(patmos/noise,lw=3,ls="--",color=cm.Greys(.9),label="PATMOSx")
    
def all_seasons_signal_to_noise_histograms(typ,xy='lat',project_on="1pctCO2" ):
    seasons = ["DJF","MAM","JJA","SON"]
    for i in range(4):
        season= seasons[i]
        plt.subplot(2,2,i+1)
        signal_to_noise_histograms(season,typ,xy=xy,project_on=project_on)
        plt.title(season)
        
def noise_projection(season,typ,xy='lat',project_on="1pctCO2"):
    
    Snew = noise(season,typ,xy=xy)
    solver = get_eof(season,typ,project_on,xy=xy)
    Snewf = cmip5.cdms_clone(Snew.filled(),Snew)
    proj = solver.projectField(da.concatenate_this(Snewf))[:,0]
    proj = MV.masked_where(np.abs(proj)>1.e10,proj)
    return proj
    

  
    #solver = get_eof(season,typ,project_on,xy=xy)

def plot_seasonal_projections(typ,xy="lat",project_on = "1pctCO2"):
    seasons = ["DJF","MAM","JJA","SON"]
    for i in range(4):
        season= seasons[i]
        plt.subplot(2,2,i+1)
        P = projections(season,typ,xy=xy,project_on = project_on)
        datasets = eval(P.getAxis(0).models)
        time_plot(P[0],label=datasets[0],color=get_seasoncolor(season),lw=3)
        x = cmip5.get_plottable_time(P[0])
        y = P[0].asma()
        p = np.ma.polyfit(x,y,1)
        plt.plot(x,np.polyval(p,x),"-",color=get_seasoncolor(season),lw=3)
        time_plot(P[1],label=datasets[1],color=get_seasoncolor(season),lw=3,ls="--")

        x = cmip5.get_plottable_time(P[1])
        y = P[1].asma()
        p = np.ma.polyfit(x,y,1)
        plt.plot(x,np.polyval(p,x),"--",color=get_seasoncolor(season),lw=3)
        plt.legend(loc=0)
        plt.title(season)
        

def plot_eof(season,typ,**kwargs):
    seasons = ["DJF","MAM","JJA","SON"]
    experiment = kwargs.pop("experiment","historical-rcp85")
    S = get_mma_array(season,typ,experiment=experiment)
    if experiment == "historical-rcp85":
        start = "1981-1-1"
        stop = "2009-12-31"
    else:
        start = cmip5.start_time(S)
        stop = cmip5.stop_time(S)
    Sclim = MV.average(MV.average(S(time=(start,stop)),axis=1),axis=0)
    

    Sa = cmip5.time_anomalies(S,start=start,stop=stop)
    mma = MV.average(Sa,axis=0)[1:-1]
    solver = Eof(mma)
    ax1 = plt.subplot(211)
    fac = cmip5.get_orientation(solver)
    eof1 = fac*solver.eofs()[0]
    plt.xlim(-.01,.07)
    plt.ylim(-60,60)
    
    #plt.plot(Sclim.asma(),eof1.asma(),"o-",**kwargs)
    for i in range(3):
        x = seasons.index(season)*.02
        y = float(Sclim[i])
        print "y = "+str(y)
        dy = eof1[i]*20
        print dy
        dx = 0
        width=0.01
        ax1.add_patch(patches.Arrow(x,y,dx,dy,width=width,color=get_seasoncolor(season)))
        
    plt.axhline(0,c="k",ls=":")

    plt.subplot(212)
    time_plot(fac*solver.pcs()[:,0],color=get_seasoncolor(season))

    print solver.varianceFraction()[0]
    return solver
def barplot_trends(typ,start='1900-1-1',stop='2100-1-1',experiment="historical-rcp85"):
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
        
        S = get_mma_array(season,typ,experiment=experiment)
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
        
    
        





def write_obs():
    T2 = TropicTracker()
    rawstart = cmip5.start_time(T2.ISCCP_raw)
    rawstop = cmip5.stop_time(T2.ISCCP_raw)
    
    CLT = MV.zeros(((2,)+T2.ISCCP.shape))
    CLT_RAW = MV.zeros((2,)+T2.ISCCP_raw.shape)
    PR = MV.zeros(((2,)+T2.GPCP.shape))

    CLT[0] = T2.ISCCP
    CLT[1]=T2.PATMOSX
    cltax = cdms.createAxis(np.arange(2))
    cltax.models = str(["ISCCP","PATMOSX"])
    CLT.setAxis(0,cltax)
    CLT.setAxis(1,T2.ISCCP.getTime())
    CLT.setAxis(2,T2.ISCCP.getLatitude())
    CLT.id = "clt"
    f = cdms.open("ZonalMeanData/clt.OBS.zonalmean.nc","w")
    f.write(CLT)
    f.close()

    
    CLT_RAW[0] = T2.ISCCP_raw(time=(rawstart,rawstop))
    CLT_RAW[1]=T2.PATMOSX(time=(rawstart,rawstop))
    clt_rawax = cdms.createAxis(np.arange(2))
    clt_rawax.models = str(["ISCCP_raw","PATMOSX_raw"])
    CLT_RAW.setAxis(0,clt_rawax)
    CLT_RAW.setAxis(1,T2.ISCCP_raw.getTime())
    CLT_RAW.setAxis(2,T2.ISCCP_raw.getLatitude())
    CLT_RAW.id = "clt"
    f = cdms.open("ZonalMeanData/clt.OBS_RAW.zonalmean.nc","w")
    f.write(CLT_RAW)
    f.close()

    PR[0] = T2.GPCP
    PR[1]=T2.CMAP
    prax = cdms.createAxis(np.arange(2))
    prax.models = str(["ISCCP","PATMOSX"])
    PR.setAxis(0,prax)
    PR.setAxis(1,T2.GPCP.getTime())
    PR.setAxis(2,T2.GPCP.getLatitude())
    PR.id = "pr"
    f = cdms.open("ZonalMeanData/pr.OBS.zonalmean.nc","w")
    f.write(PR)
    f.close()


def make_composites(season,typ,xy='lat',flavor='la nina'):
    Snew = noise(season,typ,xy=xy,anom=False)
    fssts = cdms.open("piC_SSTS.nc")
    sst=fssts("sst")
    cdutil.setTimeBoundsMonthly(sst)
    djf = cdutil.DJF(sst,criteriaarg=(1,None))
    if season != "DJF":
        djf = djf[:,:-1]
    djfmeans = np.ma.average(djf.asma(),axis=1)
    djfsig = np.ma.std(djf.asma(),axis=1)
    if flavor == 'la nina':
        boundary = (djfmeans-2*djfsig)
        criterion_for_rejection = djf.asma()>boundary[:,np.newaxis]
    else:
        boundary = (djfmeans+2*djfsig)
        criterion_for_rejection = djf.asma()<boundary[:,np.newaxis]
        
    ensomask = np.repeat(criterion_for_rejection[:,:,np.newaxis],3,axis=-1)
    return MV.average(MV.masked_where(ensomask,Snew),axis=1)

def la_nina_fingerprints_all_four_seasons(typ,xy='lat',flavor='la nina'):
    seasons = ["DJF","MAM","JJA","SON"]
    ax = plt.subplot(111)
    x=0
    for season in seasons:
        Snew = noise(season,typ,xy=xy,anom=False)
        clim = MV.average(Snew,axis=1)
        clim_mma = MV.average(clim,axis=0)
        comp = make_composites(season,typ,xy=xy,flavor=flavor)
        mma = MV.average(comp-clim,axis=0)
        for j in range(3):
            y = clim_mma[j]
            dy = mma[j]*10
            print y
            print dy
            dx=0
            #plt.arrow(x, y, 0, dy)
            width=0.2
            ax.add_patch(patches.Arrow(x,y,dx,dy,width=width,color=get_seasoncolor(season)))
        x+=1
    plt.xlim(-1,5)
    plt.ylim(-30,30)
        
def get_composites():
    f = cdms.open("/Users/kmarvel/Google Drive/SEASONAL/ENSO_indices.nc")
    ENSO = f("sst")
    models = ['ACCESS1-0',\
        'ACCESS1-3',\
        'BNU-ESM',\
        'CCSM4',\
        'CESM1-CAM5',\
        'CSIRO-Mk3-6-0',\
        'CanESM2',\
        'FGOALS-s2',\
        'GISS-E2-H',\
        'GISS-E2-H-CC',\
        'GISS-E2-R-CC',\
        'IPSL-CM5A-MR',\
        'IPSL-CM5B-LR',\
        'MIROC-ESM',\
        'MIROC-ESM-CHEM',\
        'MIROC5',\
        'MPI-ESM-LR',\
        'MPI-ESM-MR',\
        'NorESM1-M',\
        'NorESM1-ME']
    seasons =["DJF","MAM","JJA","SON"]
    EL_NINO = MV.zeros((len(models),72,4))
    LA_NINA = MV.zeros((len(models),72,4))
    CLIMATOLOGY = MV.zeros((len(models),72,4))
    fz = cdms.open("ZonalMeanData/clt.piControl.zonalmean.nc")
    zonal = fz("clt")
    piCfnames = eval(zonal.getAxis(0).models)
    
    for model in models:
        modi = models.index(model)
        print model
        piCfnamei = int(np.where(np.array([x.find(model+".piControl.r1i1p1")>=0 for x in piCfnames]))[0])
        modzonal = zonal[piCfnamei]
        
        an = ENSO[modi].anom()
        sig = np.ma.std(an)
        I_en = np.where(an>2.*sig)[0]
        I_ln = np.where(an<-2.*sig)[0]

        for season in seasons:
            seasoni = seasons.index(season)
            seasonal_data = getattr(cdutil,season)(modzonal,criteriaarg=(1,None))
            
            LN = cmip5.cdms_clone(MV.average(seasonal_data.asma()[I_ln],axis =0),seasonal_data[1])
            EN = cmip5.cdms_clone(MV.average(seasonal_data.asma()[I_en],axis =0),seasonal_data[1])
            
            EL_NINO[modi,:,seasoni] = EN
            LA_NINA[modi,:,seasoni] = LN
            CLIMATOLOGY[modi,:,seasoni]= MV.average(seasonal_data,axis=0)
    modax = cdms.createAxis(range(len(models)))
    modax.models = str(models)
    latax = zonal.getLatitude()

    seasonax = cdms.createAxis(range(len(seasons)))
    seasonax.seasons = str(seasons)
    EL_NINO.setAxis(0,modax)
    EL_NINO.setAxis(1,latax)
    EL_NINO.setAxis(2,seasonax)
    EL_NINO.id = "clt"
    
    LA_NINA.setAxis(0,modax)
    LA_NINA.setAxis(1,latax)
    LA_NINA.setAxis(2,seasonax)
    LA_NINA.id = "clt"

    CLIMATOLOGY.setAxis(0,modax)
    CLIMATOLOGY.setAxis(1,latax)
    CLIMATOLOGY.setAxis(2,seasonax)
    CLIMATOLOGY.id = "clt"
    
    fz.close()
    f.close()
    return CLIMATOLOGY,EL_NINO,LA_NINA
    
    
def plot_composites(model,season,C,E,L):
    models = eval(C.getAxis(0).models)
    seasons = eval(C.getAxis(2).seasons)
    modi = models.index(model)
    seasoni = seasons.index(season)
    lat_plot(C[modi,:,seasoni],c="k",label="Climatology")
    lat_plot(E[modi,:,seasoni],c="r",label="El Nino")
    lat_plot(L[modi,:,seasoni],c="b",label="La Nina")
    plt.legend(loc=0)

def plot_shifts(C,X):
      seasons = eval(C.getAxis(2).seasons)
      ENSO = MV.zeros((20,3,4))
      CLIM = MV.zeros((20,3,4))
      for i in range(20):
            new = SH_trough(MV.transpose(X[i]))
            old = SH_trough(MV.transpose(C[i]))
            ENSO[i,0]=new
            CLIM[i,0]=old

            new = center_max(MV.transpose(X[i]))
            old = center_max(MV.transpose(C[i]))
            
            ENSO[i,1]=new
            CLIM[i,1]=old

            new = NH_trough(MV.transpose(X[i]))
            old = NH_trough(MV.transpose(C[i]))
            ENSO[i,2]=new
            CLIM[i,2]=old
      return CLIM,ENSO
            
      
    
