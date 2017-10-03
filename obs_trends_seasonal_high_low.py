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
import seasonal_cycle_utils as sc
from Plotting import *
import CMIP5_tools as cmip5
import DA_tools 

class HeightObs():
    def __init__(self):
        fw = cdms.open("/Users/kmarvel/Google Drive/CLOUD_SEASONS/cloud-seasons/HEIGHT/masked_low_high.nc","r")
    
        self.isccp_low = fw("isccp_low")
        self.isccp_high = fw("isccp_high")
        self.isccp_low.description = "ISCCP low + mid-level cloud"
        self.isccp_high.description = "ISCCP high cloud"
        
        self.patmos_low = fw("patmos_low")
        self.patmos_high = fw("patmos_high")
        self.patmos_low.description = "PATMOS low + mid-level cloud"
        self.patmos_high.description = "PATMOS high cloud"
        fw.close()

        
    def plot_trends(self,season="ANNUALCYCLE",significance=None):
        func = getattr(cdutil,season).departures
        data = [self.isccp_low,self.isccp_high,self.patmos_low,self.patmos_high]
        for i in range(4):
            plt.subplot(2,2,i+1)
            X = func(data[i],criteriaarg=(.99,None))
            trends = cmip5.get_linear_trends(X,significance=significance)
            m = bmap(trends,cmap="RdBu",vmin=-4,vmax=4)
            plt.colorbar(orientation="horizontal",label="Trend (%/decade)")
            plt.title(getattr(data[i],"description"))
            m.drawcoastlines()
    def plot_zonal_trends(self,season="ANNUALCYCLE",significance=None):
        func = getattr(cdutil,season).departures
        lowdata = [self.isccp_low,self.patmos_low]
        highdata = [self.isccp_high,self.patmos_high]
        plt.subplot(212)
        for i in range(2):
            X = func(cdutil.averager(lowdata[i],axis='x'),criteriaarg=(.99,None))
            trends = cmip5.get_linear_trends(X,significance=significance)
            lat_plot(trends,label=getattr(lowdata[i],"description"),lw=3)
            plt.axhline(0,c="k",ls=":")
            plt.legend(loc=0)
        plt.title("LOW")
        plt.subplot(211)
        for i in range(2):
            X = func(cdutil.averager(highdata[i],axis='x'),criteriaarg=(.99,None))
            trends = cmip5.get_linear_trends(X,significance=significance)
            lat_plot(trends,label=getattr(highdata[i],"description"),lw=3)
            plt.axhline(0,c="k",ls=":")
            plt.legend(loc=0)
        plt.title("HIGH")
    def project_fields(self,season="YEAR"):
        fingerprints = solver("1pctCO2",season)
        projections = {}
        for dataset in ["isccp","patmos"]:
            projections[dataset]={}
            for typ in ["low","high"]:
                print dataset+" "+typ
                data = getattr(self,dataset+"_"+typ)
                zdata = cdutil.averager(data,axis='x')
                cdutil.setTimeBoundsMonthly(zdata)
                anom = getattr(cdutil,season).departures(zdata)
                anom = cmip5.cdms_clone(anom.filled(fill_value=1.e20),anom)
                
                fp = fingerprints[typ]
                fac = DA_tools.get_orientation(fp)
               
                proj = fac*fp.projectField(anom)[:,0]
                proj = MV.masked_where(np.abs(proj)>1.e10,proj)
                projections[dataset][typ]=proj
               
        return projections
                

            
            
        
    
            
        
def get_zonal_data(experiment,season="YEAR",anom=True):
    fhigh = cdms.open(glob.glob("ZONAL_HEIGHT/HIGH/*"+experiment+"*")[0])
    high = fhigh("high_clt",latitude=(-50,50))
    fhigh.close()
    fmid = cdms.open(glob.glob("ZONAL_HEIGHT/MID/*"+experiment+"*")[0])
    mid = fmid("mid_clt",latitude=(-50,50))
    fmid.close()
    flow = cdms.open(glob.glob("ZONAL_HEIGHT/LOW/*"+experiment+"*")[0])
    low = flow("low_clt",latitude=(-50,50))
    flow.close()
    if experiment == "historical":
        if season != "DJF":
            high = high(time=('1981-1-1','2004-12-31'))
            mid = mid(time=('1981-1-1','2004-12-31'))
            low = low(time=('1981-1-1','2004-12-31')) 
    low = low+mid
    cdutil.setTimeBoundsMonthly(low)
    cdutil.setTimeBoundsMonthly(high)
    if anom:
        lowanom = getattr(cdutil,season).departures(low,criteriaarg=(.99,None))
        highanom = getattr(cdutil,season).departures(high,criteriaarg=(.99,None))
    else:
        lowanom = getattr(cdutil,season)(low,criteriaarg=(.99,None))
        highanom = getattr(cdutil,season)(high,criteriaarg=(.99,None))
        
    if experiment == "piControl":
        if season !="DJF":
            lowanom = DA_tools.concatenate_this(lowanom)
            highanom = DA_tools.concatenate_this(highanom)
        else:
            lowanom = DA_tools.concatenate_this(lowanom(time=('0001-12-1','0020-11-31')))
            highanom = DA_tools.concatenate_this(highanom(time=('0001-12-1','0020-11-31')))       
    lowanom.id= "low_clt"
    highanom.id="high_clt"
    return lowanom,highanom    
#class solver(experiment):
def solver(experiment,season="YEAR"):
    fhigh = cdms.open(glob.glob("ZONAL_HEIGHT/HIGH/*"+experiment+"*")[0])
    high = fhigh("high_clt",latitude=(-50,50))
    fhigh.close()
    fmid = cdms.open(glob.glob("ZONAL_HEIGHT/MID/*"+experiment+"*")[0])
    mid = fmid("mid_clt",latitude=(-50,50))
    fmid.close()
    flow = cdms.open(glob.glob("ZONAL_HEIGHT/LOW/*"+experiment+"*")[0])
    low = flow("low_clt",latitude=(-50,50))
    flow.close()
    if experiment == "historical":
        if season != "DJF":
            high = high(time=('1981-1-1','2004-12-31'))
            mid = mid(time=('1981-1-1','2004-12-31'))
            low = low(time=('1981-1-1','2004-12-31')) 
    low = low+mid
    cdutil.setTimeBoundsMonthly(low)
    cdutil.setTimeBoundsMonthly(high)
    total_cloud = low+high
    
    lowanom = getattr(cdutil,season).departures(low,criteriaarg=(.99,None))
    highanom = getattr(cdutil,season).departures(high,criteriaarg=(.99,None))
    total_anom = getattr(cdutil,season).departures(total_cloud,criteriaarg=(.99,None))
    if experiment == "1pctCO2":
        if season != "DJF":
            lowmma = MV.average(lowanom,axis=0)
            highmma = MV.average(highanom,axis=0)
            total_mma = MV.average(total_anom,axis=0)
        else:
            
            lowmma = MV.average(lowanom,axis=0)[11:-2]
            highmma = MV.average(highanom,axis=0)[11:-2]
            total_mma = MV.average(total_anom,axis=0)[11:-2]
            print cmip5.stop_time(lowmma)
    elif experiment == "piControl":
        if season !="DJF":
            lowmma = DA_tools.concatenate_this(lowanom)
            highmma = DA_tools.concatenate_this(highanom)
            total_mma = DA_tools.concatenate_this(total_anom)
        else:
          
            lowmma = DA_tools.concatenate_this(lowanom(time=('0001-12-1','0020-11-31')))
            highmma = DA_tools.concatenate_this(highanom(time=('0001-12-1','0020-11-31')))
            total_mma = DA_tools.concatenate_this(total_anom(time=('0001-12-1','0020-11-31')))
    multisolver = MultivariateEof([lowmma,highmma])
    lowsolver= Eof(lowmma)
    highsolver = Eof(highmma)
    total_solver = Eof(total_mma)
    D={}
    D["high"]=highsolver
    D["low"]=lowsolver
    D["multi"]=multisolver
    D["total"]=total_solver
    return D

def plot_eofs(D):
    loweofs,higheofs = D["multi"].eofs()
    plt.subplot(212)
    lat_plot(DA_tools.get_orientation(D["low"])*D["low"].eofs()[0],label="Single-variable EOF")
    lat_plot(DA_tools.get_orientation(D["multi"])*loweofs[0],label="Multivariate EOF")
    plt.title("LOW")
    plt.legend(loc=0)
    plt.axhline(0,c="k",ls=":")
    
    plt.subplot(211)
    lat_plot(DA_tools.get_orientation(D["high"])*D["high"].eofs()[0],label="Single-variable EOF")
    lat_plot(DA_tools.get_orientation(D["multi"])*higheofs[0],label="Multivariate EOF")
    plt.title("HIGH")
    plt.legend(loc=0)
    plt.axhline(0,c="k",ls=":")

def plot_pcs(D):
    time_plot(DA_tools.get_orientation(D["high"])*D["high"].pcs()[:,0],label="Single-variable High")
    time_plot(DA_tools.get_orientation(D["low"])*D["low"].pcs()[:,0],label="Single-variable Low")

    time_plot(DA_tools.get_orientation(D["multi"])*D["multi"].pcs()[:,0],label="Multivariate")
    
    
    
def get_overlapping_slopes(X,L):
    trends = np.ma.array([])
    endi=len(X)-L
    i=0
    while i<=endi:
        trunc = X[i:i+L]
        trends=np.ma.append(trends,float(cmip5.get_linear_trends(trunc)))
        i+=1
    return trends

def is_detectable(season,plot=True,SN=False):

    #Get the solvers from the MMA (historical) 
    historical = solver("1pctCO2",season=season)
    #piC = solver("piControl",season=season)

    #Get concatenated piControl seasonal anomalies
    lowanom_piC,highanom_piC = get_zonal_data("piControl",season=season)

    #project piControl onto solvers
    lowproj_piC=historical["low"].projectField(lowanom_piC)[:,0]
    highproj_piC=historical["high"].projectField(highanom_piC)[:,0]
    #multiproj_piC=historical["multi"].projectField([lowanom_piC,highanom_piC])[:,0]
    
    lowanom_historical,highanom_historical = get_zonal_data("historical")
    nmod,nyears,nlat = lowanom_historical.shape

    lowproj_historical = MV.zeros(nmod)
    highproj_historical = MV.zeros(nmod)
   # multiproj_historical = MV.zeros(nmod)


    H = HeightObs()
    p = H.project_fields(season)
    for i in range(nmod):
        toproj_low = lowanom_historical[i]
        toproj_high = highanom_historical[i]
        #toproj_multi=[toproj_low,toproj_high]

        proj_low = historical["low"].projectField(toproj_low)[:,0]
        proj_high = historical["high"].projectField(toproj_high)[:,0]
        #proj_multi = historical["multi"].projectField(toproj_multi)[:,0]

        lowproj_historical[i] = float(cmip5.get_linear_trends(proj_low))
        highproj_historical[i] = float(cmip5.get_linear_trends(proj_high))
        #multiproj_historical[i] = float(cmip5.get_linear_trends(proj_multi))

    lownoise = get_overlapping_slopes(lowproj_piC,nyears)
    highnoise = get_overlapping_slopes(highproj_piC,nyears)
    #multinoise = get_overlapping_slopes(multiproj_piC,nmod)


    if plot:
        plt.subplot(121)
        if SN:
            fac = np.ma.std(lownoise)
        else:
            fac=1.
        plt.hist(lownoise/fac,alpha=.5,normed=True,color=cm.Greens(.7),ec=cm.Greens(.7))
        DA_tools.fit_normals_to_data(lownoise/fac,c=cm.Greens(.7),lw=3,label="piControl")
        plt.hist(lowproj_historical/fac,alpha=.5,normed=True,color=cm.Oranges(.7),ec=cm.Oranges(.7))
        DA_tools.fit_normals_to_data(lowproj_historical/fac,c=cm.Oranges(.7),lw=3,label="historical")
        isccp = cmip5.get_linear_trends(p["isccp"]["low"])
        patmos = cmip5.get_linear_trends(p["patmos"]["low"])
        plt.axvline(isccp/fac,color=cm.Reds(.8),label="ISCCP",lw=3)
        plt.axvline(patmos/fac,color=cm.Blues(.8),label="PATMOS-x",lw=3)
        
        
        plt.title("LOW")


        plt.subplot(122)
        if SN:
            fac = np.ma.std(highnoise)
        else:
            fac=1.
        plt.hist(highnoise/fac,alpha=.5,normed=True,color=cm.Greens(.7),ec=cm.Greens(.7))
        DA_tools.fit_normals_to_data(highnoise/fac,c=cm.Greens(.7),lw=3,label="piControl")
        plt.hist(highproj_historical/fac,alpha=.5,normed=True,color=cm.Oranges(.7),ec=cm.Oranges(.7))
        DA_tools.fit_normals_to_data(highproj_historical/fac,c=cm.Oranges(.7),lw=3,label="historical")
        isccp = cmip5.get_linear_trends(p["isccp"]["high"])
        patmos = cmip5.get_linear_trends(p["patmos"]["high"])
        plt.axvline(isccp/fac,color=cm.Reds(.8),label="ISCCP",lw=3)
        plt.axvline(patmos/fac,color=cm.Blues(.8),label="PATMOS-x",lw=3)
        plt.title("HIGH")

        # plt.subplot(133)

        # plt.hist(multinoise,alpha=.5,normed=True,color=cm.Greens(.7),ec=cm.Greens(.7))
        # DA_tools.fit_normals_to_data(multinoise,c=cm.Greens(.7),lw=3)
        # plt.hist(multiproj_historical,alpha=.5,normed=True,color=cm.Oranges(.7),ec=cm.Oranges(.7))
        # DA_tools.fit_normals_to_data(multiproj_historical,c=cm.Oranges(.7),lw=3)
        
        # plt.title("Multivariate")
    low_SN = lowproj_historical/np.ma.std(lownoise)
    high_SN = highproj_historical/np.ma.std(highnoise)
#    multi_SN = multiproj_historical/np.ma.std(multinoise)
    
    return low_SN,high_SN#,multi_SN,

def all_seasons(conf=.95):
    threshold=stats.norm.interval(conf)[1]
    print "SEASON \t LOW \t HIGH \t"
    for season in ["DJF","MAM","JJA","SON","YEAR"]:
        low,high=is_detectable(season,plot=False)
        low_detectable=float(len(np.where(low>threshold)[0]))/float(len(low))
        high_detectable = float(len(np.where(high>threshold)[0]))/float(len(high))
        print season+" \t "+str(low_detectable)+" \t "+str(high_detectable)
    
    
        
    
        
