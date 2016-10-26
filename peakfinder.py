### Import useful routines
import numpy as np
import string
import glob
import os,sys
from collections import Counter
from string import upper,lower,join
import itertools

### Import CDAT routines ###
import MV2 as MV
import cdms2 as cdms
import genutil
import cdutil
import cdtime
from eofs.cdms import Eof

### Import scipy routines for smoothing, interpolation
from scipy.interpolate import interp1d
from scipy.optimize import brentq,fminbound
import scipy.ndimage as ndimag
import scipy.stats as stats

### Import plotting routines
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap  
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

from matplotlib.patches import Ellipse,Rectangle
import matplotlib.mlab as mlab

### Import stuff I wrote

#from Helper import get_orientation,get_plottable_time,get_slopes


### Set classic Netcdf (ver 3)
cdms.setNetcdfShuffleFlag(0)
cdms.setNetcdfDeflateFlag(0)
cdms.setNetcdfDeflateLevelFlag(0)

import pickle


def spatially_smooth(data,sigma=5):
    """Spatially smooth even when spatial data is missing """
    # Convert resolution to grid cells from latitude
    resolution = np.abs(np.median(np.diff(data.getLatitude()[:])))
    
    nt = data.shape[0]
    sigma = sigma/resolution
    D = MV.zeros(data.shape)+1.e20
    
    #If there is no missing data, we can just smooth
    if not hasattr(data,"mask"):
        data.mask=False
    if data.mask.shape!= data.shape:        
        for i in range(nt):
            D[i] = ndimag.gaussian_filter(data[i],sigma=sigma)
    #Otherwise, apply smoothing filter to compressed data
    else:
        for i in range(nt):
            datai = data[i]

            if len(datai.compressed())!=0: #if not all the data is missing at that time step
                datai = MV.array(ndimag.gaussian_filter(datai.compressed(),sigma))
                datai = MV.masked_where(datai>1.e10,datai)
                J = np.where(~data[i].mask)[0]
                
                D[i,J] = datai
        # For data where all lats are masked (ie, missing time step)
        D = MV.masked_where(D>1.e10,D)

    D.setAxis(0,data.getAxis(0))
    D.setAxis(1,data.getLatitude())
    return D


#########################################
#                                       #
#            2D indices from pr        #  
#                                       #
#########################################

def prep_data(data, season=None, sigma = 5.,criteriaarg=[0.99,None]):
    """ Zonally average and smooth seasonal data """
    try:
        units = data.units
    except:
        units = None
    if len(np.where(np.isnan(data).flatten())[0])!=0:
        data = MV.masked_where(np.isnan(data),data)
    cdutil.setTimeBoundsMonthly(data)
    #If it's not already zonally averaged, do it
    if "lon" in data.getAxisIds():
        data = cdutil.averager(data,axis='x')
    #Get seasonal averages.  Require all 3 months to be present.
    if season is not None:
        func = getattr(cdutil,season)
        data = func(data,criteriaarg=criteriaarg)
    
        
   
    # Smooth by sigma
   
    SMOOTH= spatially_smooth(data,sigma=sigma)
 
    if units is not None:
        SMOOTH.units = units
    return SMOOTH



def find_all_peaks(a,return_maxmin=False):
    """Given smoothed zonal average, find peaks and troughs """
    lat = a.getLatitude()[:]
    # If there are missing data, get rid of mask
   # if hasattr(a,"mask"):
    #    mask = a.mask
     #   a = a.compressed()

      #  lat = np.ma.masked_where(mask,lat).compressed()
        
    # Use scipy interpolate to find off-grid maxima/minima
    interpl = interp1d(lat,np.array(a),"cubic")
    neginterpl = interp1d(lat,-1*np.array(a),"cubic")

    #Find local minima and maxima
    maxi= np.where((np.array(a)==ndimag.maximum_filter1d(a,4)))[0]
    mini = np.where((-1*np.array(a)==ndimag.maximum_filter1d(-1*a,4)))[0]
    
    mini = np.ma.masked_where(lat[mini]<=np.min(lat),mini).compressed()
    mini = np.ma.masked_where(lat[mini]>=np.max(lat),mini).compressed()
    maxi = np.ma.masked_where(lat[maxi]<=np.min(lat),maxi).compressed()
    maxi = np.ma.masked_where(lat[maxi]>=np.max(lat),maxi).compressed()
    
    #Get interpolated maxes and minima
    bnd =0.5*np.median(np.diff(lat))
    maxes = [fminbound(neginterpl, lat[x]-bnd,lat[x]+bnd) for x in maxi]
    mins = [fminbound(interpl, lat[x]-bnd,lat[x]+bnd) for x in mini]

    #return_maxmin: separate maxes from mins
    if return_maxmin:
        maxes = np.array(sorted(maxes))
        mins = np.array(sorted(mins))
        return maxes,mins,interpl(maxes),interpl(mins)
    else:
        extrema =  np.array(sorted(np.append(maxes,mins)))
        return extrema, interpl(extrema)


def get_major_features(data, \
              midpoints = False,\
              SHDryBounds = [-45,-15],\
              NHDryBounds = [15,45], \
              SHStormBounds = [-70,-40], \
              NHStormBounds = [40,70], \
              TropicBounds = [-15,15]):

    """ Peak detection algorithm.  Specify latitude bounds for storm tracks, dry zones, and tropics, or use defaults"""

    #If there's no data for this time step (as in obs), return masked array with everything masked.
    if hasattr(data,"mask"):
        if len(np.where(data.mask==True)[0]) == len(data):
            if midpoints:
                checked_lats = MV.ones(11)+1.e20
            else:
                checked_lats = MV.ones(5)+1.e20
            checked_vals = MV.ones(5)+1.e20
            checked_lats = np.ma.masked_where(checked_lats>1.e10,checked_lats)
            checked_vals = np.ma.masked_where(checked_vals>1.e10,checked_vals)
            return checked_lats,checked_vals
    
    #Otherwise, calculate peaks and troughs
    lat = data.getLatitude()[:]
    #Arrays to hold the peak/trough lats and values
    checked_lats = np.zeros(5)
    checked_vals = np.zeros(5)

    #Find the minima
    mx,mn,imaxes,imns = find_all_peaks(data,return_maxmin=True)

    #Do they correspond to dry zones?
    SHDry = []
    NHDry = []
    ###for i in range(len(mn)):
    for minlat in mn:
        
        ###minlat = lat[mn[i]]
        
        if SHDryBounds[0] <= minlat <= SHDryBounds[1]:
            ###SHDry+=[mn[i]]
            SHDry += [minlat]
        if NHDryBounds[0] <= minlat <= NHDryBounds[1]:
            ###NHDry+=[mn[i]]
            NHDry += [minlat]

    #If exactly one latitude is a trough between -45 and -15, it's the SH dry zone
    if len(SHDry)== 1:
        SHDry = SHDry[0]
        checked_lats[1] = SHDry####lat[SHDry]
        checked_vals[1] = imns[np.argwhere(mn == SHDry)]####data[SHDry]
    else:
        checked_lats[1] = 1.e20
        checked_vals[1] = 1.e20
    
       #If exactly one latitude is a trough between 15 and 45, it's the NH dry zone
    if len(NHDry)== 1:
        NHDry = NHDry[0]
        checked_lats[3] = NHDry###lat[NHDry]
        checked_vals[3] =imns[np.argwhere(mn == NHDry)] ###data[NHDry]
    else:
        checked_lats[3] = 1.e20
        checked_vals[3] = 1.e20

     #Find the maxima
    ###mx=detect_local_minima(-data)[0]
    
    #Do they correspond to storm tracks or the tropics?
   
    Tropics = []
    SHStorm = []
    NHStorm = []
    ###for i in range(len(mx)):
    for maxlat in mx:
       ### maxlat = lat[mx[i]]

        if SHStormBounds[0] <= maxlat <= SHStormBounds[1]:
            ###SHStorm+=[mx[i]]
            SHStorm += [maxlat]
        if NHStormBounds[0] <= maxlat <= NHStormBounds[1]:
            ###NHStorm+=[mx[i]]
            NHStorm += [maxlat]
        if TropicBounds[0] <= maxlat <= TropicBounds[1]:
            ###Tropics+=[mx[i]]
            Tropics +=[maxlat]

     #If exactly one latitude is a peak between -70 and -40, it's the SH storm track
    if len(SHStorm)== 1:
        SHStorm = SHStorm[0]
        checked_lats[0] = SHStorm ###lat[SHStorm]
        checked_vals[0] = imaxes[np.argwhere(mx == SHStorm)] ###data[SHStorm]
    else:
        checked_lats[0] = 1.e20
        checked_vals[0] = 1.e20
    
       #If exactly one latitude is a peak between 40 and 70, it's the NH storm track
    if len(NHStorm)== 1:
        NHStorm = NHStorm[0]
        checked_lats[4] = NHStorm 
        checked_vals[4] =  imaxes[np.argwhere(mx == NHStorm)]
    else:
        checked_lats[4] = 1.e20
        checked_vals[4] = 1.e20
     
       #If exactly one latitude is a peak between -15 and 15, it's the tropical peak
    if len(Tropics)== 1:
        Tropics = Tropics[0]
        checked_lats[2] = Tropics###lat[Tropics]
        checked_vals[2] =  imaxes[np.argwhere(mx == Tropics)]###data[Tropics]
    else:
        checked_lats[2] = 1.e20
        checked_vals[2] = 1.e20

    #Now mask
    checked_lats = np.ma.masked_where(checked_lats>1.e10,checked_lats)
    checked_vals = np.ma.masked_where(checked_lats>1.e10,checked_vals)
    
    #If midpoints is False, return two 5d arrays
    if not midpoints:
        return checked_lats,checked_vals
    else:
        interpl = interp1d(lat,np.array(data),"cubic")
        halfmax = np.ma.zeros(6)
        vals2= np.append(np.append(data[0],checked_vals),data[-1])
        locs2 = np.append(np.append(data.getLatitude()[0],checked_lats),data.getLatitude()[-1])
        mid = 0.5*(vals2[1:]+vals2[:-1])
        for i in range(len(mid)):

            if mid[i]>1.e10:
                halfmax[i] = 1.e20
            else:
                closest = lambda x: np.abs(interpl(x) - mid[i])
                halfmax[i] = fminbound(closest,locs2[i],locs2[i+1])
                
        halfmax = np.ma.masked_where(halfmax>1.e10,halfmax)
        ret_locs = np.ma.zeros(11)
        ret_locs[::2] = halfmax
        ret_locs[1::2] = checked_lats
        return ret_locs,checked_vals

def thermo_and_dynamic(a,debug = False,own_bounds = True):
    """ Calculate thermodynamic and dynamic indicators for each time step. If own_bounds is True, calculate reasonable bounds on the location of wet and dry zones using model climatology"""
    
    nt = a.shape[0]
    timeax = a.getTime()
    loc = np.ma.zeros((nt,5))
    vals = np.ma.zeros((nt,5))
    mx,mn,imaxes,imns = find_all_peaks(MV.average(a,axis=0),return_maxmin=True)
    
    if own_bounds:
        #Are there exactly 3 maxima?
        if len(mx)!=3:
            if debug:
                print "there are not 3 maxima"
            #own_bounds = False
            mx = [-60,0,60]
        #get rid of inflection points at poles
        index_min = np.argmin(mn)
        if mn[index_min]<-85.:
            if debug:
                print "got rid of clim min at "+str(mn[index_min])
            mn = np.delete(mn,index_min)
        index_max = np.argmax(mn)
        if mn[index_max] > 85.:
            if debug:
                print "got rid of clim min at "+str(mn[index_max])
            mn = np.delete(mn,index_max)
        #Are there exactly 2 minima?
        if len(mn)!=2:
            if debug:
                print "not 2 minima"
            mn = [-20,20]

    #Check that they occur in the right order
    if own_bounds:
        xavg = np.zeros(5)
        xavg[1::2] = mn
        xavg[::2] = mx
        if len(np.where(xavg == sorted(xavg))[0]) != len(xavg):
            own_bounds = False
            if debug:
                print "not in the right order"

    
    if own_bounds:
        #Bounds are halfway between climatological peaks and troughs
        diffs = np.diff(np.append(np.append(-90.,xavg),90.))/2.
        bounds=lambda i:[xavg[i]-diffs[i],xavg[i]+diffs[i+1]]

        #Bounds on reasonable latitudes are determined by model climatology
        SHStormBounds=bounds(0)
        SHDryBounds = bounds(1)
        TropicBounds = bounds(2)
        NHDryBounds = bounds(3)
        NHStormBounds = bounds(4)
        
    else:
        #If there are not exactly 5 local extrema in the climatology, extract the ones we think may correspond to reasonable wet/dry zones
        
        SHDryBounds = [-45,-15]
        NHDryBounds = [15,45]
        SHStormBounds = [-70,-40]
        NHStormBounds = [40,70]
        TropicBounds = [-15,15]
       
    

    if debug:
        print "own bounds is "+str(own_bounds)
        print "SHStormBounds = " + str(SHStormBounds)
        print "SHDryBounds = " + str(SHDryBounds)
        print "TropicBounds = " + str(TropicBounds)
        print "NHDryBounds = " + str(NHDryBounds)
        print "NHStormBounds = " + str(NHStormBounds)

    for i in range(nt):
        
        peaks,values = get_major_features(a[i],midpoints=False,\
                                 SHStormBounds=SHStormBounds,\
                                 SHDryBounds = SHDryBounds,\
                                 TropicBounds =TropicBounds,\
                                 NHDryBounds =NHDryBounds ,\
                                 NHStormBounds =NHStormBounds)
        loc[i]=peaks
        vals[i] = values
        

    lax = cdms.createAxis(np.ma.average(loc,axis=0))
    lax.designateLatitude()

    lax2 = cdms.createAxis(np.ma.average(loc,axis=0))
    lax2.designateLatitude()

    loc = MV.array(loc)
    loc.setAxis(0,timeax)
    loc.setAxis(1,lax)
    loc.id = "dynamic"
    loc.units = a.getLatitude().units
    loc.name = loc.id

    vals = MV.array(vals)
    vals.setAxis(0,timeax)
    vals.setAxis(1,lax2)
    vals.id = 'thermodynamic'
    vals.name = vals.id
    if "units" in a.attributes.keys():
        vals.units = a.units

    return loc,vals

def get_2d_writename(fname,season,forcing=None):
   
    if fname.find("observations")<0:
        expt = fname.split(".")[2]
        if forcing is None:
            prdir = "/kate/PR/"+expt+"/"
        else:
            prdir = "/kate/PR/"+expt+"/"+forcing+"/"

        procdir = fname.split("cmip5.")[0]
        writefname = fname.replace(procdir,prdir).replace(".xml","."+season+".nc")
    else:
        expt = "OBS"
        if forcing is None:
            prdir = "/kate/PR/"+expt+"/"
        else:
            prdir = "/kate/PR/"+expt+"/"+forcing+"/"
        shortname = fname.split("/")[-1]
        writefname = prdir + shortname.replace("nc",season+".nc")
    
    if not expt in os.listdir("/kate/PR/"):
        cmd = "mkdir  /kate/PR/"+expt+"/"
        os.system(cmd)

    if forcing is not None:
        if not forcing in os.listdir("/kate/PR/"+expt+"/"):
            cmd = "mkdir "+prdir
            os.system(cmd)
    
    return writefname

def get_2d_indicators(fname,season,criteriaarg=[.99,None],plot=False,forcing = None):
    """ Given a filename, obtain D,T indicators and save them """
    writefname = get_2d_writename(fname,season,forcing=forcing)
    if writefname in glob.glob(writefname.split("cmip5.")[0]+"*"):
        print "removing "+writefname
        cmd = "rm "+writefname
        os.system(cmd)
    fw = cdms.open(writefname,"w")
    

    f = cdms.open(fname)
    rawdata = f("pr")
    f.close()
    data = prep_data(rawdata,season,criteriaarg=criteriaarg)
   
    D,T = thermo_and_dynamic(data)
        
    fw.write(D)
    fw.write(T)
    fw.close()
    if plot:
        ax,cax = zonalplot(data)
        nt = data.shape[0]
        for i in range(nt):
            x = D[i].asma()
            y = T[i].asma()
            ax.plot(x,y,"o",color = cm.RdYlBu(i/float(nt)))
        if fname.find("observations")<0:
            expt = fname.split(".")[2]
        else:
            expt = "OBS"
        picdir = "/export/marvel1/HistoricalMisc/IndicatorPix/"

        if expt != "OBS":
            ax.set_title(fname.split("cmip5.")[1].split(".mo")[0])
        else:
            ax.set_title(fname.split("_")[1])
        ax.set_ylabel("Smoothed total cloud fraction")
        picname = writefname.split("/")[-1].replace(".nc",".png")
        plt.savefig(picdir+picname)
        plt.clf()
    

    
def generate_all_2d(experiment,season,search_string="*",plot=False):
    basedir =  "/work/cmip5/"+experiment+"/atm/mo/pr/"
    badfiles = open("BAD_2D_030314_"+experiment+".txt","w")
    models = sorted(np.unique([x.split(".")[1] for x in os.listdir(basedir)]))
    for mo in models:
        ens = sorted(get_ensemble(basedir,mo,search_string))
        for fname in ens:
            print fname
            try:
                get_2d_indicators(fname,season,plot = plot)
            except:
                print "*************PROBLEM************************"
                print "*     skipping"+fname+"*"
                print "*******************************************"

                badfiles.write(fname+"\n")
    badfiles.close()
                
                
#########################################
#                                       #
#         Visual checks                 #  
#                                       #
#########################################
def zonalplot(data,cmap = cm.RdYlBu,colorbar=True,ax=None):
    if ax is None:
        ax = plt.subplot(111)
    years = [x.year for x in data.getTime().asComponentTime()]
    nt = data.shape[0]
    lat = data.getLatitude()[:]
    for i in range(nt):
        ax.plot(lat,data[i].asma(),color = cmap(i/float(nt)),linewidth=3)
    ax.set_xlabel("Latitude")
    ax.set_xlim(-90,90)
    ax.set_xticks(np.arange(-90,90+30,30))
    if colorbar:
        #ax = plt.gca()
        boundaries=years
        divider = make_axes_locatable(ax)
        cax=divider.append_axes("right",size="10%",pad=0.4)
        orient = "vertical"
        norm = mpl.colors.Normalize(vmin=min(years),vmax = max(years))
        cb1 = mpl.colorbar.ColorbarBase(cax,cmap = cm.RdYlBu,norm=norm,boundaries = boundaries,orientation=orient)
        cax.set_title("Year")
        return ax,cax
    else:
        return ax

def individual_plot(data,**kwargs):
    yrs = [x.year for x in data.getTime().asComponentTime()]
    for i in range(data.shape[1]):
        plt.subplot(data.shape[1],1,i+1)
        plt.plot(yrs,data[:,i].asma(),**kwargs)



###### Special Collections of forcings ######        
        
Ant = ["*GISS*.r*i1p109.*",\
           "*GISS*.r*i1p309.*",\
           "*CCSM4*.r[1,4,6]i1p11.*",\
           "*CESM1-CAM5*.r[1,4,6]i1p11.*"]

AA = ["*GISS*.r[1-5]i1p310.*","*CCSM4*.r[1,4,6]i1p10.*","*CESM1-CAM5*.r[1,4,6]i1p10.*"]
      
       
#### Single forcing experiments ######            
            
Vl = [ "*GISS-E2*.r[1-5]i1p[13]03.*","*CCSM4*.r[1,4,6]i1p17.*","*CESM1-CAM5*.r[1,4,6]i1p17.*"]

Sl = [ "*GISS*r[1-5]i1p[13]02*","*CCSM4*.r[1,4,6]i1p16.*","*CESM1-CAM5*.r[1,4,6]i1p17.*"]

Oz = ["*GISS*i1p105*","*CCSM4*.r[1,4,6]i1p14.*","*CESM1-CAM5*.r[1,4,6]0i1p14.*"]

LU = [ "*GISS*.r[1-5]i1p104.*","*CCSM4*.r[1,4,6]i1p13.*","*CESM1-CAM5*.r[1,4,6]i1p13.*"]


def historicalMiscIndicators():

    
###### Special Collections of forcings ######        
    
    Ant = ["*GISS*.r*i1p109.*",\
               "*GISS*.r*i1p309.*",\
               "*CCSM4*.r[1,4,6]i1p11.*",\
               "*CESM1-CAM5*.r[1,4,6]i1p11.*"]
    
    AA = ["*GISS*.r*i1p107.*","*GISS*.r[1-5]i1p310.*","*CCSM4*.r[1,4,6]i1p10.*","*CESM1-CAM5*.r[1,4,6]i1p10.*"]
    
       
#### Single forcing experiments ######            
            
    Vl = [ "*GISS-E2*.r[1-5]i1p[13]03.*","*CCSM4*.r[1,4,6]i1p17.*","*CESM1-CAM5*.r[1,4,6]i1p17.*"]
    
    Sl = [ "*GISS*r[1-5]i1p[13]02*","*CCSM4*.r[1,4,6]i1p17.*","*CESM1-CAM5*.r[1,4,6]i1p17.*"]
    
    Oz = ["*GISS*i1p105*","*CCSM4*.r[1,4,6]i1p14.*","*CESM1-CAM5*.r[1,4,6]i1p14.*"]
    
    LU = [ "*GISS*.r[1-5]i1p104.*","*CCSM4*.r[1,4,6]i1p13.*","*CESM1-CAM5*.r[1,4,6]i1p13.*"]
    ForcingExperiments = [Vl,Oz,Sl,LU,AA,Ant]
    ForcingLabels = ["Vl","Oz","Sl","LU","AA","Ant"]
    path = "/work/cmip5/historicalMisc/atm/mo/pr/"
    i=0
    for forcing in ForcingExperiments:
        for mo in forcing:
            ens = glob.glob(path+mo)
            for fname in ens:
                get_2d_indicators(fname,"DJF",plot=True,forcing=ForcingLabels[i])
        i+=1

def historicalIndicators():
    histpath = "/work/cmip5/historical/atm/mo/pr/"
    natpath = "/work/cmip5/historicalNat/atm/mo/pr/"
    ghgpath =  "/work/cmip5/historicalGHG/atm/mo/pr/"
    models = ["GISS-E2-H","GISS-E2-R","CCSM4"]#,"CESM1-CAM5"]
    ss = ["*i1p[1,3].*","*i1p[1,3].*","*","*"]
    for model,search_string in zip(models,ss):
        hist_ens = get_ensemble(histpath,model,search_string=search_string)
        for fname in hist_ens:
            print fname
            get_2d_indicators(fname,"DJF",plot=True)
        nat_ens = get_ensemble(natpath,model,search_string=search_string)
        for fname in nat_ens:
            print fname
            get_2d_indicators(fname,"DJF",plot=True)
        ghg_ens = get_ensemble(ghgpath,model,search_string=search_string)
        for fname in ghg_ens:
            print fname
            get_2d_indicators(fname,"DJF",plot=True)
 
def all_historical():
    #histpath = "/work/cmip5/historical/atm/mo/pr/"
    cpath = "/work/cmip5/piControl/atm/mo/pr/"
    #fnames = glob.glob(histpath+"*")
    #models = np.unique([x.split(".")[1] for x in fnames]).tolist()
    models = ["GISS-E2-H","GISS-E2-R","CCSM4"]
    #for x in ["GISS-E2-H","GISS-E2-R","CCSM4"]:#already done
       # models.remove(x)
    for mo in models:
     #   hist_ens = get_ensemble(histpath,mo)
      #  for fname in hist_ens:
       #     print fname
        #    get_2d_indicators(fname,"DJF",plot=False)
        c_ens = get_ensemble(cpath,mo)
        for fname in c_ens:
            print fname
            get_2d_indicators(fname,"DJF",plot=True)

def ensemble_average(path,start=None,stop=None,anoms = True):
    if start is None:
        start = cdtime.comptime(1900,1,1)
    if stop is None:
        stop = cdtime.comptime(2005,12,31)
    fnames = glob.glob(path+"*")
    models = np.unique([x.split(".")[1] for x in fnames]).tolist() 
    L = len(models)
    nt = (stop.year - start.year)+1
    Da = MV.zeros((L,nt,5))+1.e20
    Ta =  MV.zeros((L,nt,5))+1.e20
    j=0
    for mo in models:
        print mo
        ens = get_ensemble(path,mo)
        Le = len(ens)
        De = MV.zeros((Le,nt,5))+1.e20
        Te = MV.zeros((Le,nt,5))+1.e20
        i = 0
        for fname in ens:
            print fname
            f = cdms.open(fname)
            Ddata= f("dynamic",time=(start,stop))
            Tdata = f("thermodynamic",time=(start,stop))
            print Ddata.shape
#            cdutil.setTimeBoundsYearly(Ddata)
 #           cdutil.setTimeBoundsMy(Tdata)
            if anoms:
               
                Ddata = cdutil.DJF.departures(Ddata)
               
                Tdata = cdutil.DJF.departures(Tdata)
               
            if Ddata.shape[0] == nt:
                De[i]=Ddata
            else:
                print nt
                print Ddata.shape
            if Tdata.shape[0] == nt:
                Te[i] = Tdata
            f.close()
            i+=1
        Da[j] = MV.average(De,axis=0)
        Ta[j] = MV.average(Te,axis=0)
        j+=1
    Da = MV.masked_where(Da>1.e10,Da)
    print "masked d"
    Ta = MV.masked_where(Ta>1.e10,Ta)
    print "masked T"
    D = MV.average(Da,axis=0)
    print "averaged d"
    T = MV.average(Ta,axis=0)
    print "averaged T"
    D.setAxis(0,Ddata.getTime())
    T.setAxis(0,Tdata.getTime())
    return D,T


def get_dt_ensemble(path,start=None,stop=None,anoms = True,ccsm_giss=False):
    if start is None:
        if path.find("hist")>=0:
            start = cdtime.comptime(1900,1,1)
        else:
            start = 0
    if stop is None:
        if path.find("hist")>=0:
            stop = cdtime.comptime(2005,12,31)
        else:
            stop = 200
    fnames = glob.glob(path+"*")
    if ccsm_giss:
        models = ["GISS-E2-H p1","GISS-E2-R p1","CCSM4"]
    else:
        models = np.unique([x.split(".")[1] for x in fnames]).tolist() 
    L = 0
    for mo in models:

        ens = get_ensemble(path,mo)
        L+=len(ens)
    print L
    if type(start) == type(cdtime.comptime(1,1,1)):
        nt = (stop.year - start.year)+1
    else:
        nt = stop-start
    Da = MV.zeros((L,nt,5))+1.e20
    Ta =  MV.zeros((L,nt,5))+1.e20
    j=0
    labels = []
    
    for mo in models:
        print mo
        ens = get_ensemble(path,mo)
        
        for fname in ens:
            rip = fname.split(".")[3]
            labels +=[mo+"."+rip]
            print mo+"."+rip
            f = cdms.open(fname)
            if type(start) == type(cdtime.comptime(1,1,1)):
                Ddata= f("dynamic",time=(start,stop))
                Tdata = f("thermodynamic",time=(start,stop))
            else:
                Ddata = f("dynamic")[start:stop]
                Tdata = f("thermodynamic")[start:stop]
            f.close()
#            cdutil.setTimeBoundsYearly(Ddata)
 #           cdutil.setTimeBoundsMy(Tdata)
            if anoms:
               
                Ddata = cdutil.DJF.departures(Ddata)
               
                Tdata = cdutil.DJF.departures(Tdata)
               
            if Ddata.shape[0] == nt:
                Da[j]=Ddata
            else:
                print "NOT THE RIGHT SHAPE"
                print nt
                print Ddata.shape
            if Tdata.shape[0] == nt:
                Ta[j] = Tdata
            

            j+=1
    Da = MV.masked_where(Da>1.e10,Da)
    
    Ta = MV.masked_where(Ta>1.e10,Ta)
    
    modaxis = cdms.createAxis(np.arange(Da.shape[0]))
    modaxis.models = str(labels)
    modaxis.id = 'models'
    
    Da.setAxis(0,modaxis)
    Ta.setAxis(0,modaxis)

    Da.setAxis(1,Ddata.getTime())
    Ta.setAxis(1,Tdata.getTime())
    return Da,Ta

def project_historical(D,T,De,Te):
    dsolver = Eof(D)
   
    tsolver = Eof(T)
    
    P_d = MV.zeros(De.shape[:-1])+1.e20
    P_t = MV.zeros(Te.shape[:-1])+1.e20
    
    to_del = []
   
    for i in range(De.shape[0]):
        Dproj = MV.array(De[i].filled(fill_value=0.))
        Tproj = MV.array(Te[i].filled(fill_value=0.))
       
        Dproj.setAxisList(De[i].getAxisList())
        Tproj.setAxisList(Te[i].getAxisList())
        

        P_d[i] = (dsolver.projectField(Dproj)*get_orientation(dsolver))[:,0]
        if len(np.where(P_d[i] ==0.)[0]) == len(P_d[i]):
            to_del +=[i]
        P_t[i] = (tsolver.projectField(Tproj)*get_orientation(tsolver))[:,0]


    #Delete to_del, create new axes
    P_d = np.delete(P_d,to_del,axis=0)
    P_t = np.delete(P_t,to_del,axis=0)
    modaxis = cdms.createAxis(np.arange(P_d.shape[0]))
    labels =  np.delete(eval(De.getAxis(0).models),to_del).tolist()
    modaxis.models = str(labels)
    modaxis.id = 'models'
    
    P_d.setAxis(0,modaxis)
    P_t.setAxis(0,modaxis)
    
    P_d.setAxis(1,De.getTime())
    P_t.setAxis(1,Te.getTime())
    
   # P_d.setAxisList(De.getAxisList()[:-1])
    #P_t.setAxisList(Te.getAxisList()[:-1])
    
    return P_d,P_t

def project_control(D,T,Dc,Tc):
    Pdc,Ptc = project_historical(D,T,Dc,Tc)
    models = eval(Pdc.getAxis(0).models)

    nmods,nyears = Pdc.shape
    to_del=[]
    
    for i in range(nmods):
        data = Pdc[i]
        zrs = np.where(data==0.)[0]
        if len(zrs) == len(data):
            to_del+=[i]
        
    mods = np.delete(models,to_del).tolist()
    Pdc = np.delete(Pdc,to_del,axis=0).flatten()
    Ptc = np.delete(Ptc,to_del,axis=0).flatten()
    tax = cdms.createAxis(np.arange(14,14+Pdc.shape[0]*365,365.))
    tax.designateTime()
    tax.setCalendar(4113)
    tax.units= 'days since 0001-01-01'
    tax.models = str(mods)

    Pdc.setAxis(0,tax)
    Ptc.setAxis(0,tax)
    return Pdc,Ptc
    

# def project_control(D,T,Dc,Tc):
#     dsolver = Eof(D)
   
#     tsolver = Eof(T)
    
    

#     P_d = MV.zeros(Dc.shape[:-1])+1.e20
#     P_t = MV.zeros(Tc.shape[:-1])+1.e20
    
   
#     for i in range(Dc.shape[0]):
#         Dproj = MV.array(Dc[i].filled(fill_value=0.))
#         Tproj = MV.array(Tc[i].filled(fill_value=0.))
#         Dproj.setAxisList(Dc[i].getAxisList())
#         Tproj.setAxisList(Tc[i].getAxisList())
       

#         P_d[i] = (dsolver.projectField(Dproj)*get_orientation(dsolver))[:,0]
#         P_t[i] = (tsolver.projectField(Tproj)*get_orientation(tsolver))[:,0]
#     P_d = P_d.reshape(P_d.shape[0]*P_d.shape[1])
#     P_t = P_t.reshape(P_t.shape[0]*P_t.shape[1])

#     tax = cdms.createAxis(np.arange(14,14+P_d.shape[0]*365,365.))
#     tax.designateTime()
#     tax.setCalendar(4113)
#     tax.units= 'days since 0001-01-01'
    
#     P_d= MV.masked_where(P_d==0,P_d)
#     P_t = MV.masked_where(P_t ==0, P_t)
    
#     P_d.setAxis(0,tax)
#     P_t.setAxis(0,tax)
    
#     return P_d,P_t

def justsignal(Pd,Pt,Pdc,Ptc,start = None,stop=None):
    if start is None:
        start = Pd.getTime().asComponentTime()[0]
    if stop is None:
        stop = Pd.getTime().asComponentTime()[-1]
    
    Pd = Pd(time=(start,stop))
    Pt = Pt(time=(start,stop))
    
    nyears = Pd.shape[1]

    Dsig = genutil.statistics.linearregression(Pd,axis=1,nointercept=1)*3650.
    Tsig = genutil.statistics.linearregression(Pt,axis=1,nointercept=1)*3650.
    
    
    Dsn = MV.array(Dsig)
    Tsn = MV.array(Tsig)

    Dsn.setAxis(0,Pd.getAxis(0))
    Tsn.setAxis(0,Pt.getAxis(0))
    return Dsn, Tsn

def signal_to_noise(Pd,Pt,Pdc,Ptc,start = None,stop=None):
    if start is None:
        start = Pd.getTime().asComponentTime()[0]
    if stop is None:
        stop = Pd.getTime().asComponentTime()[-1]
    
    Pd = Pd(time=(start,stop))
    Pt = Pt(time=(start,stop))
    
    nyears = Pd.shape[1]

    Dsig = genutil.statistics.linearregression(Pd,axis=1,nointercept=1)*3650.
    Tsig = genutil.statistics.linearregression(Pt,axis=1,nointercept=1)*3650.
    
    Dnoise = float(np.ma.std(get_slopes(Pdc,nyears)))
    Tnoise = float(np.ma.std(get_slopes(Ptc,nyears)))

    Dsn = MV.array(Dsig/Dnoise)
    Tsn = MV.array(Tsig/Tnoise)

    Dsn.setAxis(0,Pd.getAxis(0))
    Tsn.setAxis(0,Pt.getAxis(0))
    return Dsn, Tsn
    
def histogram_hist_check(Dsn,Tsn):
    normed = False
    
    models = eval(Dsn.getAxis(0).models)
    giss_r = np.where([x.find("GISS-E2-R.")>=0 for x in models])[0]
    giss_h = np.where([x.find("GISS-E2-H.")>=0 for x in models])[0]
    ccsm4 = np.where([x.find("CCSM4.")>=0 for x in models])[0]
    
    #Remove missing data (missing bc historical run too short)
    I=np.where(Dsn==0)[0]
    Dsni = np.delete(Dsn,I)

    plt.subplot(121)
    freq,bins,patches = plt.hist(Dsni,25,normed = normed,zorder=1,linewidth=4,color = cm.RdYlBu(0.),label="CMIP5",ec = cm.RdYlBu(0.))#, weights = np.zeros_like(Dsni)+1./Dsni.size)

    ccsm4data = Dsn.asma()[ccsm4]
    plt.hist(Dsn.asma()[ccsm4],bins=bins,normed=normed,color = cm.RdYlBu(.3),ec=cm.RdYlBu(.3),label = "CCSM4")
    
    ghdata = Dsn.asma()[giss_h]
    plt.hist(ghdata,normed=normed,bins=bins,color = cm.RdYlBu(.7),ec=cm.RdYlBu(.7),label="GISS-E2-H")#, weights = np.zeros_like(ghdata)+1./ghdata.size)

    grdata = Dsn.asma()[giss_r]
    plt.hist(Dsn.asma()[giss_r],bins=bins,normed=normed,color = cm.RdYlBu(.9),ec=cm.RdYlBu(.9),label="GISS-E2-R")#, weights = np.zeros_like(grdata)+1./grdata.size)
    leg=plt.legend(loc=0)
    leg.set_frame_on(False)
    plt.title("Dynamic")
    plt.xlabel("S/N")
    plt.ylabel("Number of models")


    I=np.where(Tsn==0)[0]
    Tsni = np.delete(Tsn,I)
    plt.subplot(122)
    freq,bins,patches = plt.hist(Tsni,25,normed = normed,zorder=1,linewidth=4,color = cm.RdYlBu(0.),ec = cm.RdYlBu(0.))
    plt.hist(Tsn.asma()[ccsm4],bins=bins,normed=normed,color = cm.RdYlBu(.3),ec=cm.RdYlBu(.3))
    
    plt.hist(Tsn.asma()[giss_h],bins=bins,normed=normed,color = cm.RdYlBu(.7),ec=cm.RdYlBu(.7))
    plt.hist(Tsn.asma()[giss_r],bins=bins,normed=normed,color = cm.RdYlBu(.9),ec=cm.RdYlBu(.9))
    plt.title("Thermodynamic")
    plt.xlabel("S/N")
    plt.ylabel("Number of models")
    
    

class SN():
    def __init__(self,D=None,T=None,Dc=None,Tc=None,Da=None,Ta = None):
        if D is None:
            D,T = ensemble_average("/kate/PR/historical/")
        self.D = D
        self.T = T
        if Da is None:
            self.Da,self.Ta = get_dt_ensemble("/kate/PR/historical/")
        self.dsolver = Eof(D)
        self.tsolver = Eof(T)
        if Dc is None:
            self.Dc,self.Tc = get_dt_ensemble("/kate/PR/piControl/",start =0,stop=200)
        else:
            self.Dc=Dc
            self.Tc=Tc
        self.sn_dict = None
        self.Pdc,self.Ptc=project_control(self.D,self.T,self.Dc,self.Tc)

    def get_sn(self,forcing,start=None,stop=None,ccsm_giss=True):
        if forcing.find("Nat")>=0:
            hpath = "/kate/PR/historicalNat/"
        elif forcing.find("hist")>=0:
            hpath = "/kate/PR/historical/"
        elif forcing.find("GHG")>=0:
            hpath = "/kate/PR/historicalGHG/"
        else:
            hpath = "/kate/PR/historicalMisc/"+forcing+"/"
        print hpath
        De,Te = get_dt_ensemble(hpath,ccsm_giss=ccsm_giss,start=start,stop=stop)
        print De.shape
        Pd,Pt=project_historical(self.D,self.T,De,Te)
        Dsn,Tsn = signal_to_noise(Pd,Pt,self.Pdc,self.Ptc)
        return Dsn,Tsn
    def compare_all_sn(self):
        forcings = ["hist","Nat","GHG","AA","Ant","Oz","LU","Sl","Vl"]
        
        sn_dict = {}
        for forcing in forcings:
            forcingdir = {}
            Dsn,Tsn = self.get_sn(forcing)
            models = eval(Dsn.getAxis(0).models)
            giss_r = np.where([x.find("GISS-E2-R")>=0 for x in models])[0]
            giss_h = np.where([x.find("GISS-E2-H")>=0 for x in models])[0]
            ccsm4 = np.where([x.find("CCSM4")>=0 for x in models])[0]

            forcingdir["GISS-E2-R"] = [Dsn.asma()[giss_r],Tsn.asma()[giss_r]]
            forcingdir["GISS-E2-H"] = [Dsn.asma()[giss_h],Tsn.asma()[giss_h]]
            forcingdir["CCSM4"] = [Dsn.asma()[ccsm4],Tsn.asma()[ccsm4]]
            forcingdir["ALL"] = [Dsn.asma(),Tsn.asma()]
            sn_dict[forcing] = forcingdir
        self.sn_dict = sn_dict
    

def get_stats(sn_dict,model,typ="D",return_interval=None):
    forcings = ["AA","GHG","LU","Oz","Sl","Vl"]

    if lower(typ).find("d")==0:
        idx = 0
    else:
        idx = 1

    the_mean = np.average(sn_dict["hist"][model][idx]) - sum([np.average(sn_dict[forcing][model][idx]) for forcing in forcings])

    N = [len(sn_dict["hist"][model][idx])]+ [len(sn_dict[forcing][model][idx]) for forcing in forcings]
    
    df = sum(N)-(len(forcings)+1)
    
    prefac = np.sqrt(np.sum([1/float(N[i]) for i in range(len(N))]))
    samps = [sn_dict["hist"][model][idx]]+ [sn_dict[fcg][model][idx] for fcg in forcings]
    STD = prefac*pooled_std(samps)
    
    
    tval = the_mean/STD
    pval = stats.t.sf(np.abs(tval),df)*2

    
    
    confi = STD*np.array(stats.t.interval(return_interval,df))+the_mean
    
    if return_interval is None:
        return the_mean,STD,the_mean/STD,pval
    else:
        return the_mean,confi

def get_ztest(sn_dict,model,typ="D"):
    sn_dict = sn_dict.copy()
    xs = 0
    forcings = ["AA","GHG","LU","Oz","Sl","Vl"]
    
    if lower(typ).find("d")==0:
        idx = 0
    else:
        idx = 1
    #print "Using "+str(idx)
        

    samps =  [sn_dict[fcg][model][idx] for fcg in forcings]
    trunc = min([len(x) for x in samps])
    xh = sn_dict["hist"][model][idx]
   

    xs = samps[0].copy()[:trunc]
            
   # print "samps[0] = "+str(xs)
    for i in range(len(samps))[1:]:
        xs += samps[i][:trunc]
    
    nf = float(len(samps))
    ns = float(len(xs))
    nh = float(len(xh))
    print "(Nf,nh,ns)="+str((nf,nh,ns))
  
    xhbar = np.average(xh)
    xsbar = np.average(xs)
    
    num = xhbar - xsbar
    print "xhbar = "+str(xhbar)
    print "xsbar = "+str(xsbar)
    
    denom = np.sqrt(1./nh + nf/ns)
    print "denom=" +str(denom)
    z = num/denom

    pval = stats.norm.sf(np.abs(z))*2
    return z,pval

def get_bootstrap(sn_dict,model,typ="D",dists=False):
    forcings = ["AA","GHG","LU","Oz","Sl","Vl"]

    if lower(typ).find("d")==0:
        idx = 0
    else:
        idx = 1
    samps =  [sn_dict[fcg][model][idx] for fcg in forcings]
    test = list(itertools.product(*samps))
    
    sumdist=[sum(x) for x in test]
    histdist=list(sn_dict["hist"][model][idx])
    
    if dists:
        return sumdist,histdist
    else:
        tstat,pval=stats.ttest_ind(sumdist,histdist,equal_var=False)
        
        return np.average(sumdist),np.average(histdist),float(tstat),pval
    
def stats_tables(sn_dict,typ,func = get_stats, f=None):
    for model in ["CCSM4","GISS-E2-H","GISS-E2-R","ALL"]:
        
        X = np.round(func(sn_dict,model,typ=typ),3).tolist()
        #t = np.round(get_stats(sn_dict,model,typ="t"),3).tolist()
        if f == None:
            print model +" & "+ join(map(str,X)," &")+"\\\\"
        else:
            f.write( model +" & "+ join(map(str,X)," &")+"\\\\")

def two_sided_t_test(sn_dict,model):
    pass

def make_colormap_historicalMisc():
    Lu = cm.BrBG(0.1)
    Vl = cm.Reds(.9)
    GHG = cm.Greens(.9)
    CO2 = cm.Greens(.5)
    CH4 = cm.spring(0.1)
    Sl = cm.RdYlBu(.4)
    Oz = cm.Blues(.8)
    CFCs= cm.Blues(1.)
    Ant = cm.Accent(.7)
    hist = "k"
    Nat = cm.Greens(.3)
    AA = cm.YlOrRd(.6)
    SUM = cm.Purples(.8)
    DIFF = "m"
    clist = [AA,GHG,Lu,Oz,Sl,Vl,hist,Nat,Ant,SUM,DIFF,CFCs,CO2,CH4]
    cmap = mpl.colors.ListedColormap(clist)
    return cmap
def get_color(forcing):
    cmap = make_colormap_historicalMisc()
    keys = ["AA","GHG","LU","Oz","Sl","Vl","hist","Nat","Ant","SUM","DIFF","CFCs","CO2","CH4"]
    L = float(len(keys))
    i = map(lower,keys).index(lower(forcing))
    return cmap(i/L)

def get_marker(forcing):
    keys = ["AA","GHG","LU","Oz","Sl","Vl","hist","Nat","Ant","SUM","DIFF"]
    markers = ["o","s","d",  "*","h", "8",  "p",  "v", "^",   "H", "+"]
    i = map(lower,keys).index(lower(forcing))
    return markers[i]
    
def vectorplot(sn_dict,model,show_all=False,forcings = None,legend=True):
    if forcings is None:
        forcings = sn_dict.keys()
    for forcing in forcings:
        dummy = []

        color = get_color(forcing)
        if forcing != "SUM":
            D,T = sn_dict[forcing][model]
            if show_all:
                for dx,dy in zip(D,T):    
                    plt.plot([0,dx],[0,dy],color=color,linestyle="-",linewidth=1,zorder=1)
        else:
            histf=["AA","GHG","LU","Oz","Sl","Vl"]
            D = sum([np.average(sn_dict[fcing][model][0]) for fcing in histf])
            T =  sum([np.average(sn_dict[fcing][model][1]) for fcing in histf])

        dx = np.average(D)
        dy = np.average(T)



        plt.arrow(0,0,dx,dy,color=color,width=.2,length_includes_head=False,zorder=2)
        #plt.plot([dx,dx],[0,dy],"--",color = color,linewidth=3)
        #plt.plot([0,dx],[dy,dy],"--",color = color,linewidth=3)
        item,=plt.plot([0],[0],"s",color=color,mec =color, label=forcing)

        dummy +=[item]
    plt.xlim(-5,5)
    plt.ylim(-5,5)
    if legend:
        leg = plt.legend(numpoints=1,markerscale=2,ncol=1,loc=2)
    [x.set_visible(False) for x in dummy]
    plt.xlabel("Dynamic S/N")
    plt.ylabel("Thermodynamic S/N")
    plt.title(model)
    plt.axhline(0,linestyle=":",color="k")
    plt.axvline(0,linestyle=":",color="k")
    conf = 1.96
    C  = plt.Rectangle((-10,-conf),width=20,height=2*conf,color=cm.Reds(.8),alpha=.2)
    plt.gca().add_patch(C)
    C2  = plt.Rectangle((-conf,-10),height=20,width=2*conf,color=cm.Blues(.5),alpha=.2)
    plt.gca().add_patch(C2)

def plot_three_models(sn_dict,compare_all=True,forcings = None,conf=.9):
    if forcings is None:
        if compare_all:
            forcings=["AA","Sl","Oz","LU","Vl","GHG","hist","SUM"]
        else:
            forcings = ["Ant","Nat","hist"]
    leg=[True,False,False]
   # models = ["CCSM4", "GISS-E2-H","GISS-E2-R"]#,"ALL"]
    models = ["GISS-E2-H","GISS-E2-R","CCSM4"]
    letters = ["(a): ", "(b): ","(c): ","(d): "]
    for i in range(3):
        plt.subplot(3,1,i+1)
        forcingplot(sn_dict,models[i],show_all=True, forcings=forcings,legend=leg,letter = letters[i],conf=conf)

def eofs_with_color(X,cmap=cm.RdBu,extend = "both",title=""):
    solver = Eof(X)
#    tsolver = Eof(T)
    
    eof = solver.eofs()[0]*get_orientation(solver)


    f = cdms.open("/kate/MMA/climatologies/cmip5.MMA.piControl.pr.with_axes.nc")
    reference=MV.average(f("pr"),axis=0)*60*60*24.
    smooth = MV.array(ndimag.gaussian_filter(reference,sigma=5./2.))
    smooth.setAxisList(reference.getAxisList())
    reference = smooth
    lats = reference.getLatitude()[:]
    x,y=find_all_peaks(smooth)

    plt.plot(lats,smooth.asma(),"k")
    f.close()

    mx =1.# np.max(eof)
    mn = -1.#np.min(eof)
    col = lambda i: cmap((eof[i]-mn)/(mx-mn))
    norm = mpl.colors.Normalize(vmin = mn,vmax = mx)
    j = 0
    for j in range(5):
        xx = x[j]
        yy = y[j]
        plt.plot([xx],[yy],"o",markersize=30,color=col(j))#,mec = col(j))
    
    plt.xticks(np.arange(-90,90+30,30))
    plt.ylim(0,6.5)
    plt.ylabel("Smoothed P (mm/day)")
    neg_to_S(plt.gca())
    plt.title(title)

    divider = make_axes_locatable(plt.gca())
    ax2 = divider.append_axes("bottom",size="10%",pad=0.7)
    cb1 = mpl.colorbar.ColorbarBase(ax2,cmap=cmap,norm=norm,orientation='horizontal',extend=extend)
    ax2.set_title("EOF Loading")

def all_eofs(D,T):
    T = T*60*60*24.
    plt.subplot(221)
    eofs_with_color(D,title="(a): EOF1 (dynamic)")
    
    plt.subplot(223)
    dsolver = Eof(D)
    pcd = dsolver.pcs()[:,0]*get_orientation(dsolver)
    t = get_plottable_time(pcd)
    plt.plot(t,pcd.asma(),color=cm.gray(.8),linewidth=1)
    plt.plot(t,pcd.asma(),"o",color=cm.gray(.4),mec=cm.gray(.4))
    plt.xlabel("Year")
    plt.ylabel("Temporal Amplitude")
    plt.title("(c): PC1 (dynamic)")
    plt.xlim(1900,2006)
    plt.gca().set_yticks(np.arange(-.5,2.5,.5))

    plt.subplot(222)
    eofs_with_color(T,extend="neither",title="(b): EOF1 (thermodynamic)")
   
    plt.subplot(224)
    tsolver = Eof(T)
    pct = tsolver.pcs()[:,0]*get_orientation(tsolver)
    t = get_plottable_time(pct)
    plt.plot(t,pct.asma(),color=cm.gray(.8),linewidth=1)
    plt.plot(t,pct.asma(),"o",color=cm.gray(.4),mec=cm.gray(.4))
    plt.xlabel("Year")
    plt.ylabel("Temporal Amplitude")
    plt.title("(d): PC1 (thermodynamic)")
    plt.xlim(1900,2006)
    plt.gca().set_yticks(np.arange(-.1,.3,.1))
    
def neg_to_S(ax=None,**kwargs):
    if ax is None:
        ax = plt.gca()
    tix = ax.get_xticks()
    labels = []
    for lat in tix:
        if lat<0:
            labels+=[str(abs(int(lat)))+r'$^{\circ}$S']
        elif lat > 0:
            labels+=[str(int(lat))+r'$^{\circ}$N']
        else:
            labels+=[str(int(lat))+r'$^{\circ}$']
    ax.set_xticklabels(labels,**kwargs)
    plt.draw()
  
def plot_eofs(D,T):

    #REPLACE THIS WITH PR
#    f = cdms.open("/kate/MMA/climatologies/piControl.DJF.MMA.reference.nc")
    f = cdms.open("/kate/MMA/climatologies/cmip5.MMA.piControl.pr.with_axes.nc")
    reference=MV.average(f("pr"),axis=0)*60*60*24.
    smooth = MV.array(ndimag.gaussian_filter(reference,sigma=5./2.))
    smooth.setAxisList(reference.getAxisList())
    reference = smooth
    lats = reference.getLatitude()[:]
   # plt.subplot(211)
    plt.plot(lats,smooth.asma(),"k")
    x,y=find_all_peaks(smooth)

    
    dsolver = Eof(D)
    tsolver = Eof(T)
    
    deof = dsolver.eofs()[0]*get_orientation(dsolver)
    teof = tsolver.eofs()[0]*get_orientation(tsolver)


    for i in range(len(deof)):
        #color = cm.RdYlBu(i/4.)
        color="k"
        dx = np.ma.array(deof)[i]*20
        plt.arrow(x[i],y[i],dx,0,color=color,head_width=.2,head_length=6,length_includes_head=False)#,width=.5,length_includes_head=False)
#        plt.arrow(x[i],y[i],dx,0, head_width=np.abs(dx)*.3,color = "k")
    plt.xlim(-90,90)
    plt.ylim(0,6.7)
    #plt.subplot(212)
    #plt.plot(lats,smooth.asma(),"k")
    for i in range(len(teof)):
        dy = np.ma.array(teof)[i]#*5#*10
        plt.arrow(x[i],y[i],0,dy,color="r",head_width=6,head_length=.2,length_includes_head = False)#,width=1,head_length=.2,head_width=6,length_includes_head=False)#,width=.5,length_includes_head=False)
        
#        plt.arrow(x[i],y[i],0,dy,color=color,head_width=np.abs(dy)*.6,head_length=np.abs(dy)*.1)
    plt.xlim(-90,90)
    plt.ylim(0,6.7)
    #plt.plot(x,y,"o",markersize=20,color=cm.bone(.5),mec=cm.bone(.5))
    
    

    
    #plt.plot(x,y,"ro")
   
    #[x.set_visible(False) for x in dummy]
#if __name__=='__main__':
 #   print "hello!"
#    historicalIndicators()
  #  all_historical()



def pooled_std(samples):
    num = 0.
    denom = 0.
    k = 0.
    for samp in samples:
        s2 = np.var(samp)
        n = len(samp)
        num+=(n-1)*s2
        denom+=n
        k+=1
    return np.sqrt(num/float(denom-k))


def forcingplot(sn_dict,model,show_all=False,forcings = None,legend=True,letter="",conf=.95,all_forcings = True):
    if forcings is None:
        forcings = sn_dict.keys()
    for forcing in forcings:
        dummy = []

        color = get_color(forcing)
        if forcing in ["AA","GHG","LU","Oz","Sl","Vl","Ant","Nat","hist"]:
            D,T = sn_dict[forcing][model]
            
                #for dx,dy in zip(D,T):    
                 #   plt.plot([dx],[dy],"o",color=color,mec=color,markersize=10,zorder=2)
            xa = np.average(D)
            ya = np.average(T)
            df = len(D)-1
            confi = stats.t.interval(conf,df)[1]
            xs = confi*np.std(D)*2
            ys = confi*np.std(T)*2
                #ell = Ellipse((xa,ya),width=xs,height=ys,color=color,alpha=.4,ec=color,zorder=2)
            if show_all:
                ell = Rectangle((xa-xs/2.,ya-ys/2.),width=xs,height=ys,color=color,alpha=1,lw=2,fill=None,ec=color,zorder=2)
                plt.gca().add_patch(ell)
        else:
            if all_forcings:
                histf=["AA","GHG","LU","Oz","Sl","Vl"]
            else:
                histf = ["Ant","Nat"]

            if forcing == "SUM":
                sampsd =  [sn_dict[fcg][model][0] for fcg in histf]            
                sampst =  [sn_dict[fcg][model][1] for fcg in histf]
         
                xa = sum([np.average(x) for x in sampsd])
                ya = sum([np.average(x) for x in sampst])

                Nd =  [len(sn_dict[fcg][model][0]) for fcg in histf]
                Nt =  [len(sn_dict[fcg][model][1]) for fcg in histf]
            elif forcing == "DIFF":
                sampsd =  [sn_dict[fcg][model][0] for fcg in histf]+[sn_dict["hist"][model][0]]         
                sampst =  [sn_dict[fcg][model][1] for fcg in histf]+[sn_dict["hist"][model][1]]
         
                xa = np.average(sn_dict["hist"][model][0]) - sum([np.average(x) for x in sampsd[:-1]])
                ya = np.average(sn_dict["hist"][model][1]) -sum([np.average(x) for x in sampst[:-1]])

                Nd =  [len(sn_dict[fcg][model][0]) for fcg in histf]+[len(sn_dict["hist"][model][0])]
                Nt =  [len(sn_dict[fcg][model][1]) for fcg in histf]+[len(sn_dict["hist"][model][1])]
                
                
                

            if show_all:

                df = sum(Nd)-(len(histf)+1)
                confi = stats.t.interval(conf,df)[1]
                
                xs = confi* np.sqrt(np.sum([1/float(Nd[i]) for i in range(len(Nd))]))*pooled_std(sampsd)*2

               
                df = sum(Nt)-(len(histf)+1)
                confi = stats.t.interval(conf,df)[1]

                ys = confi* np.sqrt(np.sum([1/float(Nt[i]) for i in range(len(Nt))]))*pooled_std(sampst)*2
                
                
                #ell = Ellipse((xa,ya),width=xs,height=ys,color=color,alpha=.4,ec=color,zorder=2)
                ell = Rectangle((xa-xs/2.,ya-ys/2.),width=xs,height=ys,color=color,alpha=1,ec=color,lw=2,zorder=2,fill=None)
                plt.gca().add_patch(ell)
              

#        dx = np.average(D)
 #       dy = np.average(T)



       
        plt.plot([xa],[ya],"o",color=color,mec=color,markersize=20,zorder=2)
        
        plt.plot([xa],[ya],"o",color=color,mec=color,markersize=5, label=forcing)
       
        #item,=plt.plot([0],[0],"s",color=color,mec =color, label=forcing)

        #dummy +=[item]
    
    plt.xlim(-5,5)
    plt.ylim(-5,5)
    
    [x.set_visible(False) for x in dummy]
    plt.xlabel("Dynamic S/N")
    plt.ylabel("Thermodynamic S/N")
    plt.title(letter+model)
   # plt.axhline(0,linestyle=":",color="k")
    #plt.axvline(0,linestyle=":",color="k")
    confi_norm = stats.norm.interval(conf)[1]
    pct = str(int(100*conf))
    C  = plt.Rectangle((-10,-confi_norm),width=20,height=2*confi_norm,fill=None,hatch="\\\\",zorder=1,linewidth=0,color="k",alpha=.4,label = pct+"% CI for T noise")
    plt.gca().add_patch(C)

    C2  = plt.Rectangle((-confi_norm,-10),height=20,width=2*confi_norm,hatch="//",fill=None,zorder=1,linewidth=0,alpha=.4,label=pct+"% CI for D noise")
    plt.gca().add_patch(C2)
    
    plt.xlim(-6,6)
    plt.ylim(-6,6)
    plt.axvline(0,color="k",zorder=1)
    plt.axhline(0,color="k",zorder=1)
    if legend:
        leg = plt.gca().legend(numpoints=1,markerscale=2,ncol=1,loc=0,fontsize=10)
        #leg.set_frame_on(False)

    #Ca  = plt.Rectangle((-conf,-conf),width=2*conf,height=2*conf,color=cm.Purples(.8))
    #plt.gca().add_patch(Ca)
    
   

def histogram_plot2(Dc,Tc,Da,Ta,D,T):
    colc = cm.Blues(.6)
    colh = cm.Reds(.6)
    CONTROL = project_control(D,T,Dc,Tc)
    
    HISTORICAL = project_historical(D,T,Da,Ta)
    nyears = Da.shape[0]
    legs = []
    for i in range(2):

        plt.subplot(2,1,i+1)
        C = get_slopes(CONTROL[i],nyears)
        #print C
        H = genutil.statistics.linearregression(HISTORICAL[i],axis=1,nointercept=1)*3650.
        if i == 1:
            H = H*60.*60.*24.
            C = C*60.*60.*24.
        #print H
        nc,bc,patchc=plt.hist(C,histtype="step", normed = True,color=colc,linewidth=2,alpha=.5,label = "CMIP5 piControl trend distribution")
        nh,bh,patchh = plt.hist(H,histtype="step",normed = True,color=colh,linewidth=2,alpha=.5,label = "CMIP5 historicalr trend distribution")

        #        xc = 0.5*(bc[1:]+bc[:1])
        a = max([np.max(C),np.max(H)])
        a = a + 0.5*a
        delta = a/25.
        xc = np.arange(-a,a+delta,delta)

        muc = np.ma.average(C)
        sigc = np.ma.std(C)
        fac = 1./sigc
        pdfc = mlab.normpdf(xc,muc,sigc)
        plt.plot(xc,pdfc,color=colc,linestyle="--",linewidth=3,label = "Gaussian fit to CMIP5 control")
        
        muh = np.ma.average(H)
        sigh = np.ma.std(H)
        
        pdfh = mlab.normpdf(xc,muh,sigh)
        plt.plot(xc,pdfh,color=colh,linestyle="--",linewidth=3,label = "Gaussian fit to CMIP5 historical trends")

        #Indicate GISS and CCSM4 historical trends 

        models = eval(H.getAxis(0).models)
        
        giss_r = np.where([x.find("GISS-E2-R.")>=0 for x in models])[0]
        giss_h = np.where([x.find("GISS-E2-H.")>=0 for x in models])[0]
        ccsm4 = np.where([x.find("CCSM4.")>=0 for x in models])[0]
        
        gauss = interp1d(xc,pdfh)
        mods = [ccsm4,giss_h,giss_r]
        labels = ["CCSM4 piControl trends","GISS-E2-H piControl trends","GISS-E2-R piControl trends"]
        for j in range(3):
            test = np.array(H)[mods[j]]
            for slp in test:
                
                if slp != test[0]:
                    plt.plot([slp,slp],[0,gauss(slp)],color= cm.Purples((j+1)/3.),linewidth=3,zorder=1)
                else:
                    plt.plot([slp,slp],[0,gauss(slp)],color= cm.Purples((j+1)/3.),linewidth=3,zorder=1,label = labels[j])


        #Indicate GISS and CCSM4 control trends 
        cmods = eval(CONTROL[i].getAxis(0).models)
       
        giss_r = np.where([x.find("GISS-E2-R.")>=0 for x in cmods])[0]
        giss_h = np.where([x.find("GISS-E2-H.")>=0 for x in cmods])[0]
        ccsm4 = np.where([x.find("CCSM4.")>=0 for x in cmods])[0]
        mods = [ccsm4,giss_h,giss_r]
        labels = ["CCSM4 historical trends","GISS-E2-H historical trends","GISS-E2-R historical trends"]
        gauss = interp1d(xc,pdfc)
        for j in range(3):
            Cr = CONTROL[i].reshape(len(CONTROL[i])/200,200)
            SLOPES = MV.array(np.array(Cr)[mods[j]].flatten())
           
            tax = CONTROL[i].getTime()
            newtime = cdms.createAxis(tax[:len(SLOPES)])
            newtime.designateTime()
            newtime.units = tax.units
            newtime.setCalendar(tax.getCalendar())
            SLOPES.setAxis(0,newtime)
           
            test = get_slopes(SLOPES,nyears)
            
            if i == 1:
                 test = test*60*60.*24.
            print test
            for slp in test:
                if slp != test[0]:
                    plt.plot([slp,slp],[0,gauss(slp)],color= cm.Greens((j+1)/3.),linewidth=3,zorder=1)
                else:
                    plt.plot([slp,slp],[0,gauss(slp)],color= cm.Greens((j+1)/3.),linewidth=3,zorder=1,label = labels[j])
        
        legs+=[plt.legend(loc=0,fontsize=10)]
    axes = plt.gcf().axes
    axes[0].set_xlim(-0.3,0.3)
    axes[0].set_xlabel(r"106-year dynamic trend ($^{\circ}$Lat/decade)")
    axes[0].set_ylabel("Normalized frequency")
    
    axes[1].set_xlim(-0.05,0.05)
    axes[1].set_xlabel(r"106-year thermodynamic trend (mm/day/decade)")
    axes[1].set_ylabel("Normalized frequency")
            
        


       #np.array(H)[giss_r]
        
       
    
    
    

    

def scratch():
    path = "/work/cmip5/piControl/atm/mo/pr/"
    fnames = glob.glob(path+"*")
    models = np.unique([x.split(".")[1] for x in fnames])
    L = 0
    modlabs = []
    for mo in models:
        ens = get_ensemble(path,mo)
        L+=len(ens)
        for fname in ens:
            modlabs += [mo+"."+fname.split(".")[3]]
    modax = cdms.createAxis(range(L))
    modax.models = str(modlabs)
    return modax


def global_average_noise(vari="tas",season="DJF",trunc=True):
    files = glob.glob("/kate/T/GLOBAL_MEAN/piControl/*"+season+"*")
    f = cdms.open(files[0])
    if vari != "resid":
        data = f(vari)
    else:
        model = files[0].split(".")[1]
        sens = hydrological_sensitivity(model,season=season)
        data = f("pr") - sens*f("tas")
    #get calendar
    old_time = data.getTime()
    delta = np.median(np.diff(old_time))
    
    if trunc:
        data = data[:200]
    data = data.anom(axis=0)
    f.close()
    for fil in files[1:]:
       # print fil
        f = cdms.open(fil)
        if vari != "resid":
            if trunc:
                newdata = f(vari)[:200].anom(axis=0)
            else:
                newdata = f(vari).anom(axis=0)
        else:
            model = fil.split(".")[1]
            sens = hydrological_sensitivity(model,season=season)
            newdata = f("pr")-sens*f("tas")
            if trunc:
                newdata = newdata[:200].anom(axis=0)
            else:
                newdata = newdata.anom(axis=0)
       
        data = MV.concatenate((data,newdata),axis=0)
        f.close()
    
    nt = len(data)
    
    start = old_time[0]
    
    stop = old_time[0]+delta*nt
    new_time = cdms.createAxis(np.arange(start,stop,delta))
    new_time.designateTime()
    for att in old_time.attributes:
        setattr(new_time,att,old_time.attributes[att])
   
   
    data.setAxis(0,new_time)
    return data
    
        
def global_average_trends(vari = "tas",season="DJF",start = None,stop=None,noise=None):
    if noise is None:
        noise = global_average_noise(vari,season)
   
    tdict = {}
    if start is None:
        start = cdtime.comptime(1900,1,1)
    if stop is None:
        
        stop = cdtime.comptime(2005,12,31)
    daydiff=stop.torel("days since 1900-1-1").value - start.torel("days since 1900-1-1").value
    nyear = int(np.round(daydiff/365,2))
    ns = np.std(get_slopes(noise,nyear)) # will be in K per decade units
    
    keys = ['AA', 'Ant', 'Vl', 'LU', 'hist', 'GHG', 'Oz', 'Nat', 'Sl',"CO2","CFCs","CH4"]
    prefix = "/kate/T/GLOBAL_MEAN/"
    sn_dict={}

    for k in keys:
       # print k
        model_dict={}
        if k == "Nat":
            direc = "historicalNat"
        elif k == "GHG":
            direc = "historicalGHG"
        elif k == "hist":
            direc = "historical"
        else:
            direc = k
        path = "/kate/T/GLOBAL_MEAN/"+direc+"/"
        models = ["CCSM4","GISS-E2-R*p1","GISS-E2-H*p1","GISS-E2-R*p3","GISS-E2-H*p3"]
        all_slopes = []
        for m in models:
            slopes = []
            files = sorted(glob.glob(path+"*"+m+"*"+season+"*") )
            
            for fil in files:
                f = cdms.open(fil)
                if vari != "resid":
                    data = f(vari,time = (start,stop))
                else:
                    sens = hydrological_sensitivity(m,season=season)
                    data = f("pr",time = (start,stop)) - sens*f("tas",time = (start,stop))
                slopes += [float(genutil.statistics.linearregression(data,nointercept=1))*3650.]
                all_slopes += [float(genutil.statistics.linearregression(data,nointercept=1))*3650.]
                f.close()
            model_dict[m] = np.array(slopes)/ns
        model_dict["ALL"] = np.array(all_slopes)/ns
        sn_dict[k] = model_dict
    return sn_dict
    

def TP_dict(season,start = None,stop=None):
    tnoise = global_average_noise("tas",season)
    pnoise = global_average_noise("pr",season)
    rnoise = global_average_noise("resid",season)
    T_sn_dict = global_average_trends(vari="tas",season=season,start=start,stop=stop,noise=tnoise)
    P_sn_dict = global_average_trends(vari="pr",season=season, start = start, stop = stop,noise=pnoise)
    R_sn_dict = global_average_trends(vari="resid",season=season, start = start, stop = stop,noise=rnoise)
    TP = {}
    for k in T_sn_dict.keys():
        mod_d={}
        for m in T_sn_dict["hist"].keys():
            mod_d[m] = [T_sn_dict[k][m],P_sn_dict[k][m], R_sn_dict[k][m]]
        TP[k] = mod_d
    return TP

def TP_dict_decadal(season = "YEAR",increment = 5,trend_length = 10):
    TP = {}
    start = cdtime.comptime(1900,1,1)
    stop = start.add(trend_length,cdtime.Years)
    end = cdtime.comptime(2005,12,31)
    while stop.cmp(end)<=0:
        print start
        key = str(start.year)+"_"+str(stop.year)
        TP_dec = TP_dict(season,start=start,stop=stop)
        TP[key] = TP_dec
        start = start.add(increment,cdtime.Years)
        stop = start.add(trend_length,cdtime.Years)
    f = open("TP_dictionary_trendlength"+str(trend_length)+"increment"+str(increment)+".dat","w")
    pickle.dump(TP,f)
    f.close()
    return TP
    
def plot_global_average_trends(TP,AntNat=False):
    axes = [plt.subplot(2,2,i) for i in np.arange(4)+1]
    models = ["CCSM4","GISS-E2-H","GISS-E2-R","ALL"]
    if AntNat:
        forcings = ["Ant","hist","Nat"]
    else:
        forcings = ["AA","GHG","LU","Oz","Sl","Vl","hist"]
    for k in forcings:
        for i in range(4):
            model = models[i]
            ax = axes[i]
            t=TP[k][model][0]
            p = TP[k][model][1]
            ax.plot(t,p,get_marker(k),color=get_color(k),mec = get_color(k),markersize=10,alpha=.5,label=k)
            plt.legend(loc=0,numpoints=1)
    for ax in axes:
        ax.set_xlim(-15,15)
        ax.set_ylim(-15,15)
        ax.axvline(0,color="k",linestyle=":")
        ax.axhline(0,color="k",linestyle=":")


def get_sum(TP,model,conf=.95,AntNat = False):
    if AntNat:
        hforcings = ["Ant","Nat"]
    elif model.find("*p3")>=0:
         hforcings = ["AA","CH4","CO2","LU","CFCs","Sl","Vl"]
    else:

        hforcings = ["AA","GHG","LU","Oz","Sl","Vl"]
    
    t_mean = sum([np.average(TP[forcing][model][0]) for forcing in hforcings])
    p_mean = sum([np.average(TP[forcing][model][1]) for forcing in hforcings])

    N = [len(TP[forcing][model][0]) for forcing in hforcings]
    
    df = sum(N)-(len(hforcings)+1)
    
    prefac = np.sqrt(np.sum([1/float(N[i]) for i in range(len(N))]))
    t_samps = [TP[fcg][model][0] for fcg in hforcings]
    p_samps = [TP[fcg][model][1] for fcg in hforcings]
    t_std = prefac*pooled_std(t_samps)
    p_std = prefac*pooled_std(p_samps)

    t_conf = np.array(stats.t.interval(conf,df))*t_std
    p_conf = np.array(stats.t.interval(conf,df))*p_std
    
    return t_mean+t_conf,p_mean+p_conf

def get_t_conf(X,conf=.95):
    X_mean  = np.average(X)
    N = len(X)
    df = N-1
    std_err = np.std(X)/np.sqrt(float(N))
    X_conf = np.array(stats.t.interval(conf,df))*std_err
    return X_mean+X_conf

def pretty_plot_TP(TP,AntNat=False,conf=.95,forcings=None,models=None,legend=True,white=False,shade=True,alpha=1.,nrow=2,ncol=2):
   # axes = [plt.subplot(2,2,i) for i in np.arange(4)+1]
    

    if models is None:
    
        models = ["CCSM4","GISS-E2-H","GISS-E2-R","ALL"]
    nmod = len(models)
    axes = [plt.subplot(nrow,ncol,i) for i in np.arange(nmod)+1]
    if forcings is None:
        if AntNat:
            forcings = ["Ant","hist","Nat"]
            
        else:
            forcings = ["AA","GHG","LU","Oz","Sl","Vl","hist","CH4","CO2","CFCs"]
        sumlabel = "Sum of individual forcings"
            
    else:
        sumlabel=''

    for i in range(nmod):
        model = models[i]
        ax = axes[i]
        tnoise = global_average_noise("tas","YEAR")
        pnoise = global_average_noise("pr","YEAR")
        Nt = np.std(get_slopes(tnoise,TP["trend_length"]))
        Np = np.std(get_slopes(pnoise,TP["trend_length"]))
        sens = hydrological_sensitivity(model)
        slp = sens/Np*Nt
        ax.plot(np.arange(-12,13),slp*np.arange(-12,13),"k")
        for k in forcings:
            t=TP[k][model][0]
            if len(t) ==0:
                continue
            p = TP[k][model][1]
            ta = np.average(t)
            pa = np.average(p)
            #ax.plot([ta],[pa],"o",markersize=7,color=get_color(k),mec=get_color(k),zorder=3)
            t_conf = get_t_conf(t,conf=conf)
            p_conf = get_t_conf(p,conf=conf)
            w = t_conf[1]-t_conf[0]
            h = p_conf[1]-p_conf[0]
            ll = (t_conf[0],p_conf[0])
            if white:
                color="w"
                zorder=3
                #alpha=1.
                ec = get_color(k)
            else:
                color = get_color(k)
                ec = "k"
                zorder=3
                #alpha=.5
            R = Rectangle(ll,width=w,height=h,color=color,alpha=alpha,ec=ec,zorder=zorder,label=k,lw=1)
            ax.add_patch(R)
            
        t_conf,p_conf = get_sum(TP,model,conf=conf,AntNat=AntNat)
        t_mean = np.average(t_conf)
        p_mean = np.average(p_conf)
        if sumlabel != "":
            #ax.plot([t_mean],[p_mean],"o",markersize=7,color=get_color("SUM"),mec=get_color("SUM"),zorder=3)

            w = t_conf[1]-t_conf[0]
            h = p_conf[1]-p_conf[0]
            ll = (t_conf[0],p_conf[0])
            if white:
                color="w"
                zorder=2
                ec = get_color("SUM")
                #alpha=1.
            else:
                color = get_color("SUM")
                zorder=2
                ec="k"
                #alpha=.5
            R = Rectangle(ll,width=w,height=h,color=color,ec = ec,alpha=alpha,label=sumlabel,lw=1,zorder=zorder)
            ax.add_patch(R)
        ax.set_title(model)
#    for ax in axes:
        ax.set_xlim(-14,14)
        ax.set_ylim(-14,14)
        ax.axvline(0,color="k",linestyle=":")
        ax.axhline(0,color="k",linestyle=":")
        ax.set_ylabel("P Signal/Noise")
        ax.set_xlabel("T Signal/Noise")
        ax.set_aspect("equal")
        
        
        a,b=stats.norm.interval(conf)
       # Rt = Rectangle((a,a-14),width=b-a,height=28+(b-a),hatch="//",fill=None,zorder=1,alpha=.5)
       # Rp = Rectangle((a-14,a),width=(b-a)+28,height=(b-a),hatch="\\\\",fill=None,zorder=1,alpha=.5)
        Rt = Rectangle((a,a-14),width=b-a,height=28+(b-a),zorder=1,alpha=.5,color=cm.gray(.5))
        Rp = Rectangle((a-14,a),width=(b-a)+28,height=(b-a),color=cm.gray(.5),zorder=1,alpha=.5)
        if shade:
            ax.add_patch(Rt)
            ax.add_patch(Rp)
        if legend:
           # if i == 0:
            leg=ax.legend(loc=2,numpoints=1,ncol=3,fontsize=8)
            leg.set_frame_on(False)
        
def get_TP_stats(sn_dict,model,typ="t",return_interval=None,AntNat=False):
    if AntNat:
        forcings = ["Ant","Nat"]
    else:
        forcings = ["AA","GHG","LU","Oz","Sl","Vl"]

    if lower(typ).find("t")==0:
        idx = 0
    elif lower(typ).find("p")==0:
        idx = 1
    else:
        idx = 2

    the_mean = np.average(sn_dict["hist"][model][idx]) - sum([np.average(sn_dict[forcing][model][idx]) for forcing in forcings])

    N = [len(sn_dict["hist"][model][idx])]+ [len(sn_dict[forcing][model][idx]) for forcing in forcings]
    
    df = sum(N)-(len(forcings)+1)
    
    prefac = np.sqrt(np.sum([1/float(N[i]) for i in range(len(N))]))
    samps = [sn_dict["hist"][model][idx]]+ [sn_dict[fcg][model][idx] for fcg in forcings]
    STD = prefac*pooled_std(samps)
    
    
    tval = the_mean/STD
    pval = stats.t.sf(np.abs(tval),df)*2

    
    
    
    
    if return_interval is None:
        return the_mean,STD,the_mean/STD,pval,df
    else:
        confi = STD*np.array(stats.t.interval(return_interval,df))+the_mean
        
        return the_mean,confi

def time_evolution_diff(TP,confi=.95,cmap=cm.RdYlBu,markersize=10,AntNat=False):
    keys = sorted(TP.keys())
    yrs = np.array([int(x.split("_")[0]) for x in keys])
    ax1 = plt.subplot(311)
    ax2 = plt.subplot(312)
    ax3 = plt.subplot(313)
    i=0
    
    for model in ["CCSM4","GISS-E2-R*p1","GISS-E2-H*p1","GISS-E2-R*p3","GISS-E2-H*p3"]:
        mns = np.array([get_TP_stats(TP[k],model,typ="t",AntNat=AntNat)[0] for k in keys])
        stds = np.array([get_TP_stats(TP[k],model,typ="t",AntNat=AntNat)[1] for k in keys])
        df = get_TP_stats(TP[keys[0]],model,typ='t',AntNat=AntNat)[-1]
        color = cmap((i+.25)/4.25)

        C = stats.t.interval(confi,df)[1]
       # ax1.fill_between(yrs,mns-C*stds,mns+C*stds,color=color,alpha=.4)
        I = np.where(mns>(C*stds))[0]
        J = np.where(mns<(-C*stds))[0]
        ax1.plot(yrs,mns,color=color,linewidth=3,label=model)
        ax1.plot(yrs[I],mns[I],"o",color=cmap(i/4.25),mec=color,markersize=markersize)
        ax1.plot(yrs[J],mns[J],"s",color=cmap(i/4.25),mec=color,markersize=markersize)

        pmns = np.array([get_TP_stats(TP[k],model,typ="p",AntNat=AntNat)[0] for k in keys])
        pstds = np.array([get_TP_stats(TP[k],model,typ="p",AntNat=AntNat)[1] for k in keys])
        df = get_TP_stats(TP[keys[0]],model,typ='p',AntNat=AntNat)[-1]
        C = stats.t.interval(confi,df)[1]
        I = np.where(pmns>(C*pstds))[0]
        J = np.where(pmns<(-C*pstds))[0]
       # ax2.fill_between(yrs,pmns-C*pstds,pmns+C*pstds,color=color,alpha=.4)
        ax2.plot(yrs,pmns,color=color,linewidth=3,label=model)
        ax2.plot(yrs[I],pmns[I],"o",color=cmap(i/4.25),mec=color,markersize=markersize)
        ax2.plot(yrs[J],pmns[J],"s",color=cmap(i/4.25),mec=color,markersize=markersize)

        rmns = np.array([get_TP_stats(TP[k],model,typ="resid",AntNat=AntNat)[0] for k in keys])
        rstds = np.array([get_TP_stats(TP[k],model,typ="resid",AntNat=AntNat)[1] for k in keys])
        df = get_TP_stats(TP[keys[0]],model,typ='resid',AntNat=AntNat)[-1]
        C = stats.t.interval(confi,df)[1]
        I = np.where(rmns>(C*rstds))[0]
        J = np.where(rmns<(-C*rstds))[0]
        #ax3.fill_between(yrs,rmns-C*rstds,rmns+C*rstds,color=color,alpha=.4)
        ax3.plot(yrs,rmns,color=color,linewidth=3,label=model)
        ax3.plot(yrs[I],rmns[I],"o",color=cmap(i/4.25),mec=color,markersize=markersize)
        ax3.plot(yrs[J],rmns[J],"s",color=cmap(i/4.25),mec=color,markersize=markersize)
        
        i+=1
    ax1.legend(loc=0,ncol=5,fontsize=8)
    ax1.set_title("T")
    ax1.axhline(0,color="k",linestyle=":")
    ax1.set_xlim(1900,1975)
    ax1.set_ylim(-4,4)

    
    ax2.set_title("P")
   # ax2.legend(loc=0)
    ax2.axhline(0,color="k",linestyle=":")
    ax2.set_xlim(1900,1975)
    ax2.set_ylim(-4,4)
    
    ax3.set_title("P-H*T")
    #ax3.legend(loc=0)
    ax3.axhline(0,color="k",linestyle=":")
    ax3.set_xlim(1900,1975)
    ax3.set_ylim(-4,4)

def hydrological_sensitivity(model,season="YEAR",in_percent=False):
    direc = "/kate/T/GLOBAL_MEAN/piControl/"
    if model.find("GISS")==0:
        if model.find("*p")>=0:
            model = model.split("*p")[0]
   # print direc
    fname = glob.glob(direc+"*."+model+".*"+season+"*")[0]
    #print fname
    f = cdms.open(fname)
    t = MV.masked_where(np.abs(f("tas"))>1.e10,f("tas"))
    pr = MV.masked_where(np.abs(f("pr"))>1.e10,f("pr"))
    
    sens = float(genutil.statistics.linearregression(pr,x=t,nointercept=True))
    f.close()
    return sens
import copy
def historicalGHG_p3_kludge(TP):
    X = copy.deepcopy(TP)
    if not ("AA" in TP.keys()):
        for k in TP.keys():
            for a in ["LU"]:
                X[k][a]["GISS-E2-H*p3"] =  TP[k][a]["GISS-E2-H*p1"]   
                X[k][a]["GISS-E2-R*p3"] =  TP[k][a]["GISS-E2-R*p1"]
    else:
        for a in ["LU"]:
            X[a]["GISS-E2-H*p3"] =  TP[a]["GISS-E2-H*p1"]   
            X[a]["GISS-E2-R*p3"] =  TP[a]["GISS-E2-R*p1"]
    return X
                                          

def get_resid(fil):
    f=cdms.open(fil)
    t=f("tas")
    p=f("pr")
    model = fil.split(".")[1]
    sens = hydrological_sensitivity(model)
    r = p-sens*t
    return r

def WTF_volcanoes(forcing="Vl"):
    direc = "/kate/T/GLOBAL_MEAN/"+forcing+"/"
    modelh = "GISS-E2-H*"
    modelr = "GISS-E2-R*"
    models = [modelh,modelr]
    cmaps = [cm.Reds(.4),cm.Reds(.9),cm.Blues(.4),cm.Blues(.9)]
    ps = ["*p1","*p3"]
    ax1 = plt.subplot(311)
    ax2 = plt.subplot(312)
    ax3 = plt.subplot(313)
    i = -1
    TAS = []
    for model in models:
        
        for p in ps:
            i+=1
            color = cmaps[i]
            files = glob.glob(direc+"*"+model+p+"*YEAR*")
            tas = MV.zeros((len(files),106))
            pr = MV.zeros((len(files),106))
            resid =  MV.zeros((len(files),106))
            #print files
            c=0
            for fil in files:
                f = cdms.open(fil)
                tas[c] = f("tas")
                pr[c] = f("pr")
                t = get_plottable_time(f("tas"))
                resid[c] = get_resid(fil)
                c+=1
                
            label = fil.split(".")[1]+" "+fil.split(".")[3].split("i1")[1]
            ax1.plot(t,MV.average(tas,axis=0).asma(),color = color,label=label,lw=3)
            TAS+=[tas]
   
            ply = np.polyfit(t,MV.average(tas,axis=0).asma(),1)
            ax1.plot(t,np.polyval(ply,t),color=color,linestyle="--")
                    
            ax2.plot(t,60.*60.*24.*MV.average(pr,axis=0).asma(),color=color,label=label,lw=3)
            ply = np.polyfit(t,60.*60.*24.*MV.average(pr,axis=0).asma(),1)
            ax2.plot(t,np.polyval(ply,t),color=color,linestyle="--")
                    
            ax3.plot(t,60.*60.*24*MV.average(resid,axis=0).asma().anom(),color=color,label=label,lw=3)
            ply = np.polyfit(t,60.*60.*24*MV.average(resid,axis=0).asma().anom(),1)
            ax3.plot(t,np.polyval(ply,t),color=color,linestyle="--")

            [ax.legend(loc=0,ncol=4,fontsize=8) for ax in plt.gcf().axes]
            [ax.set_xlim(1900,2005) for ax in plt.gcf().axes]
            ax1.set_ylabel("T (K)")
            ax2.set_ylabel(r"P (mm day$^{-1}$)")
            ax3.set_ylabel(r"$P - S T$ (mm day$^{-1}$)")

            ax1.set_title("(a) Global-average, annual-average temperature")
            ax2.set_title("(b) Global-average, annual-average precipitation")
            ax3.set_title("(c) Global-average, annual-average precipitation departure from expected T scaling")



def giss_forcing_lookup(forcing,variable,model="H",p=1):
    forcing_search = "*p"+str(p)
    if forcing.find("Nat")>=0:
        path = "/work/cmip5/historicalNat/atm/mo/"+variable+"/"
        
    elif forcing.find("hist")>=0:
        path = "/work/cmip5/historical/atm/mo/"+variable+"/"    
    elif forcing.find("GHG")>=0:
        path = "/work/cmip5/historical/atm/mo/"+variable+"/" 
    else:
      path = "/work/cmip5/historicalMisc/atm/mo/"+variable+"/"

   
    d = {}
    d["Sl"]= '02'
    d["Vl"]= '03'
    d["LU"]= '04'
    d["Oz"]= '05'
    d["AA"] = ['310','107']
    d["CH4"] = '11'
    d["CFCs"] = '12'
    d['CO2'] = '13'
    
    model_search = "GISS-E2-"+model
    if forcing in d.keys():
        if forcing != "AA":
            forcing_search = forcing_search+d[forcing]
        else:
            test = d["AA"]
            ok=test[int(np.where([x[0] == str(p) for x in test])[0])]
            forcing_search = "*p"+str(ok)
    forcing_search = forcing_search+"*"
    return get_ensemble(path,model_search,forcing_search)

def yearly_averages_spatial(fil,start = None,stop=None):
    criteriaarg=[0.99,None] 
    if start is None:
        start = cdtime.comptime(1900,1,1)
    if stop is None:
        stop = cdtime.comptime(2005,12,31)
    f = cdms.open(fil)
    variable = fil.split(".")[7]
    data = f(variable,time=(start,stop))
    cdutil.setTimeBoundsMonthly(data)
    yearly = cdutil.YEAR(data,criteriaarg=criteriaarg)
    f.close()
    return yearly
               
def WTF_spatial(forcing="Vl",variable="tas",model="R"):
    files_nint = giss_forcing_lookup(forcing,variable,model=model,p=1)
    files_int = giss_forcing_lookup(forcing,variable,model=model,p=3)
    f = cdms.open(files_nint[0])
    the_grid = f[variable].getGrid()
    nlat,nlon = the_grid.shape
   
    f.close()

    NINT = MV.zeros((len(files_nint),106,nlat,nlon))
    INT = MV.zeros((len(files_int),106,nlat,nlon))

    for i in range(len(files_nint)):
        fil = files_nint[i]
        print fil
        Yn = yearly_averages_spatial(fil)
        NINT[i] = Yn
    for i in range(len(files_int)):
        fil = files_int[i]
        Yi = yearly_averages_spatial(fil)
        INT[i] = Yi
    for i in range(len(Yi.getAxisList())):
        INT.setAxis(i+1,Yi.getAxis(i))
    for i in range(len(Yn.getAxisList())):
        NINT.setAxis(i+1,Yn.getAxis(i))
        
    return NINT, INT
    
                    
def generate_all_spatial(forcing="Vl"):
    for model in ["R","H"]:
        f = cdms.open("DATA/cmip5.GISS-E2-"+model+"."+forcing+".intvsnint.nc","w")
        for variable in ["tas","pr"]:
            NINT,INT = WTF_spatial(forcing=forcing,variable=variable,model=model)
            NINT.id = variable+"_nint"
            INT.id = variable+"_int"
            f.write(NINT)
            f.write(INT)
        f.close()
    
def bmap(X,**kwargs):
    """ quick plot of data on a lat,lon grid """
   # lon = X.getLongitude()[:]
    #lat = X.getLatitude()[:]
    lon = X.getLongitude().getBounds()[:,0]
    lat = X.getLatitude().getBounds()[:,0]
    m = Basemap(lon_0=np.median(lon))#,projection="moll")
    
        
    x,y=m(*np.meshgrid(lon,lat))
    #if vmin is None:
    m.pcolormesh(x,y,X,**kwargs)
    #else:
     #   m.pcolor(x,y,X,vmin=vmin,vmax=vmax)
    return m

def cmap(X,**kwargs):
    """ quick plot of data on a lat,lon grid """
    lon = X.getLongitude()[:]
    lat = X.getLatitude()[:]
    
    m = Basemap(lon_0=np.median(lon),projection="moll")
    
        
    x,y=m(*np.meshgrid(lon,lat))
    
    m.contourf(x,y,X)
    
    return m
def get_trends(X,p=0.05):
    """ Get trends significantly different from 0 at 1-p confidence level"""
    trends,errors,pt1,pt2,pf1,pf2 = genutil.statistics.linearregression(X,nointercept=1,error=2,probability=1)
    return MV.masked_where(pt2>p,trends)   


def plot_nint_or_int_trends(variable,model,forcing,ax=None,p = 1):
    f = cdms.open("DATA/cmip5.GISS-E2-"+model+"."+forcing+".intvsnint.nc")
    if p == 1:
        data = f(variable+"_nint")
    else:
        data = f(variable+"_int")
    f.close()
    if variable == "pr":
        fac = 60*60*24. #convert to mm/day
    else:
        fac = 1.
    return fac*data
    
    
#    slopes = genutil.stat

def find_ij(data,(xy)):
    lat = data.getLatitude()[:]
    lon = data.getLongitude()[:]
    x=xy[0]
    y=xy[1]
    i = np.argmin(np.abs(lat-x))
    j = np.argmin(np.abs(lon-y))
    return i,j
class InteractiveMap():
    def __init__(self,data,**kwargs):
        self.data = data
        self.avg = MV.average(data,axis=0)
        self.slopes = genutil.statistics.linearregression(self.avg,nointercept=1)*3650.
        a = max([np.abs(np.min(self.slopes)),np.abs(np.max(self.slopes))])
        self.m = bmap(self.slopes,vmin=-1*a,vmax=a)
        self.m.drawcoastlines()
        self.fig = plt.gcf()
        self.lat = data.getLatitude()
        self.lon = data.getLongitude()
        cid = self.fig.canvas.mpl_connect('button_press_event',self.onclick)

    def onclick(self,event):
        xy = (event.ydata,event.xdata)
        t = get_plottable_time(self.data)
        i,j = find_ij(self.data,xy)
        X,Y = self.m(self.lon[j],self.lat[i])
        self.m.plot(X,Y,"y*",markersize=10)
        f2 = plt.figure()
        ax2 = f2.add_subplot(111)
        for mod in range(self.data.shape[0]):
            
            plt.plot(t,self.data[mod,:,i,j].asma(),color=cm.gray(.5))
            plt.plot(t,MV.average(self.data,axis=0)[:,i,j].asma(),color="k")
            plt.title("("+str(self.lat[i])+","+str(self.lon[j])+")")
        
