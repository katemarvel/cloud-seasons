import glob
import sys
import cdms2 as cdms
import numpy as np
import MV2 as MV
import difflib
import scipy.stats as stats
global crunchy
import socket
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

from scipy import signal

def FourierPlot(tas):
    """Plot the Fourier power spectrum as a function of Fourier period (1/frequency)"""
    detrend = signal.detrend(tas)
    L = len(tas)
    freqs = np.fft.fftfreq(L)
    tas_fft = np.fft.fft(detrend)
    R = tas_fft.real
    Im = tas_fft.imag
    mag = np.sqrt(R**2+Im**2)
    plt.plot(1/freqs,mag)
def annual_cycle_dominant(tas):
    """Check to see whether the annual cycle is dominant"""
    detrend = signal.detrend(tas)
    L = len(tas)
    freqs = np.fft.fftfreq(L)
    tas_fft = np.fft.fft(detrend)
    R = tas_fft.real
    Im = tas_fft.imag
    mag = np.sqrt(R**2+Im**2)
    the_period = 1./np.abs(freqs[np.argmax(mag)])
    return the_period

def get_dominant_cycle(tas):
    """For a 2D variable, calculate the period of the dominant Fourier mode everywhere """
    nt,nlat,nlon = tas.shape
    to_mask = MV.zeros((nlat,nlon))
    for i in range(nlat):
        for j in range(nlon):
            to_mask[i,j]=annual_cycle_dominant(tas[:,i,j])
    to_mask.setAxisList(tas.getAxisList()[1:])
    return to_mask
def mask_cycle_subdominant(tas,period = 12):
    #find the closest Fourier frequency to 1/12 months
    L = len(tas)
    freqs = np.fft.fftfreq(L)
    closest = np.abs(freqs-1./period)
    i = np.argmin(closest)
    cutoff = 1/freqs[i]
    to_mask = get_dominant_cycle(tas)
    return to_mask != cutoff


def get_cycle(tas,period=12,return_complex=False):
    """Return the Fourier magnitude and phase for a given period (default 12 months)"""
    L = len(tas)
    freqs = np.fft.fftfreq(L)
    closest = np.abs(freqs-1./period)
#    i = np.where(freqs == 1./period)[0]
    i = np.argmin(closest)
    #print 1/freqs[i]
    tas_fft = np.fft.fft(tas)/L
    R = tas_fft.real
    Im = tas_fft.imag
    if return_complex:
        return R[i],Im[i]
    else:
        mag = np.sqrt(R**2+Im**2)
        phase = np.arctan2(Im,R)
        return mag[i],phase[i]

def get_tan_phase(tas,period=12):
    """Return the tangent of the phase associated with Fourier mode"""
    L = len(tas)
    freqs = np.fft.fftfreq(L)
    closest = np.abs(freqs-1./period)
#    i = np.where(freqs == 1./period)[0]
    i = np.argmin(closest)
    #print 1/freqs[i]
    tas_fft = np.fft.fft(tas)/L
    R = tas_fft.real
    Im = tas_fft.imag
    mag = np.sqrt(R**2+Im**2)
    phase = Im/R
    return phase[i]

def get_semiannual_cycle(tas):
    """Helper function: get the magnitude and phase for the mode with period 6 """
    return get_cycle(tas,period=6)

def get_cycle_map(tas,period = 12):
    ntime,nlat,nlon = tas.shape
    AMP = MV.zeros((nlat,nlon))
    PHI = MV.zeros((nlat,nlon))
    for i in range(nlat):
        for j in range(nlon):
            mag,phase = get_cycle(tas[:,i,j],period=period)
            AMP[i,j] = mag
#            PHI[i,j]= phase_to_day(phase)
            PHI[i,j] = phase
    AMP.setAxis(0,tas.getLatitude())
    AMP.setAxis(1,tas.getLongitude())
    PHI.setAxis(0,tas.getLatitude())
    PHI.setAxis(1,tas.getLongitude())
    return AMP,PHI

def get_zonal_cycle(tas,period=12.):
    if "longitude" in tas.getAxisIds():
        tas = cdutil.averager(tas,axis='x')
    AMP = MV.zeros(tas.shape[1])
    PHI = MV.zeros(tas.shape[1])
    for j in range(tas.shape[1]):
        mag,phase = get_cycle(tas[:,j],period=period)
        AMP[j] = mag
        PHI[j] = phase
        
    AMP.setAxis(0,tas.getLatitude())
    PHI.setAxis(0,tas.getLatitude())
    return AMP,PHI




def get_annual_cycle_trends(X):
    """Given a variable X, get amplitude gain and phase lag time series """
    nt,nlat,nlon=X.shape
    time = X.getTime()[::12]
    tax = cdms.createAxis(time)
    tax.designateTime()
    atts = X.getTime().attributes
    for k in atts:
        setattr(tax,k,atts[k])
    axlist = [tax,X.getLatitude(),X.getLongitude()]
    nyears = nt/12
    yrs = X.reshape((nyears,12,nlat,nlon))
    for i in [2,3]:
        yrs.setAxis(i,X.getAxis(i-1))
    AMP = MV.zeros((nyears,nlat,nlon))
    PHI = MV.zeros((nyears,nlat,nlon))

    for y in range(nyears):
        amp,phi = get_cycle_map(yrs[y])
        AMP[y] = amp
        PHI[y] = phi
    amp_solar,phase_solar = get_insolation()
    AMP = AMP/amp_solar
    AMP.setAxisList(axlist)
    
    PHI.setAxisList(axlist)
    PHI = correct_phase(PHI)
    return AMP,PHI


    
def phase_to_day(phase):
    """Convert phase to day of the year """
    if phase < 0:
        phase += 2*np.pi
    return phase*(365./(2*np.pi))
     
 
def subtract_months(P, reference):
    """ Write P in terms of days relative to solar insolation phase.  Return "forward" phases (leads solar insolation) and "backward" (lags solar insolation)"""
    ref = phase_to_day(reference)
    phase = phase_to_day(P)
   

    if ref < phase:
        backward = phase - ref #Ref lags phase
        forward = -ref - (365 - phase) #Move into the next year
       
    else:
        forward = phase + 365 - ref
        backward = phase - ref


    return forward, backward
        
        
import sys
sys.path.append("")
from Helper import cdms_clone,get_plottable_time,get_orientation,get_slopes



def get_extremum(P,func = np.argmin):
    return np.unravel_index(func(P),P.shape)
    
def merge(forward,backward):
    X = np.ma.zeros(forward.shape)
    mask = np.abs(forward) > np.abs(backward)
    
    I = np.where(mask.flatten())[0]
    J = np.where(~mask.flatten())[0]
    Xf = X.flatten()
    Xf[I] = np.array(backward).flatten()[I]
    Xf[J] = np.array(forward).flatten()[J]
    
    Xt = Xf.reshape(X.shape)
    if 1:
        return Xt
    Diff = np.diff(Xf.reshape(X.shape),axis=0)
    bad = np.where(np.abs(Diff)>365./2)
    for ibad in range(len(bad[0])):
        i= bad[0][::2][ibad]+1
        j = bad[1][::2][ibad]
        k = bad[2][::2][ibad]
        if Xt[i,j,k] == forward[i,j,k]:
            Xt[i,j,k] = backward[i,j,k]
        else:
            Xt[i,j,k] = forward[i,j,k]

    
    return Xt

def check_for_badness(t):
    Diff = np.diff(t)
    if len(np.where(Diff>365/2.)[0]) >0:
        return True
    else:
        return False

def where_bad(P):
    nt,nlat,nlon = P.shape
    X = MV.zeros((nlat,nlon))
    for i in range(nlat):
        for j in range(nlon):
            X[i,j] = check_for_badness(P[:,i,j])
    return X

def fix_bad(t,debug = False,discont=365./2.):
    bad = t.copy()
    too_big = []
    too_small = []
    smallest = 0
    for i in range(len(t)):
        if (bad[i] - bad[smallest]) >= discont:
            if np.abs(365 - bad[i]) <= 365: #Is it an allowed value?
                if np.abs( (bad[i]-365.) - bad[smallest]) < np.abs(bad[i]-bad[smallest]): #Does it make things better?
                    bad[i] = bad[i] - 365.
                    too_big += [i]
        elif (bad[i] - bad[smallest]) <= -discont:
            if np.abs(bad[i]+ 365) <= 365: #Is it an allowed value?
                if np.abs( (bad[i]+365.) - bad[smallest]) < np.abs(bad[i]-bad[smallest]): #Does it make things better?
                    bad[i] = 365 + bad[i]
                    too_small += [i]

    #Need to ensure that no points in the time series are more than discont away from each other
    if debug:
        return bad ,np.array(too_big),np.array(too_small)
    else:
        return bad

def fix_all_bad(P):
    Fix = P.copy()
    nt,nlat,nlon = P.shape
    for j in range(nlat):
        for k in range(nlon):
            time_series = Fix[:,j,k]
            if check_for_badness(time_series):
                Fix[:,j,k] = fix_bad(time_series)
        
    return Fix

def correct_phase(P,reference = None):
    if reference is None:
        amp_solar, phase_solar = get_insolation()
        if phase_solar.shape != P.shape:
            grid = P.getGrid()
            phase_solar = phase_solar.regrid(grid,regridTool='regrid2')
        reference = phase_solar
    
    #print "got insolation"
    Convert2Day = np.vectorize(subtract_months)
    #print "vectorized"
    forward,backward = Convert2Day(P,reference)
    #print "converted to days"
    Merged = merge(forward,backward)
    #print "Merged"
    Pnew = fix_all_bad(Merged)
    #print "fixed discontinuities"
    Pnew = MV.array(Pnew)
    Pnew.setAxisList(P.getAxisList())
    return Pnew



def phase_angle_form(obs,period=12):
    mag,phase = get_cycle(obs,period = period)
    ###TEST
    mag=2*mag
    return mag*np.cos(2*np.pi/period*np.arange(len(obs))+phase)

def var_expl_by_annual_cycle(obs,period = 12,detrend = True):
   # mag,phase = get_cycle(obs)
    if detrend:
        obs = signal.detrend(obs)
    recon = phase_angle_form(obs,period=period)
    return np.corrcoef(recon,obs)[0,1]**2

def variance_map(X,period = 12,detrend = True):
    if len(X.shape )==3:
        nt,nlat,nlon = X.shape
        V = MV.zeros((nlat,nlon))
        for i in range(nlat):
            for j in range(nlon):
                V[i,j]=var_expl_by_annual_cycle(X[:,i,j],period = period,detrend=detrend)
    elif len(X.shape)==2:
        nt,nlat = X.shape
        V = MV.zeros((nlat))
        for i in range(nlat):
            V[i]=var_expl_by_annual_cycle(X[:,i],period = period,detrend=detrend)
    V.setAxisList(X.getAxisList()[1:])
    return V
                  

##### Vectorizing stuff #####
def broadcast(fvec):
    def inner(vec, *args, **kwargs):
        if len(vec.shape) > 1:
            return np.array([inner(row, *args, **kwargs) for row in vec])
        else:
            return fvec(vec, *args, **kwargs)
    return inner

def fast_annual_cycle(X,debug=False,semiann=False,zonal_average=False):
    if len(X.shape)==4:
        nt,nlat,nlon,nmod = X.shape
        has_models=True
    elif len(X.shape)==3:
        nt,nlat,nlon = X.shape
        has_models = False
    elif len(X.shape)==2:
        nt,nlat = X.shape
        has_models=False
    if zonal_average:
        if 'lon' in X.getAxisIds():
            X = cdutil.averager(X,axis='x')
    nyears = nt/12

    newshape = (nyears, 12) + X.shape[1:]
    yrs = X.reshape(newshape)

    if semiann:
        vec_cycle=broadcast(get_semiannual_cycle)
    else:
        vec_cycle = broadcast(get_cycle)
   # print "vectorized"
    apply_everywhere = np.apply_along_axis(vec_cycle,1,yrs)
    
    R = MV.array(apply_everywhere[:,0])
    P = MV.array(apply_everywhere[:,1])
    #print "got R and P"
    
    if debug:
        return R,P
    time = X.getTime()[::12]
    tax = cdms.createAxis(time)
    tax.designateTime()
    tax.id = "time"
    atts = X.getTime().attributes
    for k in atts:
        setattr(tax,k,atts[k])
    print "got new time"
    axlist = [tax]+X.getAxisList()[1:]
    R.setAxisList(axlist)
    P.setAxisList(axlist)
    
    return R,P
    

def decade_fast_annual_cycle(X,debug=False,semiann=False,zonal_average=False,return_Pdays=True):
    if len(X.shape)==4:
        nt,nlat,nlon,nmod = X.shape
        has_models=True
    elif len(X.shape)==3:
        nt,nlat,nlon = X.shape
        has_models = False
    if zonal_average:
        X = cdutil.averager(X,axis='x')
    nyears = nt/60

    newshape = (nyears, 60) + X.shape[1:]
    print newshape
    yrs = X.reshape(newshape)

    if semiann:
        vec_cycle=broadcast(get_semiannual_cycle)
    else:
        vec_cycle = broadcast(get_cycle)
   # print "vectorized"
    apply_everywhere = np.apply_along_axis(vec_cycle,1,yrs)
    
    R = MV.array(apply_everywhere[:,0])
    P = MV.array(apply_everywhere[:,1])
    #print "got R and P"
    
    if debug:
        return R,P
    time = X.getTime()[::60]
    tax = cdms.createAxis(time)
    tax.designateTime()
    tax.id = "time"
    atts = X.getTime().attributes
    for k in atts:
        setattr(tax,k,atts[k])
    print "got new time"
    axlist = [tax]+X.getAxisList()[1:]
    R.setAxisList(axlist)
    P.setAxisList(axlist)
    if return_Pdays is False:
        return R,P
    #Pnew = P.copy()
    #Pold = P.copy()
    Pnew = MV.zeros(P.shape)
    if has_models:
        
        for mod_i in range(nmod):
            chunk = P[:,:,:,mod_i]
            Pnew[:,:,:,mod_i] = correct_phase(chunk)
            
            
    else:
        chunk = P[:,:,:]
        Pnew = correct_phase(chunk)
    Pnew = MV.array(Pnew)
    Pnew.setAxisList(axlist)            
   
    return R,P,Pnew


def mask_data(data,basicmask):
    dim = len(data.shape)
    if dim == 2:
        mask= basicmask
    elif dim ==3:
        nt = data.shape[0]
        mask= np.repeat(basicmask[np.newaxis,:,:],nt,axis=0)
    elif dim ==4:
        nmod,nt,nx,ny = data.shape
        mask= np.repeat(np.repeat(basicmask[np.newaxis,:,:],nt,axis=0)[np.newaxis],nmod,axis=0)
    return MV.masked_where(mask,data)

def get_variance_maps_models(variable="pr",models=None,cmip_dir = None,period=12):
    """ find latitudes in each model where the annual cycle is not dominant"""
    if models is None:
        f = cdms.open("/work/marvel1/SEASONAL/MMA/cmip5.ZONALMMA.historical-rcp85.rip.mo.atm.Amon.pr.ver-1.AmpPhase.nc")
        phase = f("phase")
        models = eval(phase.getAxis(0).models)
        f.close()
    if cmip_dir is None:
        cmip_dir = "/work/cmip5/historical-rcp85/atm/mo/"+variable+"/"
    fobs = cdms.open("/work/marvel1/SEASONAL/OBS/GPCP.precip.mon.mean.nc")
    the_grid = fobs("precip").getGrid()
    nlat,nlon=the_grid.shape
    fobs.close()
    VarianceMaps = np.zeros((len(models),nlat))+1.e20
    counter=0
    for model in models:
        print model
        try:
            fname = sorted(glob.glob(cmip_dir+"*."+model+".*"))[0]
            fp = cdms.open(fname)
            prtest = fp(variable,time=("1979-1-1","2014-12-31")).regrid(the_grid,regridTool='regrid2')
            przonal = cdutil.averager(prtest,axis='x')
            dom = variance_map(przonal,period=period)
            VarianceMaps[counter]=dom
            fp.close()
            counter+=1
        except:
            continue
        
    modax = cdms.createAxis(range(len(models)))
    modax.id = "model"
    modax.models = str(models)
    VarianceMaps = MV.array(VarianceMaps)
    VarianceMaps.setAxis(0,modax)
    VarianceMaps.setAxis(1,the_grid.getLatitude())
    return MV.masked_where(VarianceMaps>1.e10,VarianceMaps)

