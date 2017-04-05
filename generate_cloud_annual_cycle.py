#!/usr/local/anaconda2/envs/latest/bin/python
import glob
import os
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
import seasonal_cycle_utils as sc
import CMIP5_tools as cmip5

from Plotting import *


def write_amplitude_phase(experiment,variable,search_string = "*",overwrite=False):
    path = "/work/cmip5/"+experiment+"/atm/mo/"+variable+"/"
    writepath = "/kate/CLT_ANNUALCYCLE/"+experiment+"/"
    os.system("mkdir "+writepath)
    #Files in the ensemble
    files = cmip5.only_most_recent(glob.glob(path+search_string))
    bad = open("failed_files.dat","w")
    for fname in files:
        f = cdms.open(fname)
        X = f(variable)
        #Need to ensure that X contains full annual cycle (HADLEY CENTER MODELS START IN DECEMBER AND ARE THE WORST BOOOOO HADLEY CENTER)
        if cmip5.start_time(X).month != 1: #if it doesn't start in January, make it start in January
            starttime = cdtime.comptime(cmip5.start_time(X).year+1,1,1)
        else:
            starttime = cmip5.start_time(X)
        if cmip5.stop_time(X).month != 12: #if it doesn't stop in December, make it.
            stoptime = cdtime.comptime(cmip5.stop_time(X).year-1,12,31)
        else:
            stoptime = cmip5.stop_time(X)
        X = X(time=(starttime,stoptime))
            
        fshort = fname.split("/")[-1]
        fwrite = fshort.replace("xml","nc")
        fwrite = fwrite.replace(variable,variable+"AmpPhase")
        if overwrite is False:
            if fwrite in os.listdir(writepath):
                continue
        writefile = cdms.open(writepath+fwrite,"w")
        try:
            R,P = sc.fast_annual_cycle(X)
        except:
            bad.write(fname+"\n")
            
            continue
        historical = experiment.find("historical")>=0
        if historical:
            relative_to = "1996-2009"
        else:
            relative_to = str(cmip5.start_time(X).year)+"-"+str(cmip5.stop_time(X).year)
        Pa = sc.get_phase_anomalies(P,historical=historical)
        if historical:
            Rclim = MV.average(R(time=('1996-1-1','2009-12-31')),axis=0) #this anomaly period chosen to overlap ISCCP/PATMOS most recent time period with all data present
        else:
            Rclim = MV.average(R,axis=0)
        Ra = R-Rclim
        

        R.id = "amp"
        if hasattr(X,"umits"):
            R.units = X.units
        R.long_name = "Amplitude of the "+variable+" annual cycle"
        writefile.write(R)

        P.id = "phase"
        P.units = "radians"
        P.long_name = "Phase of the "+variable+" annual cycle"
        writefile.write(P)
        
        Ra.id = "amp_anom"
        if hasattr(X,"umits"):
            Ra.units = X.units
        Ra.long_name = "Amplitude anomaly relative to "+relative_to
        writefile.write(Ra)
        
        Pa.id = "phase_anom"
        Pa.units = "days"
        Pa.long_name = "Phase anomaly relative to "+relative_to
        writefile.write(Pa)
        
        f.close()
        writefile.close()
    bad.close()
    

    
    
if __name__ == "__main__":
    experiments = ["historical","rcp85","1pctCO2","piControl"]
    variable = "clt"
    for experiment in experiments:
        
        write_amplitude_phase(experiment,variable)
