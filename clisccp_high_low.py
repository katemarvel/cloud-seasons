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

def low_mid_high(fname):
    fgrid = cdms.open("~/precip.mon.mean.nc")
    the_grid = fgrid["precip"].getGrid()
    f = cdms.open(fname)
    cl = f("clisccp")
    #Sum over optical depth bins
    tau_ax=cl.getAxisIds().index("tau")
    cl_tau = MV.sum(cl,axis=tau_ax)
    #Now sum over height to get low,mid,high
    levax=cl_tau.getAxisIds().index('plev')
    low = MV.sum(cl_tau(level=(1000*100,680*100)),axis=levax)
    mid = MV.sum(cl_tau(level=(680*100,440*100)),axis=levax)
    high = MV.sum(cl_tau(level=(440*100,50*100)),axis=levax)

    low = low.regrid(the_grid,regridTool='regrid2')
    low.id = "low_clt"
    mid = mid.regrid(the_grid,regridTool='regrid2')
    mid.id = "mid_clt"
    high = high.regrid(the_grid,regridTool='regrid2')
    high.id = "high_clt"
    fgrid.close()
    f.close()
    return low,mid,high


def write_all_low_mid_high(experiment):    
    sorted_files = sorted(cmip5.get_datafiles(experiment,"clisccp"))
    for fname in sorted_files[0:1]:
        try:
            low,mid,high = low_mid_high(fname)
            trunc_fname = fname.split("/")[-1]
            
            os.system("mkdir "+ "/kate/CLISCCP/LOW/"+experiment+"/")
            writename_low = "/kate/CLISCCP/LOW/"+experiment+"/"+trunc_fname.replace("clisccp","low").replace(".xml",".nc")
            flow = cdms.open(writename_low,"w")
            flow.write(low)
            
            os.system("mkdir "+ "/kate/CLISCCP/LOW/"+experiment+"/")
            writename_mid = "/kate/CLISCCP/MID/"+experiment+"/"+trunc_fname.replace("clisccp","mid").replace(".xml",".nc")
            fmid = cdms.open(writename_mid,"w")
            fmid.write(mid)

            os.system("mkdir "+ "/kate/CLISCCP/LOW/"+experiment+"/")
            writename_high = "/kate/CLISCCP/HIGH/"+experiment+"/"+trunc_fname.replace("clisccp","high").replace(".xml",".nc")
            fhigh = cdms.open(writename_high,"w")
            fhigh.write(high)
        except:
            print "*****************"
            print fname+" IS BAD"
            print "*****************"
            

if __name__ == "__main__":
    experiments = ["piControl","historical","1pctCO2","amip"]
    for experiment in experiments:
        write_all_low_mid_high(experiment)
