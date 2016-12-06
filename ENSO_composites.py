import glob
import sys
import cdms2 as cdms
import numpy as np
import MV2 as MV
import difflib
import pickle

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

import CMIP5_tools as cmip5
### Set classic Netcdf (ver 3)
cdms.setNetcdfShuffleFlag(0)
cdms.setNetcdfDeflateFlag(0)
cdms.setNetcdfDeflateLevelFlag(0)

fnew = cdms.open("cltcontrol.nc")
S = fnew("clt")
models = eval(S.getAxis(0).models)
fnew.close()
L = len(models)
SSTS = MV.zeros((L,250*12))+1.e20

NINO34 = cdutil.region.domain(latitude=(-5.,5.),longitude=(190.,240.))
i = 0
for model in models:
    candidates = glob.glob(model.split("ver")[0].replace("clt","tas")+"*")
    tasmodel = cmip5.get_latest_version(candidates)
    try:
        f = cdms.open(tasmodel)
    
        ssts=cdutil.averager(f("tas",NINO34),axis='xy')[:12*250]
        SSTS[i] = ssts
    except:
        print "PROBLEM WITH "+model
        continue
    i+=1
SSTS = MV.masked_where(SSTS>1.e10,SSTS)
SSTS.id = 'sst'
SSTS.setAxis(1,ssts.getTime())
SSTS.setAxis(0,S.getAxis(0))

