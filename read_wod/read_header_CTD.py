# Read CTD header file

import os
import numpy as np
import sys
import importlib
import binascii
import struct

#sys.path.append('/data/home004/dmitry.dukhovskoy/python_codes/MyPython')
#import mod_utils_fig
#from mod_utils_fig import bottom_text
#from mod_utils_fig import plot_fld2D

prb = 'CTD' # Probe name
krcrd = 150
pthdata = '/data/ncei2/OCL/data/'+prb+'/'

fhdr = pthdata+'headers'

# Check system endian
print('Endian sysetm: '+sys.byteorder)

import mod_utils_wod
importlib.reload(mod_utils_wod)
from mod_utils_wod import get_obsparam, get_obsparam_code
from mod_utils_wod import get_obsparam_name, probe_obs
#PARAM = get_obsparam()
#vcode = get_obsparam_code("Temperature")
#vname = get_obsparam_name(25)
LCd_obs, LNm_obs = probe_obs(prb) # list of probe obs codes and names
nvar = len(LCd_obs)   # # of observed variables

class ProbePntr:
  def __init__(self,tmp,nvar,LNm_obs):
    self.CC     = tmp[0].decode()  # convert binary --> string
    self.lat    = tmp[1]
    self.lon    = tmp[2]
    self.yr     = tmp[3]
    self.mo     = tmp[4]
    self.md     = tmp[5]
    self.gtm    = tmp[6]
    self.cid    = tmp[7]      # Cruise ID
    self.nobs   = tmp[8]
    self.nstd   = tmp[9]
    self.pdepth = tmp[10]   # pointer - depth
    self.VarNm  = LNm_obs
    p1          = 11
    p2          = p1+nvar-1
    self.pobs   = tmp[p1:p2]
    p1          = p2+1
    p2          = p1+nvar-1
    self.pstd   = tmp[p1:p2]
#
# Format for main header information:
# CC code - 2 char
# lat, lon - 32-bit float little endian
# year, mo, day - 32-bit int
# GMT time - 32-bit float
# Cruise ID - int
# obs - int, # of observed levels in profile
# nstd -  # of st levels in profile - int
# pointers = depth, var1 obs levels, var1 st level, ... 
fmt = '2sffiiifiiii'
for ii in range(nvar):
  fmt = fmt+'ii'

rec_size = struct.calcsize(fmt)
fid = open(fhdr,'rb')
try: 
  fid.seek(0)
  bd = fid.read(rec_size)
  tmp = struct.unpack(fmt,bd)

  PRB = ProbePntr(tmp,nvar,LNm_obs)

finally:
  fid.close()


