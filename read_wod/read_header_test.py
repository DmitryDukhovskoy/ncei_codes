# Read T profile
from IPython import get_ipython
get_ipython().magic('reset -sf')

import os
import numpy as np
import matplotlib.pyplot as plt
#import torch
import sys
import pdb
#import netCDF4
import importlib
#from netCDF4 import Dataset as ncFile
import timeit
import binascii
import struct

#sys.path.append('/data/home004/dmitry.dukhovskoy/python_codes/MyPython')
#import mod_utils_fig
#from mod_utils_fig import bottom_text
#from mod_utils_fig import plot_fld2D

prb = 'XBT' # Probe name
pthdata = '/data/ncei2/OCL/data/'+prb+'/'

fhdr = pthdata+'headers'

# Check system endian
print('Endian sysetm: '+sys.byteorder)

#
# Format for main header information:
# CC code - 2 char
# lat, lon - 32-bit float little endian
# year, mo, day - 32-bit int
# GMT time - 32-bit float
# Cruise ID - int
# obs - int, # of observed levels in profile
# nstd -  # of st levels in profile - int
# pointers = depth, T obs levels, T st level (in XBT)
fmt = '2sffiiifiiiiii'
rec_size = struct.calcsize(fmt)
fid = open(fhdr,'rb')
try: 
  fid.seek(0)
  bd = fid.read(rec_size)
  tmp = struct.unpack(fmt,bd)

  CC  = tmp[0].decode()  # convert binary --> string
  lat = tmp[1]
  lon = tmp[2]
  yr  = tmp[3]
  mo  = tmp[4]
  md  = tmp[5]
  gtm = tmp[6]
  cid = tmp[7]
  nobs= tmp[8]
  nstd= tmp[9]
  p1  = tmp[10]
  p2  = tmp[11] 
  p3 =  tmp[12]

finally:
  fid.close()


