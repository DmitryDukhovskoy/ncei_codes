# Read XBT 
# 
# Imporant note on pointers for 1 variable Probe:
# if there is only one variable, there is only one observed level pointer 
# which is used for both depth and the variable.   So, instead of 
# pointer(depth), pointer(variable_one_observed), pointer(variable_one_standard),
#  you have pointer(depth/variable_one_observed), pointer(variable_one_standard).
#
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


sys.path.append('/data/home004/dmitry.dukhovskoy/python_codes/MyPython')
import mod_utils_fig
from mod_utils_fig import bottom_text

#
# Reloading modules:
def reload_mod(__name__):
  try:
      del sys.modules[__name__]
  except AttributeError:
      pass
  import __name__

import mod_utils_wod
importlib.reload(mod_utils_wod)
#reload_mod('mod_utils_wod')
from mod_utils_wod import get_probe_info, read_header_1var, read_obsdepths_1var
from mod_utils_wod import read_varobs_1var, read_varstd_1var

prb   = 'XBT' # Probe name
krcrd = 1434200     # record to read
varnm = 'Temperature' # for XBT only 1 type of var

PRBINFO = get_probe_info(prb)
HEAD  = read_header_1var(PRBINFO,krcrd)

Zobs = read_obsdepths_1var(PRBINFO,HEAD)        # observed depth levels
if np.min(Zobs) >= 0.:
  Zobs=-Zobs

Tobs = read_varobs_1var(PRBINFO,HEAD,varnm)    # read observation
Tstd = read_varstd_1var(PRBINFO,HEAD,varnm)

btx = 'read_XBT.py'
f_plt=1
if f_plt == 1:
  plt.ion()
  fig1 = plt.figure(1,figsize=(6,8), constrained_layout=False)
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.5, 0.8])

  plt.plot(Tobs,Zobs,'.-')

  cll = '{0} Var1_obs_lev, rec={1}'.format(prb,krcrd)
  ax1.set_title(cll)

  CC  = HEAD.get("CC")
  lat = HEAD.get("lat")
  lon = HEAD.get("lon")
  yr  = HEAD.get("yr")
  mo  = HEAD.get("mo")
  md  = HEAD.get("md")
  gtm = HEAD.get("gtm")
  cid = HEAD.get("cid")

  Txt = ['Cruise# {1}{0}'.format(cid,CC)]
  Txt.append('lat={0:.2f}N'.format(lat))
  Txt.append('lon={0:.2f}E'.format(lon))
  Txt.append('{0:02d}/{1:02d}/{2}:{3:.2f}'.format(md,mo,yr,gtm))
  Txt.append('Record={0}'.format(krcrd))
  nt=len(Txt)

  x0 = 1
  y0 = 1
  ax2 = plt.axes([0.7,0.6,0.25,0.1])
  ax2.set(xlim=(x0,4),ylim=(y0,y0+nt))
  for itt in range(nt):
    ssl=Txt[itt]
    x1 = x0
    y1 = y0+itt
    ax2.text(x1,y1,ssl)
  ax2.axis('off')
#  ax2.axes.xaxis.set_visible(False)
#  ax2.axes.yaxis.set_visible(False)

  bottom_text(btx)

  
  
    


