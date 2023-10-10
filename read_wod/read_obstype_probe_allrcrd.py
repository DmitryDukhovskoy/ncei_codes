# Loop through all records in file
# Read specified type of observation from a probe
# 
# Imporant note on pointers for 1 variable Probe:
# One variable probe:
# 1 pointer: depth/variable_one_observed
# 1 pointer: variable_one_standard
# 1 pointer: second_header
#
# multi-variable probe:
# 1 pointer: depth
# N pointers: variables_[1 to N]_observed
# N pointers: variabels_{1 to N]_standard
# 1 pointer: second_header
#
# Standard level depths are listed in OCL/sys.inf/standard_depths_{LEVS].dat 
# where LEVS is either
# orig (main database) or full (databases used for climatologies - w21).
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
#import binascii
import struct
sys.path.append('/data/home004/dmitry.dukhovskoy/python_codes/MyPython')
import mod_utils_fig
from mod_utils_fig import bottom_text

#===========================================================
# Choose probe type, variable type (if supported in the probe)
# record # 
#
#===========================================================
#prbnm = 'CTD'         # Probe type
#prbnm = 'DRB'         # Probe type
#prbnm = 'PFL'
#prbnm = 'OSD'
#prbnm = 'GLD'
prbnm = 'XBT'
#prbnm = 'SUR'
#prbnm = 'UOR'
#prbnm = 'MRB'
varnm = 'Temperature'  # read obs type
#varnm = 'Salinity'  # read obs type
#krcrd = 349420  # sequential recrod# to read = 1, ... 
#===========================================================


import mod_utils_wod
importlib.reload(mod_utils_wod)
#reload_mod('mod_utils_wod')
from mod_utils_wod import get_obsparam_name, probe_obs
from mod_utils_wod import read_header, read_obsdepths, standard_depths
from mod_utils_wod import find_var_indx, get_uniqnmb_cast
from mod_utils_wod import read_varobs, read_varstd
LCd_obs, LNm_obs = probe_obs(prbnm) # list of probe obs codes and names
nvar = len(LCd_obs)   # # of observed variables

krcrd = 0
f_err = 0
while f_err == 0:
  krcrd += 1
  print('Reading {0} seq.record={1}'.format(prbnm,krcrd))
  PRBINFO, PRB = read_header(prbnm,krcrd)
  print('Info for extracted cast, probe {0}'.format(prbnm))
  print(PRB.__dict__)

  # Derive depths arrays
  Zobs = read_obsdepths(PRBINFO,PRB)        # observed depth levels
  nstd = PRB.nstd
  Zstd = standard_depths(nzlev=nstd)    # stand.depths, get only needed 
  if Zstd.size == 0:
    Zstd = 0

  if np.min(Zobs) >= 0.:
    Zobs=-Zobs
  if np.min(Zstd) >= 0.:
    Zstd=-Zstd


  # Read observed variable
  Vobs = read_varobs(PRBINFO,PRB,varnm)    # obs variab. at obs. depths
#  Vstd = read_varstd(PRBINFO,PRB,varnm)    # obs. variab. at stand. depths
  if np.size(Vobs) == 0:
#    raise NameError(varnm+' pointer <0: no observation ')
    print('WARNING: {0} pointer <0: no observations'.format(varnm))
    continue

  # Get unique cast number:
  uniqnmb = get_uniqnmb_cast(prbnm,krcrd)

  #
  # Some missing values (?) in stand-depth variables -1e+10
  #Vstd[np.where(np.abs(Vstd)>1.e4)]=np.nan

btx = 'read_obstype_probe.py'
f_plt=0
f_std=0  # plot standard depth ontop
if f_plt == 1:
  plt.ion()
  fig1 = plt.figure(1,figsize=(6,8), constrained_layout=False)
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.5, 0.8])

  plt.plot(Vobs,Zobs,'.-')

  cll = '{0} {1} obs depths, rec={2}'.format(prbnm,varnm,krcrd)
  if f_std==1:
    plt.plot(Vstd,Zstd,'.-')
    cll = '{0} {1} obs & std depths, rec={2}'.format(prbnm,varnm,krcrd)

  plt.grid()
  ax1.set_title(cll)

  CC  = PRB.CC
  lat = PRB.lat
  lon = PRB.lon
  yr  = PRB.yr
  mo  = PRB.mo
  md  = PRB.md
  gtm = PRB.gtm
  cid = PRB.cid

  Txt = ['Cruise# {1}{0}'.format(cid,CC)]
  Txt.append('lat={0:.2f}N, lon={0:.2f}E'.format(lat,lon))
  Txt.append('{0:02d}/{1:02d}/{2}:{3:.2f}'.format(md,mo,yr,gtm))
  Txt.append('UnqNmb={0}'.format(uniqnmb))
  Txt.append('Seq.Recrd={0}'.format(krcrd))
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

  bottom_text(btx, pos=[0.0, 0.05], fsz=9)

  
  
    


