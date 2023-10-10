"""
  Plot statistics
"""
import os
import numpy as np
import matplotlib.pyplot as plt
#import torch
import sys
sys.path.append('/data/home004/dmitry.dukhovskoy/python_codes/MyPython')
import mod_utils_fig
from mod_utils_fig import bottom_text


def plot_2Dboxes(ST,alfa=0.05,ifig=1, nsize=50):
  """
    Plot Statistics in the boxes
    ST - p-values from normality tests,
        p-val<alfa - rejet the H0 (data Normally distributed)
    using WOA topo interpolated from ETOPO
  """
  from copy import copy
  import matplotlib.colors as colors
  from matplotlib.patches import Rectangle

  pthtopo = '/data/ncei2/work/dmitry.dukhovskoy/topo/woa/'
  pthinp  = '/data/home004/hernan.garcia/work/FOR/WODQC/'
  ftopo = pthtopo+'etopo2woa025_bilinear.dat'

#  fid = open(ftopo,'rb')
#  HH = np.fromfile(fid, dtype='float', count=idmW*jdmW)
#  fid.close()
  dtype = np.dtype('float32')
  flon = pthtopo + 'lonWOA025.dat'
  with open(flon,'rb') as ff:
    LON = np.fromfile(ff,dtype)

  flat = pthtopo + 'latWOA025.dat'
  with open(flat,'rb') as ff:
    LAT = np.fromfile(ff,dtype)

  idm = LON.shape[0]
  jdm = LAT.shape[0]

  with open(ftopo,'rb') as ff:
    HH = np.fromfile(ff, dtype).reshape((jdm,idm))


#
# Read box coordinates:
  fbx = pthinp + 'outboxes.dat'
  print('Reading box coords')
  fid = open(fbx,'r')
  cc = -1
  XY = np.empty((0,4),float)
  for line in fid:
    aa = line.strip()
    cll = line.split()
    xb1 = float(cll[3])
    xb2 = float(cll[4])
    yb1 = float(cll[1])
    yb2 = float(cll[2])

    cc += 1
    XY = np.append(XY, np.array([[xb1,xb2,yb1,yb2]]),axis=0)

  fid.close()



  clrmp = copy(plt.cm.winter)
  clrmp.set_over(color=[0.5, 0.5, 0.5])

  print('Plotting Stat Map')
  plt.ion()
  fig1 = plt.figure(ifig, figsize=(8,8), constrained_layout=False)
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  cl1=-8000.
  cl2=0.

  im1 = plt.pcolormesh(LON,LAT,HH,shading='nearest',cmap=clrmp,\
															norm=colors.Normalize(vmin=cl1, vmax=cl2))
  im1.set_clim(cl1,cl2)
  ax1.axis('scaled')
  ax1.set(xlim=(-180,180),ylim=(-90,90))

  ax1.autoscale(enable=True, axis='both', tight=True)

#
# Plot boxes
  nbx = XY.shape[0]
  for ibx in range(nbx):
    x1 = XY[ibx,0]
    x2 = XY[ibx,1]
    y1 = XY[ibx,2]
    y2 = XY[ibx,3]
    dx = np.abs(x2-x1)
    dy = np.abs(y2-y1)
    pv = ST[ibx]

    if x1>180.:
      x1=x1-360.
    if x2>180.:
      x2=x2-360.

    if pv < alfa:
      nrm = 0
    else:
      nrm = 1
      clr = [0.8,0.6,0]

    xv=np.array([x1,x2,x2,x1,x1])
    yv=np.array([y1,y1,y2,y2,y1])
    lclr=[0.2,0.2,0.2]
    if nrm==0:
#      ax1.plot(xv,yv,'-',color=lclr)
      bb=0.
#      ax1.plot(x1,y1,'k.')
    else:
      pp1 = plt.Rectangle((x1,y1),dx,dy,color=clr,alpha=0.8)
      ax1.add_patch(pp1)

  ctl = 'Normally Distr. Data, Shap-Wilk, alpha={0:4.3f}, M-CarloSize={1}'.\
               format(alfa,nsize)
  ax1.set_title(ctl)
  fig1.colorbar(im1,ax=ax1,orientation='horizontal')

  plt.show()

  return

def plot_hist(Din,nbins=0, ifig=1, ctl='hist', btx='plot_fields.py'):
  if nbins==0:
    hst, bin_edges = np.histogram(Din,density=True)
  else:
    hst, bin_edges = np.histogram(Din,bins=nbins,density=True)

#  hst = hst/np.sum(hst)
  Nb = bin_edges.shape[0]-1
  xbin = np.zeros((Nb))
  db = bin_edges[1]-bin_edges[0]
  for ib in range(Nb):
    xbin[ib] = bin_edges[ib]+db/2.

  fig1 = plt.figure(ifig, figsize=(9,7), constrained_layout=False)
  plt.clf()
  ax1 = plt.axes([0.1, 0.2, 0.8, 0.6])
  ax1.bar(xbin,hst,width=db,color=(0.,0.5,0.9),edgecolor=(0.3,0.3,0.3)) 

  ax1.set_title(ctl)
  bottom_text(btx)  

  plt.show()

  return    



