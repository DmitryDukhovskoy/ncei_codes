"""
	Read HYCOM *.[ab] files
	grid, topo, 
	
	Dmitry Dukhovskoy, FSU COAPS
"""
import os
from netCDF4 import Dataset as ncFile
import numpy as np
import sys

def read_grid_topo(pthtopo,ftopo,fgrid,f_lon=180,pthgrid=None, dpth_neg=True):
# Get topo and lon/lat (plon/plat) for plotting 
# To get other lon/lat (ulon, ulat, vlon, vlat) - need to 
# read full grid file
#
# f_lon = 180 - convert lon to [-180, 180]
#       = 360 - convert lon to [0, 360]
#       = else - leave as in the grid file
#
# read HYCOM grid and topo files *.[ab]
	if pthgrid is None: 
	  pthgrid = pthtopo
	fltopoa = pthtopo+ftopo+'.a'
	fltopob = pthtopo+ftopo+'.b'
	flgrida = pthgrid+fgrid+'.a'
	flgridb = pthgrid+fgrid+'.b'

# Read dim, lon, lat
	fgb = open(flgridb,'r')
	fgb.seek(0)
	data = fgb.readline().split()
	IDM = int(data[0])
	data = fgb.readline().split()
	JDM = int(data[0])
	IJDM = IDM*JDM
	fgb.close()

	npad =4096-IJDM%4096

	print('Reading HYCOM grid/topo:{0} {1} '.format(fgrid,ftopo))
	print('Grid: IDM={0}, JDM={1}'.format(IDM,JDM))

# Read direct access binary file
	fga = open(flgrida,'rb')
	fga.seek(0)
	plon = np.fromfile(fga,dtype='>f',count=IJDM)

	if f_lon == 180:
	  while any(plon > 180.):
	    plon = np.where(plon<=180.,plon,plon-360.)

	  while any(plon < -180.):
	    plon = np.where(plon >= -180., plon, plon+360.)

	if f_lon == 360:
	  while any(plon < 0.):
	    plon = np.where(plon >= 0., plon, plon+360.)
	  while any(plon > 360.):
	    plon = np.where(plon <= 360., plon, plon-360.)

	plon = plon.reshape((JDM,IDM))


	fga.seek(4*(npad+IJDM),0)
	plat = np.fromfile(fga, dtype='>f',count=IJDM)
	plat = plat.reshape((JDM,IDM))

	fga.close()

	print('Min/max lon = {0}, {1}, Min/max lat = {2}, {3}'.format(np.min(plon),\
	      np.max(plon),np.min(plat), np.max(plat)))

# Read bottom topography:
# Big endian float 32
	fbt = open(fltopoa,'rb')
	fbt.seek(0)
	HH = np.fromfile(fbt, dtype='>f', count=IJDM)
	HH = HH.reshape((JDM,IDM))
	fbt.close()

	if dpth_neg:
		HH[HH<1.e10] = -1.*HH[HH<1.e10]
		HH[HH>1.e10] = 100.

	print('Min/max Depth = {0}, {1}'.format(np.min(HH),np.max(HH)))

	return plon,plat,HH


def read_topo(pthtopo,ftopo,IDM,JDM,dpth_neg=True):
# Get topo only need to know I and J dimensions
#
# f_lon = 180 - convert lon to [-180, 180]
#       = 360 - convert lon to [0, 360]
#       = else - leave as in the grid file
#
# read HYCOM grid and topo files *.[ab]
	fltopoa = pthtopo+ftopo+'.a'
	fltopob = pthtopo+ftopo+'.b'

	IJDM = IDM*JDM
	npad =4096-IJDM%4096

	print('Reading HYCOM topo {0} '.format(ftopo))
	print('Grid: IDM={0}, JDM={1}'.format(IDM,JDM))

# Read bottom topography:
# Big endian float 32
	fbt = open(fltopoa,'rb')
	fbt.seek(0)
	HH = np.fromfile(fbt, dtype='>f', count=IJDM)
	HH = HH.reshape((JDM,IDM))
	fbt.close()

	if dpth_neg:
		HH[HH<1.e10] = -1.*HH[HH<1.e10]
		HH[HH>1.e10] = 100.

	print('Min/max Depth = {0}, {1}'.format(np.min(HH),np.max(HH)))

	return HH


def readnc_grid_topo(pthgrid,flgrid,f_lon=180,dpth_neg=True):
	"""
	Read grid topo from NetCDF 
	If dpth_neg: make depths <0 and land 100 (default)
	If not dpth_neg: depths>0, land - huge 
	"""
	from netCDF4 import Dataset as NetCDFFile

	hg = 1.e30
	nc = NetCDFFile(pthgrid+flgrid)
	lonc = nc.variables['Longitude'][:]
	latc = nc.variables['Latitude'][:]
	HH   = nc.variables['Bathymetry'][:]
	print('Reading ',pthgrid+flgrid)

	if f_lon == 180:
	  while any(lonc.flatten() > 180.):
	    lonc = np.where(lonc<=180.,lonc,lonc-360.)

	  while any(lonc.flatten() < -180.):
	    lonc = np.where(lonc >= -180., lonc, lonc+360.)

	if f_lon == 360:
	  while any(lonc.flatten() < 0.):
	    lonc = np.where(lonc >= 0., lonc, lonc+360.)
	  while any(lonc.flatten() > 360.):
	    lonc = np.where(lonc <= 360., lonc, lonc-360.)

	if dpth_neg:
		if np.min(HH) > 0.:
			HH[HH<1.e10] = -1.*HH[HH<1.e10]
			HH[HH>1.e10] = 100.
	else:
		if np.min(HH) < 0.:
			HH[HH>=0] = hg
			HH[HH<0]  = -1.*HH[HH<1.e10]

	print('Min/max lon = {0}, {1}, Min/max lat = {2}, {3}'.format(np.min(lonc),\
				np.max(lonc),np.min(latc), np.max(latc)))
	print('Min/max Depth = {0}, {1}'.format(np.min(HH),np.max(HH)))

	return lonc,latc,HH

####
def land3d_pmask(fina,finb, fland=0, fsea=1, fld='temp'):
	"""
	Create 3D land mask for all model layers 
	for p-point variables
	default: land = 0, sea = 1
	"""
 	
	print('land3d_pmaks:  Creating 3D Land mask p-points')
	F, idm, jdm, kdm = read_hycom(fina,finb,'temp')
	F[F>1.e20] = np.nan
	
	Lmsk = np.zeros((jdm,idm,kdm))
	for kk in range(kdm):
		aa = np.squeeze(F[:,:,kk])
		aa = np.where(np.isnan(aa),fland,fsea)
		Lmsk[:,:,kk] = aa

	return Lmsk


######
def read_hycom(fina,finb,fld,Rtrc=None,rLayer=None):
	"""
	reads hycom binary archive files (model output), 
	returns specified field 'fld'
	and dimensions of the grid
	Rtrc - passive tracer # to read, if there are tracers 
	rLayer - layer number to read, otherwise all layers are read - more time/memory
					numbering is lr=1,...,nlayers in the model


	Based on matlab function
%%  Dmitry Dukhovskoy, FSU, 2010
%  2017: added options for nlayer, n tracers
% if several tracers, options are:
% read tracer N: 'r_tracer',1
% read all tracers by default
% any variable can be read in 1 layer Nl:
%       'r_layer',1
%
% If all tracers are read, recommended to specify 
% only 1 layer to read 
%
	"""
	fgb = open(finb,'r')
	fgb.seek(0)
	nl0 = 0
	while nl0 < 100:
		nl0 += 1
		data = fgb.readline().split()
		adim = data[1]
		ii = adim.find('idm')
		if ii>0: 
			break

	if ii<0:
		fgb.close()
		sys.exit('No idm found: Reading ',finb)	

	IDM = int(data[0])
	data = fgb.readline().split()
	JDM = int(data[0])
	IJDM = IDM*JDM
#  fgb.seek(0)

	npad =4096-IJDM%4096

	print('Reading HYCOM :{0} '.format(finb))
	print('Grid: IDM={0}, JDM={1}'.format(IDM,JDM))

	aa = fgb.readline().split()

# Find location of the requested field:
	cntr= 0
	nf = len(fld)
	FLOC=[]
	while cntr<1e6:
		aa = fgb.readline().split()
		if len(aa) == 0:  # end of file
			break
		cntr += 1
		aname = aa[0]
		ii = aname.find(fld)
		if ii >= 0:
			FLOC.append(cntr)

	fgb.close()

	nrec = len(FLOC)
	if nrec == 0:
		sys.exit('read_hycom: Field {0} not found in {1}'.format(fld,finb))

# N. of v. layers
	"""
 If fld = tracer and # of tracers >1
 need to distinguish # of layers 
 vs # of tracers*layers
 if strmatch(fld,'tracer')
	"""
	ll = len(FLOC)
	FLOC = np.array(FLOC)
	if ll == 1:
		nVlev = 1
		nTR = 1
	else:	
		dI = np.diff(FLOC)
		dmm = np.where(dI>1)
		dindx = dmm[0]
		nTR = dindx[0]+1         # N. of tracers in 1 layer
		nVlev = ll/nTR				# # of v. layers

#	breakpoint()
	if nTR != 1:
		print('read_hycom: Found {0} variables {1}  per layer'.format(nTR,fld))

# Find layers to read, if specified
# and tracers, if specified
	lr1=-1
	lr2=-1
	if Rtrc is not None:
		if nTR < Rtrc:
			sys.exit('Number of saved tracers {0} < requested {1}'.format(nTR,Rtrc))
		dmm = np.copy(FLOC)
		FLOC = dmm[Rtrc-1::nTR]

		if lr1 < 0 or lr2 < 0 :
			lr2 = FLOC.shape[0]
			ll = lr2

	if rLayer is not None:
		lr1 = rLayer
		lr2 = lr1

# If a layer Number not specified - read all
	if lr1 < 0 or lr2 < 0:
		lr1 = 1
		lr2 = ll

	print('Reading {0}, Layers: {1}-{2}'.format(fld,lr1,lr2))

	fga = open(fina,'rb')
	F = []
	ccL = -1
	for ii in range(lr1,lr2+1):
		fga.seek(0)
		k0 = FLOC[ii-1]-1
		fga.seek(k0*(npad+IJDM)*4,0)
		dmm = np.fromfile(fga, dtype='>f',count=IJDM) # read 1 layer
		dmm = dmm.reshape((JDM,IDM))
		ccL += 1
#		print('lr={0}, ccL={1}'.format(ii,ccL))
#		breakpoint()
		if ccL == 0:
			F = np.copy(dmm)
		else:
			F = np.dstack((F,dmm))

	if ll == 0:
		print('!!! read_hycom: {0} not found in {1} ERR'.format(fld,fina))
		print('!!! read hycom: check fields in ',finb)
		
	fga.close()

	return F, IDM, JDM, ll

		
