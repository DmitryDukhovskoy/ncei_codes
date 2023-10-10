import numpy as np
import importlib

def inpolygon(xq, yq, xv, yv):
	""" 
	Function similar to matlab inpolygon
	based on interent stackoverflow
	returns in indicating if the query points specified by xq and yq 
	are inside or on the edge of the polygon area defined by xv and yv.
	"""
	from matplotlib import path
	shape = xq.shape
	xq = xq.reshape(-1)
	yq = yq.reshape(-1)
	xv = xv.reshape(-1)
	yv = yv.reshape(-1)
	q = [(xq[i], yq[i]) for i in range(xq.shape[0])]
	p = path.Path([(xv[i], yv[i]) for i in range(xv.shape[0])])

# Alternative:
#	shape = xq.shape
#	q = np.column_stack((xq.flatten(), yq.flatten()))
#	p = path.Path(np.column_stack((xv.flatten(), yv.flatten())))

	return p.contains_points(q).reshape(shape)

def rotate_vector(uin,vin,thtd):
	"""
	Rotate vector U(uin,vin)
	by angle thtd - in degrees
	"""
	tht = thtd*np.pi/180.
	R = [np.cos(tht), -np.sin(tht), np.sin(tht), np.cos(tht)] # rotates c/clkws if tht<0
	R = np.array(R).reshape(2,2)
	UV = np.array([uin,vin]).reshape(2,1)
	UVr = R.dot(UV)

	ur = UVr[0].item()
	vr = UVr[1].item()

	"""
	nf = 3
	arrowprops = dict(color='darkorange', linewidth=2)
	ax = compass(uin, vin, arrowprops, nf)

	nf = 4
	arrowprops = dict(color='blue', linewidth=2)
	ax = compass(ur, vr, arrowprops, nf)
	"""

	return ur, vr

def dist_sphcrd(xla1,xlo1,xla2,xlo2, Req=6371.0e3, Rpl=6357.e3):
	"""
#	this procedure calculates the great-circle distance between two
#	geographical locations on an ellipse using Lambert formula
#
#	lat-lon coordinates with its appropiate trigonometric
#	signs. 
#	INPUT: xla1, xlo1 - first point coordinates (latitude, longitude)
#				 xla2, xlo2 - second point
# all input coordinates are in DEGREES: latitude from 90 (N) to -90,
# longitudes: from -180 to 180 or 0 to 360,
# LAT2, LON2 can be either coordinates of 1 point or N points (array)
# in the latter case, distances from Pnt 1 (LAT1,LON1) to all pnts (LAT2,LON2)
# are calculated 
#	OUTPUT - distance (in m)
#	R of the earth is taken 6371.0 km
#
	"""
#	print("xla1=",xla1)
# breakpoint()
	xla1 = np.float64(xla1)
	xlo1 = np.float64(xlo1)
	xla2 = np.float64(xla2)
	xlo2 = np.float64(xlo2)

	if np.absolute(xla1).max() > 90.0:
		print("ERR: dist_sphcrd Lat1 > 90")
		dist = float("nan")
		return dist
	if np.absolute(xla2).max() > 90.0:
		print("ERR: dist_sphcrd Lat2 > 90")
		dist = float("nan")
		return dist

	cf = np.pi/180.
	phi1 = xla1*cf
	phi2 = xla2*cf
	lmb1 = xlo1*cf
	lmb2 = xlo2*cf
	dphi = abs(phi2-phi1)
	dlmb = abs(lmb2-lmb1)
	
	I0s = None
	if type(dphi) in (int,float):   # scalar input
		if dphi == 0.0 and dlmb == 0.0:
			dist_lmbrt = 0.0
			return dist_lmbrt
	else:											# N dim array
		if np.min(dphi) == 0.0 and np.min(dlmb[dphi==0])==0:
			I0s = np.argwhere((dphi==0.0) & (dlmb==0.0))
#
	Eflat = (Req-Rpl)/Req  # flatenning of the Earth
# Haversine formula to calculate central angle:
	aa1 = (np.sin(dphi/2.))**2
	aa2 = np.cos(phi1)*np.cos(phi2)*(np.sin(dlmb/2.))**2
	dsgm_hv = 2.*np.arcsin(np.sqrt(aa1+aa2))  # haversine form for central angle
#
# Reduced latitides - due to flattening:
	beta1 = np.arctan((1.-Eflat)*np.tan(phi1))
	beta2 = np.arctan((1.-Eflat)*np.tan(phi2))
	PP = 0.5*(beta1+beta2)
	QQ = 0.5*(beta2-beta1)
	X = (dsgm_hv-np.sin(dsgm_hv))*( (np.sin(PP))**2*(np.cos(QQ))**2 )/( (np.cos(dsgm_hv/2.))**2 )
	Y = (dsgm_hv+np.sin(dsgm_hv))*( (np.cos(PP))**2*(np.sin(QQ))**2 )/( (np.sin(dsgm_hv/2.))**2 )
#	if np.sin(dsgm_hv/2.) == 0.0:
#		breakpoint()

	dist_lmbrt = Req*(dsgm_hv-Eflat/2.*(X+Y))

	if np.min(dist_lmbrt)<0.0:
		print('WARNING: spheric distance <0: ',np.min(dist_lmbrt))

#	breakpoint()
	
	return dist_lmbrt

