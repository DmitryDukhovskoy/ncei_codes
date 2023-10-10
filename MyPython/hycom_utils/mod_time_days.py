import datetime
from datetime import date
from datetime import datetime as dtime
from datetime import timedelta
import calendar
import numpy as np

# Reference date for calculating day numbers
# Same as in Matlab shifted 1 year because python dtime 
# is not defined for year 0
Mday = dtime(1,1,1)

def get_dtime(yr,mo=0,mday=0,jday=None):
	"""
	Find hycom date
	specify either month/month day or Julian day (year day)
	"""
	dH1 = dtime(1900,12,31)

# Define requested date
	if jday is not None:
		dj1 = dtime(yr,1,1)
		dnmb = dj1+timedelta(days=jday-1)
		yr = int(dnmb.strftime('%G'))
		mo = int(dnmb.strftime('%m'))
		mday = int(dnmb.strftime('%d'))
	else:
		dnmb = dtime(yr,mo,mday)

	dlt = dnmb-dH1
	day0 = dlt.days

# Start end dates of the current month:
	mday_last = calendar.monthrange(yr,mo)[1]
	dnmbS = dtime(yr,mo,1)
	dnmbE = dtime(yr,mo,mday_last)

	dayS = (dnmbS-dH1).days
	dayE = (dnmbE-dH1).days

	return dayS, dayE, day0


def datenum(yr,mo=0,mday=0,hr=0,mint=0, jday=None):
	"""
	Similar to matlab datnum
	converts calendar dates to sum number starting at day=t0
	"""

# Define requested date
	if jday is not None:
		dj1 = dtime(yr,1,1)
		dnmb = dj1+timedelta(days=jday-1)
		yr = int(dnmb.strftime('%G'))
		mo = int(dnmb.strftime('%m'))
		mday = int(dnmb.strftime('%d'))
		fr_day = jday-np.floor(jday)
		hr = int(np.floor(fr_day*24.))
		mint = int(np.round((fr_day*24.-hr)*60.))
	else:
		dnmb = dtime(yr,mo,mday,hr,mint)

	dlt = dnmb-Mday
	day0 = dlt.days+1
	scnd0 = dlt.seconds

	sec_day = 3600.*24.
	dnmb0 = day0+scnd0/sec_day

	return dnmb0

def datevec(dnmb):
	"""
	Similar to matlab datevec: given datenumber return vector of Yr, mo, day, ...
	"""
	yr0 = int(Mday.strftime('%G'))
	mo0 = int(Mday.strftime('%m'))
	mday0 = int(Mday.strftime('%d'))
	hr0   = 0 
	mint0 = 0

	yr = yr0
	dmm = 0
	while dmm < dnmb:
		dmm_old = dmm
		if calendar.isleap(yr):
			dmm += 366
		else:
			dmm += 365
		yr += 1

	yr = yr-1
	dmm = dmm_old

	imo = 1
	while dmm < dnmb:
		mday_last = calendar.monthrange(yr,imo)[1]
		dmm += mday_last
		imo += 1

	imo = imo-1
	if imo > 12:
		print('ERR: datevec month {0} > 12'.format(imo))
	dmm = dmm-mday_last

	mday = dnmb-dmm
	fr_day = mday-np.floor(mday)
	mday = int(np.floor(mday))
	hr = int(np.floor(fr_day*24.))
	mint = int(np.round((fr_day*24.-hr)*60.))
	
	dj1 = datenum(yr,1,1)
	jday = int(dnmb-dj1+1)

	DV = [yr,imo,mday,hr,mint,jday]

	return DV
	






