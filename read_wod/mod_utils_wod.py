""" 
  Subroutines for reading WOD headers and profiles

  Dmitry Dukhovskoy, NOAA NESDIS NCEI
  June 2022
"""

import importlib
import timeit
import binascii
import struct
import sys
import numpy as np

#
def get_obsparam():
# Get list of obs parameters and their numeric codes
  PARAM={
    "Temperature"  : 1,
    "Salinity"     : 2,
    "Oxygen"       : 3,
    "Phosphate"    : 4,
    "TotalPhos"    : 5,
    "Silicate"     : 6,
    "Nitrite"      : 7,
    "Nitrate"      : 8,
    "pH"           : 9,
    "Ammonia"      : 10,
    "Chlorophyll"  : 11,
    "Phaeophytin"  : 12,
    "PrimaryProd"  : 13,
    "Biochem"      : 14,
    "LightC14"     : 15,
    "DarkC14"      : 16,
    "Alkalinity"   : 17,
    "POC"          : 18,
    "DOC"          : 19,
    "pCO2"         : 20,
    "tCO2"         : 21,
    "XCO2sea"      : 22,
    "NO2NO3"       : 23,
    "Transmissiv"  : 24,
    "Pressure"     : 25,
    "Conductivit"  : 26,
    "CO2warm"      : 27,
    "xCO2atm"      : 28,
    "AirPress"     : 29,
    "Latitude"     : 30,
    "Longitude"    : 31,
    "JulianDay"    : 32,
    "Tritium"      : 33,
    "Helium"       : 34,
    "DeltaHe3"     : 35,
    "DeltaC14"     : 36,
    "DeltaC13"     : 37,
    "Argon"        : 38,
    "Neon"         : 39,
    "CFC11"        : 40,
    "CFC12"        : 41,
    "CFC113"       : 42,
    "Oxy18"        : 43
    }
  return PARAM 

def get_obsparam_code(varname):
# Find Code # given parameter name
  PARAM = get_obsparam()
  if varname not in PARAM:
    raise NameError(varname+' is not a valid name for observed parameters')

  list_codes = list(PARAM.values())
  list_names = list(PARAM.keys())
  ii = list_names.index(varname)
  vcode = list_codes[ii]

  return vcode  

def get_obsparam_name(vcode):
# Find obs parameter name given its code
  PARAM = get_obsparam()
  list_codes = list(PARAM.values())
  list_names = list(PARAM.keys())
  ii = list_codes.index(vcode)
  vname = list_names[ii]

  return vname

def probe_obs(prbname):
# Library of all probe types with obs. parameter codes
# Not finished - need to add all missing: DRB, MRB, Argo, etc.
  xbt = {
    "Vcode" : [1]
    }

  ctd = {
    "Vcode" : [1,2,3,8,9,11,20,24,25,26]
    }

  dmm1 = list(range(1,27))
  dmm2 = list(range(33,44))
  osd = {
    "Vcode" : dmm1+dmm2
    }

  dmm1 = list(range(1,33))
  sur = {
    "Vcode" : dmm1
    }

  mrb = {
    "Vcode" : [1,2,3,8,9,11,20,24,25,26]
    }

  pfl = {
    "Vcode" : [1,2,3,8,9,11,20,24,25,26]
    }

  uor = {
    "Vcode" : [1,2,3,11,24,25,26,30,31,32]
    }

  drb = {
    "Vcode" : [1,2,3,11,24,25,26]
    }

  gld = {
    "Vcode" : [1,2,3,11,24,25,26,30,31,32]
    }

  PRBC={
    "XBT" : xbt,
    "CTD" : ctd,
    "OSD" : osd,
    "SUR" : sur,
    "MRB" : mrb,
    "PFL" : pfl,
    "UOR" : uor,
    "DRB" : drb,
    "GLD" : gld
   }

  list_names = list(PRBC.keys())
  LCd_obs = list(PRBC.get(prbname).get("Vcode"))
  LNm_obs = []

  nvar = len(LCd_obs)
  for jj in range(nvar):
    vcode = LCd_obs[jj]
    vname = get_obsparam_name(vcode)
    LNm_obs.insert(jj,vname)
  
  return LCd_obs, LNm_obs
#
# directories, data format, etc for specified probe type
def get_probe_info(prb):
# For 1-variable case, e.g. XBT
  pthdata = '/data/ncei2/OCL/data/'+prb+'/'
  fhdr = pthdata+'headers'
#
# After Python is updated to 3.10 or later rewrite 
# using match / case algorithm
# for now, if/else/ ...
  if prb == "XBT":
    fmt = '2sffiiifiiiiii'  # format of a record in header file
    rec_size = struct.calcsize(fmt)

    PRBINFO={
      "probe_nm" : prb,
      "pth"      : pthdata,
      "fl_head"  : fhdr,
      "nobtypes" : 1,
      "obs1"     : "temp",
      "fmt_head" : fmt,
      "rec_size" : rec_size
    }

  else:
    raise NameError(prb+' is not supported in get_probe_info')

  return PRBINFO 
    
def read_header_1var(PRBINFO,krcrd):
# This function is set for 1 variable probe
# one variable probe pointers:
#1 pointer: depth/variable_one_observed
#1 pointer: variable_one_standard
#1 pointer: second_header
  fhdr = PRBINFO.get("fl_head")
  fmt  = PRBINFO.get("fmt_head")
  rec_size = PRBINFO.get("rec_size")
  nskip = rec_size*(krcrd-1)

  print("Reading "+fhdr)
  fid = open(fhdr,'rb')
  try:
    fid.seek(nskip)
    bd = fid.read(rec_size)
    tmp = struct.unpack(fmt,bd)

    HEAD={
      "CC"  : tmp[0].decode(),  # convert binary --> string
      "lat" : tmp[1],
      "lon" : tmp[2],
      "yr"  : tmp[3],
      "mo"  : tmp[4],
      "md"  : tmp[5],
      "gtm" : tmp[6],
      "cid" : tmp[7],
      "nobs": tmp[8],
      "nstd": tmp[9],
      "p1"  : tmp[10],   # pointer - depth and pointer of T obs.levels
      "p2"  : tmp[11],   # pointer - variable at standard levels
      "p3"  : tmp[12]    # second header pointer
    }

  finally:
    fid.close()

  return HEAD

def read_header(prb,krcrd):
  """
# Read header for specified multi-variable probe
    Multi-variable probe:
    1 pointer: depth
    N pointers: variables_[1 to N]_observed
    N pointers: variabels_{1 to N]_standard
    1 pointer: second_header

    One variable probes (XBT) have different pointers' order:
    The difference between one variable probes and multivariable probes 
    is the one variable probe use the same pointer for depth and variable.
    So one variable probe as follows:
      1 pointer: depth/variable_one_observed
      1 pointer: variable_one_standard
      1 pointer: second_header
  """
  pthdata = '/data/ncei2/OCL/data/'+prb+'/'
  fhdr = pthdata+'headers'

  LCd_obs, LNm_obs = probe_obs(prb) # list of probe obs codes and names
  nvar = len(LCd_obs)   # # of observed variables

  class ProbePntr:
    def __init__(self,tmp,nvar,LNm_obs,uniqnmb=0):
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
      self.uniqnmb= uniqnmb   # unique cast number
      self.VarNm  = LNm_obs
      if nvar == 1:          # special case for XBT and 1-var probes
        p1        = 10
        p2        = 11 
      else:
        p1        = 11
        p2        = p1+nvar

      self.pobs   = tmp[p1:p2] # pointer to var=1,...,nvar obs depths
      p1          = p2
      p2          = p1+nvar
      self.pstd   = tmp[p1:p2] # pointer to var=1,...,nvar stand depths
      self.shead  = tmp[p2]    # pointer to second header

# Format for main header information:
# CC code - 2 char
# lat, lon - 32-bit float little endian
# year, mo, day - 32-bit int
# GMT time - 32-bit float
# Cruise ID - int
# obs - int, # of observed levels in profile
# nstd -  # of st levels in profile - int
# pointers = depth, var1 obs levels, var1 st level, ... 
  if prb == "XBT":
    fmt = '2sffiiifiiiiii'  # format of a record in header file
  else:
    fmt = '2sffiiifiiii'
    for ii in range(nvar):
      fmt = fmt+'i'         # pointers vars obs levels
    for ii in range(nvar):
      fmt = fmt+'i'         # pointers vars std levels
    fmt = fmt +'i'     # second header pointer

  rec_size = struct.calcsize(fmt)
  nskip = rec_size*(krcrd-1)

  fid = open(fhdr,'rb')
#
# Determine the size of file
  fid.seek(0,2)
  fend = fid.tell()
  nrec_tot = fend/rec_size
#  print('Total records in the file {0:.0f}'.format(nrec_tot))
  if krcrd > nrec_tot:
    print(' Recrod to read is at the EOF - no data returned')
    print(' Recrod to read={0}, EOF record ={1}'.\
          format(krcrd,nrec_tot))
    PRBINFO = []
    PRB = []
    
    return PRBINFO, PRB
    
  try:
    fid.seek(nskip)
    bd = fid.read(rec_size)
    tmp = struct.unpack(fmt,bd)

    uniqnmb = get_uniqnmb_cast(prb,krcrd) # unique cast nmb

    PRB = ProbePntr(tmp,nvar,LNm_obs,uniqnmb=uniqnmb)
  except:
    print('ERROR: binary string ={0}'.format(bd))
    print('ERROR: converted = {0}'.format(tmp))
    print('ERROR: Cannot read {0}'.format(fhdr))
     
  finally:
    fid.close()

  nvar = len(LNm_obs)
  PRBINFO={
    "probe_nm" : prb,
    "pth"      : pthdata,
    "fl_head"  : fhdr,
    "nobtypes" : nvar,
    "obs"      : LNm_obs,
    "fmt_head" : fmt,
    "rec_size" : rec_size
  }

  return PRBINFO, PRB


def read_obsdepths(PRBINFO,PRB,fecho=1):
  pth   = PRBINFO.get("pth")
  fdpth = pth+"depths"
  nzz   = PRB.nobs
  pntrz = PRB.pdepth  # pointer for depth record
  nskip = pntrz*struct.calcsize('f')
  rec_size = nzz*struct.calcsize('f')

  fmtz = ''
  for iz in range(nzz):
    fmtz = fmtz+'f'

  zz = np.empty
  if fecho == 1:
    print("Reading "+fdpth)

  fid = open(fdpth,'rb')
  try:
    fid.seek(nskip)
    bd = fid.read(rec_size)
    zz = struct.unpack(fmtz,bd)

  finally:
    fid.close()
   
  zz =  np.asarray(zz)
  return zz


def read_obsdepths_1var(PRBINFO,HEAD):
  pth   = PRBINFO.get("pth")
  fdpth = pth+"depths"
  nzz   = HEAD.get("nobs")
  pntrz = HEAD.get("p1")  # pointer for depth record
  nskip = pntrz*struct.calcsize('f')
  rec_size = nzz*struct.calcsize('f')

  fmtz = ''
  for iz in range(nzz):
    fmtz = fmtz+'f'

  zz = np.empty
  print("Reading "+fdpth)
  fid = open(fdpth,'rb')
  try:
    fid.seek(nskip)
    bd = fid.read(rec_size)
    zz = struct.unpack(fmtz,bd)

  finally:
    fid.close()
   
  zz =  np.asarray(zz)
  return zz


def read_varobs_1var(PRBINFO,HEAD,varnm):
# one variable probe ! 
# Observations are written in this format:
# # of sig figures, in int
# observation
  pth   = PRBINFO.get("pth")+"observed/"
  nobs  = PRBINFO.get("nobtypes")
  fvar  = pth+varnm
  nzz   = HEAD.get("nobs")
#  pntrv = HEAD.get("p2")  # pointer for obs var. record
  pntrv = HEAD.get("p1")  # pointer for obs var. record = depth pointer
  nskip = pntrv*struct.calcsize('if')
  
  fmtz = ''
  for iz in range(nzz):
    fmtz = fmtz+'i'
    fmtz = fmtz+'f'

  rec_size = struct.calcsize(fmtz)

  var = np.empty
  print("Reading "+fvar)
  fid = open(fvar,'rb')
  try:
    fid.seek(nskip)
    bd  = fid.read(rec_size)
    var = struct.unpack(fmtz,bd)

  finally:
    fid.close()

  dmm = np.asarray(var)
  dmm = dmm.reshape(nzz,2)
  var = dmm[:,1]

  return var

def read_varstd_1var(PRBINFO,HEAD,varnm):
# one-variable probe
# Read variable varnm at standard depths
# Observations are written in this format:
# # of sig figures, in int
# observation
  pth   = PRBINFO.get("pth")+"standard/"
  nobs  = PRBINFO.get("nobtypes")
  fvar  = pth+varnm
  nzz   = HEAD.get("nstd")
#  pntrv = HEAD.get("p3")  # pointer for obs var. record
  pntrv = HEAD.get("p2")  # for 1var Probe
  nskip = pntrv*struct.calcsize('if')
  
  fmtz = ''
  for iz in range(nzz):
    fmtz = fmtz+'i'
    fmtz = fmtz+'f'

  rec_size = struct.calcsize(fmtz)

  var = np.empty
  print("Reading "+fvar)
  fid = open(fvar,'rb')
  try:
    fid.seek(nskip)
    bd  = fid.read(rec_size)
    var = struct.unpack(fmtz,bd)

  finally:
    fid.close()

  dmm = np.asarray(var)
  dmm = dmm.reshape(nzz,2)
  var = dmm[:,1]

  return var

def standard_depths(nzlev=1000):
# Original standard depths - main database
# From /data/ncei2/OCL/sys.inf/standard_depths_orig.dat
  ar1 = np.array([0.,10.,20.,30.])
  ar2 = np.arange(50,150,25)
  ar3 = np.arange(150,300,50)
  ar4 = np.arange(300,1500,100)
  ar5 = np.arange(1500,2000,250)
  ar6 = np.arange(2000,9001,500)

  Zstd = np.concatenate((ar1,ar2,ar3,ar4,ar5,ar6))

  nz = Zstd.size
  if nzlev < nz:
    Zstd = Zstd[0:nzlev]

  return Zstd

def find_var_indx(PRB,varnm):
# Variable indx in the list of types of variables for this probe
# vindx = [0,...] - array index
  VarNm = PRB.VarNm

  try:
    vindx = VarNm.index(varnm)
  except:
    raise NameError(varname+' is not found in this probe')
  return vindx

def read_varobs(PRBINFO,PRB,varnm,fecho=1):
# Read variable varnm at observed depths
  vindx = find_var_indx(PRB,varnm)
  pth   = PRBINFO.get("pth")+"observed/"
  nobs  = PRBINFO.get("nobtypes")
  fvar  = pth+varnm
  nzz   = PRB.nobs
  pntrv = PRB.pobs[vindx]

  if pntrv<0:
#    raise NameError(varnm+' pointer <0: no observation ')
    print(varnm+' pointer <0: no observation ')
    var = np.empty(0)
    return var

  nskip = pntrv*struct.calcsize('if')

  fmtz = ''
  for iz in range(nzz):
    fmtz = fmtz+'i'
    fmtz = fmtz+'f'

  rec_size = struct.calcsize(fmtz)
  
  var = np.empty(0)
  if fecho == 1:
    print("Reading "+fvar)

  fid = open(fvar,'rb')
  try:
    fid.seek(nskip)
    bd  = fid.read(rec_size)
    var = struct.unpack(fmtz,bd)

  finally:
    fid.close()

  dmm = np.asarray(var)
  dmm = dmm.reshape(nzz,2)
  var = dmm[:,1]

  return var

def read_varstd(PRBINFO,PRB,varnm):
# Read variable varnm at standard depths
  vindx = find_var_indx(PRB,varnm)
  pth   = PRBINFO.get("pth")+"standard/"
  nobs  = PRBINFO.get("nobtypes")
  fvar  = pth+varnm
  nzz   = PRB.nstd
  pntrv = PRB.pstd[vindx]

  if pntrv<0:
    raise NameError(varnm+' pointer <0: no observation ')

  nskip = pntrv*struct.calcsize('if')
  
  fmtz = ''
  for iz in range(nzz):
    fmtz = fmtz+'i'
    fmtz = fmtz+'f'

  rec_size = struct.calcsize(fmtz)

  var = np.empty
  vindx = find_var_indx(PRB,varnm)
  print("Reading "+fvar)
  fid = open(fvar,'rb')
  try:
    fid.seek(nskip)
    bd  = fid.read(rec_size)
    var = struct.unpack(fmtz,bd)

  finally:
    fid.close()

  dmm = np.asarray(var)
  dmm = dmm.reshape(nzz,2)
  var = dmm[:,1]

  return var

#
def get_uniqnmb_cast(prbnm,krcrd):
# Get unique cast number for given probe
# and sequential record in the probe
  pthdata = '/data/ncei2/OCL/data/'+prbnm+'/masks/other/'
  funiq = pthdata+'uniq.s'
  fmtz = 'i'
  rec_size = struct.calcsize(fmtz) 
  nskip = rec_size*(krcrd-1)

  fid = open(funiq,'rb')
  try:
    fid.seek(nskip)
    bd  = fid.read(rec_size)
    uniqnmb = struct.unpack(fmtz,bd)[0]

  finally:
    fid.close()
  
  return uniqnmb 


 
