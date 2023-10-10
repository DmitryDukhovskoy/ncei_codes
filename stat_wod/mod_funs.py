"""
  Functions and subroutines for WOD stat analysis
"""
import numpy as np

def read_NAdata_ascii(fnm):
  print('Reading NAtlantic {0}'.format(fnm))
  fid = open(fnm,'r')
  cc = -1
  for line in fid:
    aa  = line.strip()
    cll = line.split()
    x0  = float(cll[0])
    y0  = float(cll[1])
    dmm = float(cll[2])
    zz  = float(cll[3])

    cc += 1
    if cc == 0:
      X = np.array([x0])
      Y = np.array([y0])
      D = np.array([dmm])
      Z = np.array([zz])
    else:
      X = np.append(X,[x0])
      Y = np.append(Y,[y0])
      D = np.append(D,[dmm])
      Z = np.append(Z,[zz])

  print('{0} lines read'.format(cc+1))
  fid.close()

  return X,Y,D,Z

def read_data_ascii(fnm):
  print('Reading {0}'.format(fnm))
  fid = open(fnm,'r')
  cc = -1
  for line in fid:
    aa  = line.strip()
    cll = line.split()
    x0  = float(cll[0])
    y0  = float(cll[1])
    dmm = float(cll[2])

    cc += 1
    if cc == 0:
      X = np.array([x0])
      Y = np.array([y0])
      D = np.array([dmm])
    else:
      X = np.append(X,[x0])
      Y = np.append(Y,[y0])
      D = np.append(D,[dmm])

  print('{0} lines read'.format(cc+1))
  fid.close()

  return X,Y,D,cc+1

def norm_KStest(D):
  """
    Kolmogorov-Smirnov Normality test
  """
  from scipy import stats

# Generate a normal standardized distribution
# Seed for reproducibility
  rng   = np.random.default_rng(1300)
  NN = 1000
  distZ = rng.standard_normal(NN)
  chckZ = rng.standard_normal(NN)
#  ST = stats.ks_2samp(distZ,chckZ)
# Normalize tested data
  sgm = stats.tstd(D)
  mn  = np.nanmean(D)
  Dz = (D-mn)/sgm
  ST = stats.ks_2samp(distZ,Dz)

  return ST

def norm_ShapWtest(D):
  """
    Shapiro-Wilk test for normality
    works better for small sample size (<50)
  """
  from scipy import stats

# Generate a normal standardized distribution
# Seed for reproducibility
  f_chck=0
  if f_chck>0:
    rng   = np.random.default_rng(1300)
    NN = 1000
    chckZ = rng.standard_normal(NN)
    ST = stats.shapiro(chckZ)
    ST.pvalue

# Normalize tested data
  sgm = stats.tstd(D)
  mn  = np.nanmean(D)
  Dz = (D-mn)/sgm
  ST = stats.shapiro(Dz)

  return ST

def MCarlo(Din,nsize):
  """
    A simple Monte-Carlo simulation to increase sample size 
    nsize = desired sample size 
  """
  Nin = Din.shape[0]
  hst, bin_edges = np.histogram(Din)
  dh=bin_edges[1]-bin_edges[0]
  hst = hst/(np.sum(hst)*dh)
  Dout = np.zeros((nsize))

  lh = hst.shape[0]
  bin_mid = np.zeros((lh))
  for kk in range(lh):
    bin_mid[kk] = bin_edges[kk]+0.5*dh

#
# Interpolate missing probabilities:
# assumed 1st and last bin values exist
  lDv = 1000
  ddh = (bin_edges[-1]-bin_edges[0])/float(lDv)
#  Dval = np.arange(bin_edges[0],bin_edges[-1],ddh)
  Dval = np.arange(0,1,1./float(lDv))*\
         (bin_edges[-1]-bin_edges[0])+bin_edges[0]

  Prb = np.zeros(lDv)
  In = np.where(hst>0)[0]
  for kk in range(lDv):
    d0 = Dval[kk]
    ii = np.max(np.where(bin_edges<=d0))
    if hst[ii] <= 1.e-20:
      imm=In-ii
      ix1=np.max(np.where(imm<0))
      i1=In[ix1]
      i2=In[ix1+1]
      h1=hst[i1]
      h2=hst[i2]
      L = h1*(ii-i2)/(i1-i2)+h2*(ii-i1)/(i2-i1)
      Prb[kk] = L
    else:
      Prb[kk] = hst[ii]

# Smooth probability density function:
  mnwnd = round(0.05*lDv)
  mnwnd = round(mnwnd/2)*2+1  # make it odd
  mnwnd = np.max([mnwnd,3])
  Prb = runmn(Prb,mnwnd=mnwnd)
 # Prb = runmn(Prb,mnwnd=mnwnd)
  Prb = Prb/(np.sum(Prb)*ddh)  # for plotting
  Prb0 = Prb/np.sum(Prb)

  Dout = np.random.choice(Dval, size=nsize, p=Prb0)     

  f_chck=0
  if f_chck:
    plt.figure(1)
    plt.plot(bin_mid,hst)
    plt.plot(Dval,Prb) 
 
  return Dout 


def runmn(B0, ZZ=None, mnwnd=5):
  """
    Running mean, mnwnd - averaging window should be odd
    uneven intervals are allowed
    ZZ - intervals, values B0 are in the mid-points
         i.e. ZZ length = B0 length+1
    ZZ=[] if intervals are equi-distant
  """
  import math
  if math.floor((mnwnd-1)/2)*2 != mnwnd-1:
    print('Averaging window {0} should be odd, adding 1==> {1}'.
           format(mnwnd,mnwnd+1))
    mnwnd = mnwnd+1

  dk = math.floor((mnwnd-1)/2)
  nlev = B0.shape[0]
  BF = np.zeros(nlev)
  for kk in range(nlev):
    k1 = kk-dk
    k2 = kk+dk
    k1 = max(0,k1)
    k2 = min(nlev-1,k2)

    dz_sum = 0.
    b_sum = 0.
    for ksb in range(k1,k2+1):
#      print('Averaging ksb = {0}'.format(ksb))
      if ZZ is not None: 
        dz = np.abs(ZZ[ksb]-ZZ[ksb+1])
      else:
        dz=1.
      dz_sum = dz_sum+dz
      b_sum = b_sum+B0[ksb]*dz
#    print('  ')
    BF[kk] = b_sum/dz_sum

  return BF




