"""
  Functions for solving Sturm-Liouville Problem
  Dmitry Dukhovskoy, NOAA NESDIS NCEI
  July 2022
"""
import numpy as np

def form_mtrxA(Z_phi,kbtm,zbtm):
  """
    d2(Phi)/dz2 = AA*Phi(zz)
    AA is a tridiagonal matrix
    for approximating 2nd derivative - centered difference
    with uneven stepping

    AA[0,:] = subdiagonal elements
    AA[1,:] = diagonal elements
    AA[2,:] = superdiagonal elements

    Eigenfunctions are defined in the mid-points
    of the vertical layers - for centered difference

    eigenfunctions defined in the water column
    at the bottom - Phi=0

    Input:
      Z_phi - mid-point depths, np array kdm pnts
              depths are negative, indexing from surface downward
      kbtm - bottom index
  """
  
  kdm = Z_phi.shape[0]
  AA  = np.zeros((kbtm+1,3)) # index+1 because indexing starts from 0
  Dk  = np.abs(Z_phi[0])
  Dkp1= (Z_phi[0]-Z_phi[1])
  check_dz(Dk,Dkp1,0)
  AA[0,0] = -2./(Dk*Dkp1)
  AA[0,1] = 2./(Dkp1*(Dk+Dkp1))
  for kk in range (1,kbtm):
    Dk   = (Z_phi[kk-1]-Z_phi[kk])
    Dkp1 = (Z_phi[kk]-Z_phi[kk+1])
    check_dz(Dk,Dkp1,kk)

    AA[kk,0] = 2./(Dk*(Dkp1+Dk))
    AA[kk,1] = -2./(Dk*Dkp1)
    AA[kk,2] = 2./(Dkp1*(Dk+Dkp1))
    
# Bottom BC - phi(kbtm+1) - at the bottom = 0
# note that dist from z(kbtm) to btm is only half of the grid thickness
  kk   = kbtm
  Dk   = Z_phi[kk-1]-Z_phi[kk]
  Dkp1 = Z_phi[kk]-zbtm
  check_dz(Dk,Dkp1,kk)
   
  AA[kbtm,0] = 2./(Dk*(Dkp1+Dk))
  AA[kbtm,1] = -2./(Dk*Dkp1)
#  AA[kbtm,2] = 0.

  return AA

def check_dz(Dk,Dkp1,kk):
  if Dk < 0.:
    print('ERR: mod_solver - negative dlt Z kk={2} Dk={0}, Dkp1={1}'.\
    format(kk,Dk,Dkp1))
    sys.exit()

  return

def eig_unfoldA(kbtm,AA):
  """
    Form explicitely tridiagonal matrix A from 3 row elements
    For testing
    elements below bottom are ignored
    Compute eigenvalues and eigenvectors w, v
  """
  kdm = AA.shape[0]
  AF = np.zeros((kbtm+1,kbtm+1))
  AF[0,0:2] = AA[0,1:3]
  for ik in range(1,kbtm):
    AF[ik,ik-1:ik+2] = AA[ik-1,:]

  AF[kbtm,kbtm-1:kbtm+1] = AA[kbtm,0:2]
  w, v = np.linalg.eig(AF)

  return w, v


def runmn(B0,ZZ,mnwnd=5):
  """
    Running mean, mnwnd - averaging window should be odd
    uneven depth levels are allowed
    ZZ - interface depths, values B0 are in the mid-points
         i.e. ZZ length = B0 length+1
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
      dz = np.abs(ZZ[ksb]-ZZ[ksb+1])
      dz_sum = dz_sum+dz
      b_sum = b_sum+B0[ksb]*dz
#    print('  ')
    BF[kk] = b_sum/dz_sum

  return BF


def eig_wkb(kbtm,N2z,ZZ,phi,mode=1):
  """
    WKB-approximation for the baroclinic Rossby radii
    outside/inside equatorial band (5-dgr) - Chelton
    numerical integration - need update better quadrature
    ZZ - interface depths, values B0 are in the mid-points
         i.e. ZZ length = N2z length+1
  """
  dZ = np.abs(np.diff(ZZ))
  nlev = N2z.shape[0]
  
  Isum = 0.
  for kk in range(nlev):
    Isum = Isum+np.sqrt(N2z[kk])*dZ[kk]

  Rearth = 6371.e3  # Earth R
  omg = 7.29e-5
  fcor = 2.*omg*np.sin(np.deg2rad(phi))
  betta = (2.*omg*np.cos(np.deg2rad(phi)))/Rearth

  if abs(phi) >= 5.0:
    Lmb = 1./(abs(fcor)*mode*np.pi)*Isum
  else:
    Lmb = np.sqrt(1./(2*betta*mode*np.pi)*Isum)

  return Lmb



def Housholder_Atridiag(A0):
  """
    Factorization A=QR using Householder transformation
    Matrix A is square tridiagonal real not required to be symmetric, 
    only nbnd elements are kept
    in each row (subdiagonal, diagonal, superdiagonal)

    For each step k - have 1 vector = kth column of A
    and submatrix sA=A(k+1:end,k+1:end) that is updated
    (see Bjork p 53 for an algorithm)

    Idea is that H*sA=(H*sA(:,i) H*sA(:,2) ... H*sA(:,end))
    denote a(j)=sA(:,j) is a vector
    Then H*a(j) = (I+alf*u*u')*a(j)= ...
    so do not need to form matrix H
    Keep u(k) vectors for reconstructing H(k) matrices
    H(k)=I+alf*u*u'
    I'm using different version of this formula:
    H(k)=I-2*u*u', where u = w/||w||2
    and w = v+gamma*e1, gamma = ||v||2
    In this case, need to keep only u
    A0(mm+1,:)=0; % extra row for u

  """
  mm = A0.shape[0]
  nbnd = A0.shape[1]

  





