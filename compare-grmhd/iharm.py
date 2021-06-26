"""

  hooks for loading iharm v3.5 data


  load_iharm(harm_dump_name, gcov, gcon): returns (rho, UU, U, B, ucon, ucov, bcon, bcov)

    If gcov/gcon are None, then the four-vector quantities are not
    computed from the primitives (ucon,ucov,bcon,bcov are None).


  get_RHP(harm_dump_name): returns (R,H,P)
    
    ... set

  get gcov(harm_dump_name): returns gcov[N1,N2,N3,4,4]




  get_gfname(dump_fname) : 
    returns string of conventional location for grid file associated
    with dump at dump_fname

  get_gcov(grid_fname) :
    returns N1xN2xN3x4x4 numpy array of gcov stored in grid file

  get_gcon
  get_alpha
  get_gdet
  get_gdet_zone
  get_r
  get_h
  get_p
  load_iharm_header
  load_iharm_Jsq
  load_thermo
  load_iharm


"""

import os
import h5py
import numpy as np

__version__ = 1.2

def get_RHP(fname,centers=False):
  """ load 3-dimensional R, H, P for FMKS grid (with or without center alignment) """
  R = None ; H = None ; P = None
  hfp = h5py.File(fname,'r')
  hdrname = 'header'
  if 'fluid_header' in hfp.keys(): hdrname = 'fluid_header'
  metric = hfp[hdrname]['metric'][()].decode('utf-8').lower() 
  if metric == "mks":
    raise NotImplementedError
  elif metric in ["fmks", "mmks"]:
    N1 = hfp[hdrname]['n1'][()]
    N2 = hfp[hdrname]['n2'][()]
    N3 = hfp[hdrname]['n3'][()]
    Rin = None ; Rout = None
    try: Rin = hfp[hdrname]['geom'][metric]['r_in'][()]
    except: Rin = hfp[hdrname]['geom'][metric]['Rin'][()]
    try: Rout = hfp[hdrname]['geom'][metric]['r_out'][()]
    except: Rout = hfp[hdrname]['geom'][metric]['Rout'][()]
    mks_smooth = hfp[hdrname]['geom'][metric]['mks_smooth'][()]
    hslope = hfp[hdrname]['geom'][metric]['hslope'][()]
    poly_alpha = hfp[hdrname]['geom'][metric]['poly_alpha'][()]
    poly_xt = hfp[hdrname]['geom'][metric]['poly_xt'][()]
    poly_norm = 0.5*np.pi*1./(1. + 1./(poly_alpha + 1.)*1./np.power(poly_xt, poly_alpha))
    x1 = np.linspace(np.log(Rin),np.log(Rout),N1+1)
    x2 = np.linspace(0,1,N2+1)
    x3 = np.linspace(0,2.*np.pi,N3+1)
    if centers:
      x1 = ( x1[1:] + x1[:-1] ) / 2.
      x2 = ( x2[1:] + x2[:-1] ) / 2. 
      x3 = ( x3[1:] + x3[:-1] ) / 2.
    r = np.exp(x1)
    hg = np.pi*x2 + (1.-hslope)*np.sin(2.*np.pi*x2)/2.
    y = 2.*x2 - 1.
    hj = poly_norm*y*(1.+np.power(y/poly_xt,poly_alpha)/(poly_alpha+1.)) + 0.5*np.pi
    X1,HG,X3 = np.meshgrid(x1,hg,x3)
    R,HJ,P = np.meshgrid(r,hj,x3)
    H = HG + np.exp(mks_smooth*(np.log(Rin) - X1))*(HJ - HG)
    pass
  else:
    raise NotImplementedError
  hfp.close()
  return R.transpose((1,0,2)),H.transpose((1,0,2)),P.transpose((1,0,2))

def gcov2gcon(gcov, is2d=None):
  """ convert gcov to gcon. makes some assumptions about gcov being 2d """
  if is2d is None:
    n3 = gcov.shape[2]
    if np.allclose(gcov[:, :, 0], gcov[:, :, n3//3]) and \
       np.allclose(gcov[:, :, 0], gcov[:, :, n3//2]):
       is2d = True
  if is2d:
    gcov2d = gcov[:, :, 0, :, :]
    gcon2d = np.linalg.inv(gcov2d)
    gcon = np.zeros_like(gcov)
    gcon[:,:,:,:,:] = gcon2d[:,:,None,:,:]
    return gcon
  return np.linalg.inv(gcov)

def gcov2gdet(gcov):
  """ convert gcov -> gdet """
  return np.sqrt(-np.linalg.det(gcov))

def get_metric_name(fname):

  hfp = h5py.File(fname,'r')
  hdrname = 'header'
  if 'fluid_header' in hfp.keys(): hdrname = 'fluid_header'
  try:
    metric = hfp[hdrname]['metric'][()].decode('utf-8').lower() 
  except:
    metric = hfp[hdrname]['metric'][()].lower()
  hfp.close()

  return metric

def load_gcov(fname):

  metric_type = get_metric_name(fname)

  if metric_type in ["mmks", "fmks"]:
    return gcov_fmks(fname)
  elif metric_type == "eks":
    return gcov_eks(fname)

  raise NotImplementedError
  return None

def gcov_eks(fname):
  # provides gcov with KS components at R
  # this function leverages the fact that phi is a cyclic coordinate

  hfp = h5py.File(fname,'r')

  hdrname = 'header'
  if 'fluid_header' in hfp.keys(): hdrname = 'fluid_header'
  metric = get_metric_name(fname)
  if metric != "eks": 
    raise NotImplementedError
    return None

  N1 = hfp[hdrname]['n1'][()]
  N2 = hfp[hdrname]['n2'][()]
  N3 = hfp[hdrname]['n3'][()]

  a = hfp[hdrname]['geom'][metric]['a'][()]
  startx1 = hfp[hdrname]['geom']['startx1'][()]
  startx2 = hfp[hdrname]['geom']['startx2'][()]
  startx3 = hfp[hdrname]['geom']['startx3'][()]
  dx1 = hfp[hdrname]['geom']['dx1'][()]
  dx2 = hfp[hdrname]['geom']['dx2'][()]
  dx3 = hfp[hdrname]['geom']['dx3'][()]

  x1 = startx1 + dx1 * np.arange(N1)
  x2 = startx2 + dx2 * np.arange(N2)
  X1,X2 = np.meshgrid(x1,x2)

  r = np.exp(x1)
  h = np.pi*x2
  R,H = np.meshgrid(r,h)

  hfp.close()

  R = R.T
  H = H.T

  gcov = np.zeros((R.shape[0],R.shape[1],4,4))
  cth = np.cos(H)
  sth = np.sin(H)
  s2 = sth*sth
  rho2 = R*R + a*a*cth*cth
  gcov[:,:,0,0] = (-1. + 2. * R / rho2)
  gcov[:,:,0,1] = (2. * R / rho2)
  gcov[:,:,0,3] = (-2. * a * R * s2 / rho2)
  gcov[:,:,1,0] = gcov[:,:,0,1]
  gcov[:,:,1,1] =  (1. + 2. * R / rho2)
  gcov[:,:,1,3] =  (-a * s2 * (1. + 2. * R / rho2))
  gcov[:,:,2,2] = rho2
  gcov[:,:,3,0] = gcov[:,:,0,3]
  gcov[:,:,3,1] = gcov[:,:,1,3]
  gcov[:,:,3,3] = s2 * (rho2 + a*a * s2 * (1. + 2. * R / rho2))

  # get transformation to eks
  dxdX = np.zeros((R.shape[0],R.shape[1],4,4))
  dxdX[:,:,0,0] = 1.
  dxdX[:,:,1,1] = R
  dxdX[:,:,2,2] = 1.
  dxdX[:,:,3,3] = 1.

  # transform and extend along x3 dimension
  gcov2d = np.einsum('abki,abkj->abij',dxdX,np.einsum('ablj,abkl->abkj',dxdX,gcov))
  gcov = np.zeros((N1,N2,N3,4,4))
  gcov[:,:,:,:,:] = gcov2d[:,:,None,:,:]

  return gcov

def gcov_fmks(fname):
  # provides gcov with fmks components for the specified harm dump file
  # this function leverages the fact that phi is a cyclic coordinate

  # get metric/coordinate parameters & dimensions
  hfp = h5py.File(fname,'r')
  hdrname = 'header'
  if 'fluid_header' in hfp.keys(): hdrname = 'fluid_header'
  metric = get_metric_name(fname)
  if metric != "mmks" and metric != "fmks": 
    raise NotImplementedError
    return None
  N1 = hfp[hdrname]['n1'][()]
  N2 = hfp[hdrname]['n2'][()]
  N3 = hfp[hdrname]['n3'][()]
  Rin = None ; Rout = None
  try: Rin = hfp[hdrname]['geom'][metric]['r_in'][()]
  except: Rin = hfp[hdrname]['geom'][metric]['Rin'][()]
  try: Rout = hfp[hdrname]['geom'][metric]['r_out'][()]
  except: Rout = hfp[hdrname]['geom'][metric]['Rout'][()]
  mks_smooth = hfp[hdrname]['geom'][metric]['mks_smooth'][()]
  hslope = hfp[hdrname]['geom'][metric]['hslope'][()]
  a = hfp[hdrname]['geom'][metric]['a'][()]
  poly_alpha = hfp[hdrname]['geom'][metric]['poly_alpha'][()]
  poly_xt = hfp[hdrname]['geom'][metric]['poly_xt'][()]
  poly_norm = 0.5*np.pi*1./(1. + 1./(poly_alpha + 1.)*1./np.power(poly_xt, poly_alpha))
  x1 = np.linspace(np.log(Rin),np.log(Rout),N1+1)
  x2 = np.linspace(0,1,N2+1)
  x1 = ( x1[1:] + x1[:-1] ) / 2.
  x2 = ( x2[1:] + x2[:-1] ) / 2. 
  X1,X2 = np.meshgrid(x1,x2)
  r = np.exp(x1)
  hg = np.pi*x2 + (1.-hslope)*np.sin(2.*np.pi*x2)/2.
  y = 2.*x2 - 1.
  hj = poly_norm*y*(1.+np.power(y/poly_xt,poly_alpha)/(poly_alpha+1.)) + 0.5*np.pi
  X1,HG = np.meshgrid(x1,hg)
  R,HJ = np.meshgrid(r,hj)
  H = HG + np.exp(mks_smooth*(np.log(Rin) - X1))*(HJ - HG)
  hfp.close()
  R = R.T ; H = H.T  ;  X1 = X1.T ; X2 = X2.T

  # calculate for ks components in ks coordinates at r,h
  gcov = np.zeros((R.shape[0],R.shape[1],4,4))
  cth = np.cos(H)
  sth = np.sin(H)
  s2 = sth*sth
  rho2 = R*R + a*a*cth*cth
  gcov[:,:,0,0] = (-1. + 2. * R / rho2)
  gcov[:,:,0,1] = (2. * R / rho2)
  gcov[:,:,0,3] = (-2. * a * R * s2 / rho2)
  gcov[:,:,1,0] = gcov[:,:,0,1]
  gcov[:,:,1,1] =  (1. + 2. * R / rho2)
  gcov[:,:,1,3] =  (-a * s2 * (1. + 2. * R / rho2))
  gcov[:,:,2,2] = rho2
  gcov[:,:,3,0] = gcov[:,:,0,3]
  gcov[:,:,3,1] = gcov[:,:,1,3]
  gcov[:,:,3,3] = s2 * (rho2 + a*a * s2 * (1. + 2. * R / rho2))

  # get transformation to fmks
  dxdX = np.zeros((R.shape[0],R.shape[1],4,4))
  dxdX[:,:,0,0] = 1.
  dxdX[:,:,1,1] = R
  dxdX[:,:,3,3] = 1.
  dxdX[:,:,2,1] = - np.exp(mks_smooth*(np.log(Rin)-X1))*mks_smooth*( \
    np.pi/2. - \
    np.pi*X2 + \
    poly_norm*(2.*X2-1.)*(1.+(np.power((-1.+2.*X2)/poly_xt,poly_alpha))/(1.+poly_alpha)) - \
    1./2.*(1.-hslope)*np.sin(2.*np.pi*X2))
  dxdX[:,:,2,2] = np.pi + (1. - hslope)*np.pi*np.cos(2.*np.pi*X2) + \
    np.exp(mks_smooth*(np.log(Rin)-X1))*( \
      -np.pi + \
      2.*poly_norm*(1.+np.power((2.*X2-1.)/poly_xt,poly_alpha)/(poly_alpha+1.)) + \
      (2.*poly_alpha*poly_norm*(2.*X2-1.)*np.power((2.*X2-1.)/poly_xt,poly_alpha-1.))/((1.+poly_alpha)*poly_xt) - \
      (1.-hslope)*np.pi*np.cos(2.*np.pi*X2))

  # extend along x3 dimension
  gcov2d = np.einsum('abki,abkj->abij',dxdX,np.einsum('ablj,abkl->abkj',dxdX,gcov))
  gcov = np.zeros((N1,N2,N3,4,4))
  gcov[:,:,:,:,:] = gcov2d[:,:,None,:,:]
  return gcov

def get_gdet_zone(grid_name):
  """ UNVERIFIED . """
  return h5py.File(grid_name,'r')["gdet_zone"][:,:,:]

def load_iharm_header(harm_dump_name):
  """ UNVERIFIED . does what you expect """

  header = {}

  hfp = h5py.File(harm_dump_name,"r")
  header["t"] = hfp["t"][()]
  header["n1"] = hfp["header"]["n1"][()]
  header["n2"] = hfp["header"]["n2"][()]
  header["n3"] = hfp["header"]["n3"][()]
  header["geom"] = {}
  header["geom"]["startx1"] = hfp["header"]["geom"]["startx1"][()]
  header["geom"]["startx2"] = hfp["header"]["geom"]["startx2"][()]
  header["geom"]["startx3"] = hfp["header"]["geom"]["startx3"][()]
  header["geom"]["dx1"] = hfp["header"]["geom"]["dx1"][()]
  header["geom"]["dx2"] = hfp["header"]["geom"]["dx2"][()]
  header["geom"]["dx3"] = hfp["header"]["geom"]["dx3"][()]
  header["metric"] = str(hfp["header"]["metric"][()]).replace("b'","").replace("'","")
  try: header["gridfile"] = hfp["header"]["gridfile"][()]
  except: pass

  if header["metric"] == "MKS3":
    header["geom"]["mks3"] = {}
    header["geom"]["mks3"]["H0"] = hfp["header"]["geom"]["mks3"]["H0"][()]
    header["geom"]["mks3"]["R0"] = hfp["header"]["geom"]["mks3"]["R0"][()]
    header["geom"]["mks3"]["MY1"] = hfp["header"]["geom"]["mks3"]["MY1"][()]
    header["geom"]["mks3"]["MY2"] = hfp["header"]["geom"]["mks3"]["MY2"][()]
    header["geom"]["mks3"]["MP0"] = hfp["header"]["geom"]["mks3"]["MP0"][()]
    header["geom"]["mks3"]["a"] = hfp["header"]["geom"]["mks3"]["a"][()]

  hfp.close()

  return header

def load_iharm_Jsq(harm_dump_name, gcov, gcon):
  """ UNVERIFIED . does what you expect """

  hfp = h5py.File(harm_dump_name,"r")
  N1 = hfp["header"]["n1"][()]
  N2 = hfp["header"]["n2"][()]
  N3 = hfp["header"]["n3"][()]

  jcon = hfp["jcon"][:,:,:,:]
  jcov = np.einsum('abcij,abci->abcj',gcov,jcon)

  U    = hfp["prims"][:,:,:,2:5]
  alpha = 1. / np.sqrt(-gcon[:,:,:,0,0])
  gamma = np.sqrt(1. + np.einsum('abci,abci->abc',np.einsum('abcij,abci->abcj',gcov[:,:,:,1:,1:],U),U))
  ucon = np.zeros((N1,N2,N3,4))
  ucon[:,:,:,1:] = U - gamma[:,:,:,None] * alpha[:,:,:,None] * gcon[:,:,:,0,1:]
  ucon[:,:,:,0] = gamma / alpha

  data = {}
  data["jcon"] = jcon
  data["jsq"]  = np.einsum('abci,abci->abc',jcon,jcov)
  data["Jsq"]  = data["jsq"] + np.einsum('abci,abci->abc',ucon,jcov)

  return data

def load_thermo(harm_dump_name):
  """ UNVERIFIED . """
  
  hfp = h5py.File(harm_dump_name,"r")
  KEL   = hfp["prims"][:,:,:,8]
  KTOT  = hfp["prims"][:,:,:,9]
  hfp.close()

  return (KTOT, KEL)

def get_thetae(UU, rho, bsq, gam, rhigh=40., game=4./3., gamp=5./3, rlow=1.):
  """ compute Thetae from dump-loaded variables given rlow/rhigh kwargs """
  mp = 1.67262171e-24
  me = 9.1093826e-28
  bsq = np.copy(bsq)
  bsq[bsq == 0.] = UU[bsq == 0.]/1.e10
  beta = (gam - 1.) * UU / bsq / 0.5
  R = (beta*beta * rhigh + rlow) / (1. + beta*beta)
  ThetaeUnit = mp/me * (game-1.) * (gamp-1.) / ( (gamp-1.) + (game-1.)*R )
  Thetae = ThetaeUnit * UU / rho
  Thetae[Thetae < 1.e-3] = 1.e-3
  return Thetae

def load_iharm_all(harm_dump_name):
  """ Loads gcov, gcon for you. """
  gcov = load_gcov(harm_dump_name)
  gcon = gcov2gcon(gcov)
  return load_iharm(harm_dump_name, gcov, gcon), gcov, gcon

def get_mdot_code(rho, ucon, gdet, dx2, dx3):
  """Returns dMact in code units."""
  dMacts = gdet * rho * ucon[:, :, :, 1]
  return dMacts[:21].sum() * dx2 * dx3 / 21.

def get_Phi_code(B1, gdet, dx2, dx3, norm=1.):
  """Returns un-normalized Phi in code units."""
  Phi_comp = 0.5 * norm * np.fabs(B1[:21]) * gdet[:21]
  return Phi_comp.sum() * dx2 * dx3 / 21.

def load_iharm(harm_dump_name, gcov, gcon, compute_b=True):
  """ does what you expect """

  hfp = h5py.File(harm_dump_name,"r")
  N1 = hfp["header"]["n1"][()]
  N2 = hfp["header"]["n2"][()]
  N3 = hfp["header"]["n3"][()]

  rho  = hfp["prims"][:,:,:,0]
  UU   = hfp["prims"][:,:,:,1]
  U    = hfp["prims"][:,:,:,2:5]
  U1   = hfp["prims"][:,:,:,2]
  U2   = hfp["prims"][:,:,:,3]
  U3   = hfp["prims"][:,:,:,4]
  B    = hfp["prims"][:,:,:,5:8]
  B1   = hfp["prims"][:,:,:,5]
  B2   = hfp["prims"][:,:,:,6]
  B3   = hfp["prims"][:,:,:,7]

  hfp.close()

  if gcov is not None:

    alpha = 1. / np.sqrt(-gcon[:,:,:,0,0])
    gamma = np.sqrt(1. + np.einsum('abci,abci->abc',np.einsum('abcij,abci->abcj',gcov[:,:,:,1:,1:],U),U))
    ucon = np.zeros((N1,N2,N3,4))
    ucon[:,:,:,1:] = U - gamma[:,:,:,None] * alpha[:,:,:,None] * gcon[:,:,:,0,1:]
    ucon[:,:,:,0] = gamma / alpha
    ucov = np.einsum('abcij,abci->abcj',gcov,ucon)

    if compute_b:

      bcon = np.zeros((N1,N2,N3,4))
      bcon[:,:,:,0] = np.einsum('abci,abci->abc',B,ucov[:,:,:,1:])
      bcon[:,:,:,1:] = ( B + ucon[:,:,:,1:] * bcon[:,:,:,0,None] ) / ucon[:,:,:,0,None]
      bcov = np.einsum('abcij,abci->abcj',gcov,bcon)

    else:
      
      bcon = None
      bcov = None

  else:

    ucon = None
    ucov = None
    bcon = None
    bcov = None

  return (rho, UU, U, B, ucon, ucov, bcon, bcov)
  

