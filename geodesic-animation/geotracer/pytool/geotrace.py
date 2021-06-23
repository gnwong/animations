import ctypes
from scipy.integrate import solve_ivp
import numpy as np
import mpmath
import os

class geotrace:
  """

    __init__(bhspin=0.)

    init_model(bhspin)

    init_XK(i,j,Xcam,fovx,fovy,X,Kcon,nx,ny)

  """

  # constants
  G = 6.6742e-8
  Msun = 1.989e33
  CL = 2.99792458e10

  # scaling
  Lunit = 0.

  dll = None
  Xcam = np.zeros(4)     # internal representation
  xcam_ks = np.zeros(4)  # ks
  fovx = 1.
  fovy = 1.
  bhspin = 0.

  def __init__(self, bhspin=0., mbh_msun=6.2e9, eps=None, rmax_geo=100.):
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "libgeotrace.so")
    self._load_dll(path=path)
    self.bhspin = bhspin
    self.dll.init_model(bhspin, rmax_geo)
    self.Lunit = self.G * mbh_msun * self.Msun / self.CL / self.CL
    if eps is not None:
      self.dll.set_eps(eps)

  def _load_dll(self,path="./libgeotrace.so"):
    self.dll = ctypes.CDLL(path)

    self.dll.calc.argtypes = [ np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                               ctypes.c_int,
                               np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS") ]
   
    self.dll.init_model.argtypes = [ ctypes.c_double, ctypes.c_double ]

    self.dll.init_XK.argtypes = [ ctypes.c_int,
                                  ctypes.c_int,
                                  ctypes.c_double,
                                  ctypes.c_double,
                                  np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                  ctypes.c_double,
                                  ctypes.c_double,
                                  np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                  np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                  np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                  ctypes.c_int,
                                  ctypes.c_int ]

    self.dll.theta_rootfind.argtypes = [ np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS") ] 
    self.dll.theta_rootfind.restype = ctypes.c_double

    self.dll.valid_geodesic.argtypes = [ np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                         np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS") ] 
    self.dll.valid_geodesic.restype = ctypes.c_int

    self.dll.bl_coord_vec.argtypes = [ np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                       np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS") ] 

    self.dll.bl_coord_vec_many.argtypes = [ np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                            np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                            ctypes.c_int ]

    self.dll.push_photon_hidden.argtypes = [ np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                             np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS") ] 

    self.dll.push_photon_hidden_n.argtypes = [ np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                               np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                               ctypes.c_int ] 
    self.dll.push_photon_hidden_n.restype = ctypes.c_int

    self.dll.push_photon_hidden_n_dl.argtypes = [ np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                                  np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                                  ctypes.c_int ] 
    self.dll.push_photon_hidden_n_dl.restype = ctypes.c_int

    self.dll.set_eps.argtypes = [ ctypes.c_double ]

    self.dll.set_dl.argtypes = [ np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS") ] 

    self.dll.set_dl_native_coordinates.argtypes = [ np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                                    np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS") ] 

    self.dll.set_gcov.argtypes = [ np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                   np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS") ]

  ### These functions all provide interfaces into the geotrace c library.
  def set_gcov(self, X, gcov):
    self.dll.set_gcov(X, gcov)

  def set_camera_location(self, tcam=0., rcam=0., thetacam=0., phicam=0.):
    """Arguments given in degrees."""
    self.xcam_ks = np.array([ tcam, rcam, thetacam/180.*np.pi, phicam/180.*np.pi ])
    self.Xcam[0] = self.xcam_ks[0]
    self.Xcam[1] = np.log(self.xcam_ks[1])
    self.Xcam[3] = self.xcam_ks[3]/180.*np.pi
    self.Xcam[2] = self.dll.theta_rootfind(self.xcam_ks)

  def set_camera_properties(self, DX=None, fov_muas=None, Dsource_cm=16.9e6*3.085678e18):
    if fov_muas is not None:
      DX = fov_muas * Dsource_cm / self.Lunit / 2.06265e11;
    self.fovx = DX / self.xcam_ks[1]
    self.fovy = self.fovx

  def init_XK(self,i,j,X,Kcon,Kcov,nx,ny,di=0.,dj=0.,fovx=None,fovy=None,Xcam=None):
    if Xcam is None: Xcam = self.Xcam
    if fovx is None: fovx = self.fovx
    if fovy is None: fovy = self.fovy
    self.dll.init_XK(i,j,di,dj,Xcam,fovx,fovy,X,Kcon,Kcov,nx,ny)

  def theta_rootfind(self,xcam_ks):
    return self.dll.theta_rootfind(xcam_ks)

  def calc(self,a,n,b):
    self.dll.calc(a,n,b)

  def valid_geodesic(self, X, Kcon):
    return self.dll.valid_geodesic(X, Kcon)

  def bl_coord(self, X, Xbl):
    if type(X) == np.ndarray:
      self.dll.bl_coord_vec_many(np.ascontiguousarray(X), np.ascontiguousarray(Xbl), X.shape[0])
    else:
      self.dll.bl_coord_vec(X, Xbl)

  def push_photon(self, X, Kcon):
    self.dll.push_photon_hidden(X, Kcon)

  def push_photon_n(self, X, Kcon, n, dl=False):
    if dl: return self.dll.push_photon_hidden_n_dl(X, Kcon, n)
    return self.dll.push_photon_hidden_n(X, Kcon, n)

  def push_photon_full(self, X, Kcon, n=100, dl=False):
    if not dl: return self._push_photon_full_nodl(X, Kcon, n=n)
    Xs = np.zeros((1,5))
    Xs[0,:4] = X
    while True:
      Xarr = np.zeros((n,5))
      Xarr[0,:] = Xs[-1]
      k = self.push_photon_n(Xarr, Kcon, n, dl=True)
      Xs = np.concatenate((Xs,Xarr[:k,:]))
      if k < n: break
    return Xs

  def _push_photon_full_nodl(self, X, Kcon, n=100):
    Xs = np.zeros((1,4))
    Xs[0,:4] = X
    while True:
      Xarr = np.zeros((n,4))
      Xarr[0,:] = Xs[-1]
      k = self.dll.push_photon_hidden_n(Xarr, Kcon, n)
      Xs = np.concatenate((Xs,Xarr[:k,:]))
      if k < n: break
    return Xs

  def set_dl(self, X, dl, ks=False):
    if ks:
      self.dll.set_dl(X, dl, len(dl))
    else:
      self.dll.set_dl_native_coordinates(X, dl, len(dl))
