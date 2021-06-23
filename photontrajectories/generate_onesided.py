import pickle
import numpy as np
import matplotlib as mpl# ; mpl.use('agg')
import matplotlib.pyplot as plt
import imageio

nframes = 320
nframe_per_second = 20
nframepost = 2*nframe_per_second
nsources = 80
np.random.seed(4321)
photon_path_length = 100
photon_path_length_dl = 3
photon_path_jump = 10
run_with_plotting = True

def get_ci_lasti(photon,i):
  ci = np.argmax(photon[0]>=i)
  lasti = np.argmax(photon[0]>=i-photon_path_length_dl)
  if ci == 0:
    ci = photon.shape[1]
    true_last_time = i
    lasti = np.argmax( photon[0] > true_last_time - photon_path_length_dl )
    if true_last_time - photon_path_length_dl > np.max(photon[0]):
      lasti = photon.shape[1]
  if lasti < 0: lasti = 0
  return ci,lasti

def genearte_photon_from_geodesic(geodesic,t0,force=False):
  dr = geodesic[1,1:] - geodesic[1,:-1]
  dh = geodesic[3,1:] - geodesic[3,:-1]
  dl = np.sqrt( np.power(dr,2.) + np.power(geodesic[1,1:]*dh,2.) )
  if geodesic.shape[1] == 0 and not force: return None
  xs = geodesic[1,:] * np.cos(geodesic[3,:])
  ys = geodesic[1,:] * np.sin(geodesic[3,:])
  if ys[0] < -40./7 or ys[0] > 40./7 and not force: return None
  ts = np.zeros(geodesic[0,:].shape)
  ts[1:] = np.cumsum(dl)
  return np.stack([ts+t0,xs,ys])

def plot_photon(ax, photon, i):
  ci, lasti = get_ci_lasti(photon, i)
  ax.plot(photon[1,lasti:ci+1],photon[2,lasti:ci+1],'-',lw=1,c="#dddddd")
  ax.plot(photon[1,ci-1:ci+1],photon[2,ci-1:ci+1],'-w',lw=1.2)

def clean_photons(photons, i):
  newphotons = []
  for photon in photons:
    ci,lasti = get_ci_lasti(photon, i)
    if lasti < photon.shape[1]: 
      if photon[1,lasti] > -10. and photon[1,lasti] < 10.: 
        if photon[2,lasti] > -40./7 and photon[2,lasti] < 40./7: 
          newphotons.append(photon)
  return newphotons

if __name__ == "__main__":

  images = []

  sourceframes = [ int(x) for x in np.random.random(nsources-1)*nframes ]
  sourceframes.append(0)

  print("loading geodesics")
  with open("data/geodesics-fixed.pickle",'rb') as handle:
    loaded_geodesics = pickle.load(handle)

  allphotons = []
  for i in range(len(loaded_geodesics)):
    lgeo = loaded_geodesics[i]
    photon = genearte_photon_from_geodesic(lgeo,0,force=True)
    if photon is not None:
      allphotons.append( photon )

  plt.close('all')
  fig = plt.figure(figsize=(7,4))
  ax = plt.subplot(1,1,1)
  for photon in allphotons[len(allphotons)//2:]:
    ax.plot(photon[1],photon[2],'-',lw=1,c="#cccccc")
  thetas = np.linspace(0,2.*np.pi,101)
  ax.plot(2.*np.cos(thetas),2.*np.sin(thetas),'-',color="#222222")
  blackhole = plt.Circle((0., 0.), 2., color='#555555')
  ax.add_artist(blackhole)
  ax.set_xlim(-10,10)
  ax.set_ylim(-40./7,40./7)
  ax.set_aspect('equal')
  ax.set_facecolor('black')
  plt.tight_layout()
  plt.box()
  plt.savefig('schwarzschild_cameraphotons.png', bbox_inches='tight', facecolor='black')




