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
  with open("data/geodesics.pickle",'rb') as handle:
    loaded_geodesics = pickle.load(handle)

  photons = []
  geodcount = 0

  allphotons = []

  for frame in range(nframes+nframepost):

    print("processing frame {0:d} of {1:d} (tracking {2:d})".format(frame+1,nframes+nframepost,len(photons)))
    ofname = 'imgs/schwarzschild_{0:04d}.png'.format(frame)

    if run_with_plotting:
      plt.close('all')
      fig = plt.figure(figsize=(7,4))
      ax = plt.subplot(1,1,1)

    photons = clean_photons(photons, frame)

    for photon in photons:
      if run_with_plotting:
        plot_photon(ax, photon, frame)

    if frame in sourceframes:
      while True:
        lgeo = loaded_geodesics[geodcount]
        geodcount += 1
        photon = genearte_photon_from_geodesic(lgeo,frame)
        if photon is not None:
          photons.append( photon )
          allphotons.append( photon )
          break
      x0 = photon[1,0]
      y0 = photon[2,0]
      while True:
        lgeo = loaded_geodesics[geodcount]
        geodcount += 1
        if lgeo.shape[1] == 0:
          geodcount -= 1
          break
        tx = lgeo[1,0]*np.cos(lgeo[3,0])
        ty = lgeo[1,0]*np.sin(lgeo[3,0])
        if tx==x0 and ty==y0:
          photon = genearte_photon_from_geodesic(lgeo,frame,force=True)
          photons.append( photon )
          allphotons.append( photon )
        else:
          geodcount -= 1
          break
      if run_with_plotting:
        ax.plot(x0,y0,'*',c='#ffffaa',ms=3)

    if run_with_plotting:
      ax.set_xlim(-10,10)
      ax.set_ylim(-40./7,40./7)
      ax.set_aspect('equal')
      ax.set_facecolor('black')

      thetas = np.linspace(0,2.*np.pi,101)
      ax.plot(2.*np.cos(thetas),2.*np.sin(thetas),'-',color="#222222")
      blackhole = plt.Circle((0., 0.), 2., color='#555555')
      ax.add_artist(blackhole)

      plt.tight_layout()
      plt.box()
      plt.savefig(ofname, bbox_inches='tight', facecolor='black')

      images.append(imageio.imread(ofname))

    if len(photons) == 0:
      break

  if run_with_plotting:
    imageio.mimsave("schwarzschild_movie.gif", images)

  allphotons = []
  for i in range(len(loaded_geodesics)):
    lgeo = loaded_geodesics[i]
    photon = genearte_photon_from_geodesic(lgeo,frame,force=True)
    if photon is not None:
      allphotons.append( photon )

  plt.close('all')
  fig = plt.figure(figsize=(7,4))
  ax = plt.subplot(1,1,1)
  for photon in allphotons:
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
  plt.savefig('schwarzschild_allphotons.png', bbox_inches='tight', facecolor='black')




