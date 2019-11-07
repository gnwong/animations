import pickle
import numpy as np
import matplotlib as mpl# ; mpl.use('agg')
import matplotlib.pyplot as plt
import imageio

nframes = 80
nframe_per_second = 20
nframepost = 2*nframe_per_second
nsources = 15
np.random.seed(4321)
photon_path_length = 100
photon_path_jump = 10

def genearte_photon_from_geodesic(geodesic,t0):
  dr = geodesic[1,1:] - geodesic[1,:-1]
  dh = geodesic[3,1:] - geodesic[3,:-1]
  dl = np.sqrt( np.power(dr,2.) + np.power(geodesic[1,1:]*dh,2.) )
  if geodesic.shape[1] == 0: return None
  xs = geodesic[1,:] * np.cos(geodesic[3,:])
  ys = geodesic[1,:] * np.sin(geodesic[3,:])
  if ys[0] < -40./7 or ys[0] > 40./7: return None
  ts = np.zeros(geodesic[0,:].shape)
  ts[1:] = np.cumsum(dl)
  return np.stack([ts+t0,xs,ys])

def plot_photon(ax, photon, i):
  ci = np.argmax(photon[0]>=i)
  lasti = ci - photon_path_length
  if lasti < 0: lasti = 0
  ax.plot(photon[1,lasti:ci+1],photon[2,lasti:ci+1],'-',lw=1,c="#dddddd")
  ax.plot(photon[1,ci-1:ci+1],photon[2,ci-1:ci+1],'-w',lw=1.2)

def clean_photons(photons, i):
  newphotons = []
  for photon in photons:
    lasti = np.argmax(photon[0]>=i) - photon_path_length
    if lasti < 0: lasti = 0
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

  for frame in range(nframes+nframepost):

    print("processing frame {0:d} of {1:d} (tracking {2:d})".format(frame+1,nframes+nframepost,len(photons)))
    ofname = 'imgs/schwarzschild_{0:04d}.png'.format(frame)

    plt.close('all')
    fig = plt.figure(figsize=(7,4))
    ax = plt.subplot(1,1,1)

    photons = clean_photons(photons, frame)

    for photon in photons:
      plot_photon(ax, photon, frame)

    if frame in sourceframes:

      for i in range(9):  #should be 9
        while True:
          lgeo = loaded_geodesics[geodcount]
          geodcount += 1
          photon = genearte_photon_from_geodesic(lgeo,frame)
          if photon is not None:
            photons.append( photon )
            break
      x0 = photon[1,0]
      y0 = photon[2,0]
      ax.plot(x0,y0,'*',c='#ffffaa',ms=3)

    ax.set_xlim(-10,10)
    ax.set_ylim(-40./7,40./7)
    ax.set_aspect('equal')
    ax.set_facecolor('black')

    plt.tight_layout()
    plt.box()
    plt.savefig(ofname, bbox_inches='tight', facecolor='black')

    images.append(imageio.imread(ofname))

    if len(photons) == 0:
      break

  imageio.mimsave("schwarzschild_movie.gif", images)

