import numpy as np
import matplotlib as mpl ; mpl.use('agg')
import matplotlib.pyplot as plt
import imageio

nframes = 80
nframe_per_second = 20
nframepost = 2*nframe_per_second
nsources = 15
np.random.seed(4321)
photon_path_length = 3

def generate_photon_flat(x0,y0,h0,t0):
  t = np.linspace(0,1.*(nframes/nframe_per_second),nframes)
  xs = t * np.cos(h0) + x0
  ys = t * np.sin(h0) + y0
  t *= nframe_per_second
  return np.stack([t+t0,xs,ys])

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
    if photon[1,lasti] > 0 and photon[1,lasti] < 1: 
      if photon[2,lasti] > 0 and photon[2,lasti] < 1: 
        newphotons.append(photon)
  return newphotons

if __name__ == "__main__":

  images = []

  sourceframes = [ int(x) for x in np.random.random(nsources-1)*nframes ]
  sourceframes.append(0)

  photons = []

  for frame in range(nframes+nframepost):

    print("processing frame {0:d} of {1:d} (tracking {2:d})".format(frame+1,nframes+nframepost,len(photons)))
    ofname = 'imgs/flat_{0:04d}.png'.format(frame)

    plt.close('all')
    fig = plt.figure(figsize=(7,4))
    ax = plt.subplot(1,1,1)

    photons = clean_photons(photons, frame)

    for photon in photons:
      plot_photon(ax, photon, frame)

    if frame in sourceframes:
      th0 = np.random.random()*2.*np.pi
      x0 = np.random.random()
      y0 = np.random.random() * 4./7
      for th in np.linspace(0,2.*np.pi,11)[1:]:
        photons.append( generate_photon_flat(x0,y0,th+th0,frame) )
      ax.plot(x0,y0,'*',c='#ffffaa',ms=3)

    ax.set_xlim(0,1)
    ax.set_ylim(0,4./7)
    ax.set_aspect('equal')
    ax.set_facecolor('black')

    plt.tight_layout()
    plt.box()
    plt.savefig(ofname, bbox_inches='tight', facecolor='black')

    images.append(imageio.imread(ofname))

    if len(photons) == 0:
      break

  imageio.mimsave("flat_movie.gif", images)

