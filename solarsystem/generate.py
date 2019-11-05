import numpy as np
import matplotlib as mpl ; mpl.use('agg')
import matplotlib.pyplot as plt
import imageio


## data from https://ssd.jpl.nasa.gov/?horizons where appropriate

# names
planets = [ "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune" ]
# in km
aphelia   = [ 69816900, 108939000, 152100000, 249200000, 816.62e6, 1514.5e6, 3008e6, 4.54e9 ]
perihelia = [ 46001200, 107477000, 147095000, 206700000, 740.52e6, 1352.6e6, 2742e6, 4.46e9 ]
# in years. second is closer to resonances for pretty plotting
orbital_periods = [ 0.240846, 0.615198, 1.000017, 1.88082, 11.862, 29.4571, 84.0205, 164.8 ]
orbital_periods = [ 0.25, 0.6666, 1.0, 1.75, 12., 30., 85., 165. ]
# mean in km
sun_radius = 695510
radii = [ 2439.7, 6051.8, 6371.0, 3389.5, 69911, 58232, 25362, 24622 ]
# in html format
colors = [ "#9B9B9F", "#C2C053", "#6E96FB", "#FC845F", "#F6D395", "#DEC077", "#CFF5F6", "#7191B8" ]


## formatting things
framesperyear = 40.
ext = aphelia[3]*1.1
plot_orbits = True


def find_synch_year(indices,tol=0.1):
  nyears = 0
  while True:
    nyears += 1
    maxdiff = 0.
    for index in indices:
      diff = (nyears / orbital_periods[index]) % 1.
      if diff > 0.5: diff = 1. - diff
      if diff > maxdiff: maxdiff = diff
    print(nyears, maxdiff)
    if maxdiff < tol: break
  return nyears

def plot_sun(ax):
  ax.plot(0,0,'o',c='#FFFF00',ms=20)

def plot_planet(ax, index, time):
  meanradius = (aphelia[index]+perihelia[index])/2.
  x = meanradius * np.cos(2.*np.pi*time/orbital_periods[index])
  y = meanradius * np.sin(2.*np.pi*time/orbital_periods[index])
  ax.plot(x,y,'o',c=colors[index],ms=5)

def plot_orbit(ax, index):
  angles = np.linspace(0,2.*np.pi,101)
  meanradius = (aphelia[index]+perihelia[index])/2.
  x = np.cos(angles) * meanradius
  y = np.sin(angles) * meanradius
  ax.plot(x,y,'-',c=colors[index],lw=0.2)


if __name__ == "__main__":

  images = []

  to_plot = [0,1,2,3]
  sync = find_synch_year(to_plot,tol=0.1)
  nframes = int(framesperyear * sync)

  print("computing for {0:d} frames or {1:f} years".format(nframes, sync))

  for i in range(nframes):

    print("processing frame {0:d} of {1:d}".format(i+1, nframes))

    plt.close('all')
    fig = plt.figure(figsize=(8,8))
    ax = plt.subplot(1,1,1)
    ax.set_facecolor('black')
    ax.set_xlim(-ext,ext)
    ax.set_ylim(-ext,ext)
    ax.set_aspect('equal')
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    time = i/framesperyear
    plot_sun(ax)
    for j in range(8):
      plot_planet(ax, j, time)
      if plot_orbits:
        plot_orbit(ax, j)
    ofname = "imgs/{0:04d}.png".format(i)
    plt.savefig(ofname, bbox_inches="tight")
    images.append(imageio.imread(ofname))

  imageio.mimsave('movie.gif', images)

