import numpy as np
import matplotlib as mpl ; mpl.use('agg')
import matplotlib.pyplot as plt
import imageio

plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=10)

nframes_per_second = 10
nframes = 100

def plot_ring(ax, x0,y0, r, c):
  phi = np.linspace(0,2.*np.pi,101)
  x = np.cos(phi) * r + x0
  y = np.sin(phi) * r + y0
  ax.plot(x,y,c,lw=1)

def add_signal_for(x0,y0,t0,x1,y1,signal):
  dt = 1./nframes_per_second
  d = np.sqrt( (x0-x1)**2 + (y0-y1)**2 )
  tidx = int(d/dt) + t0 + 20
  if tidx < len(signal): signal[ tidx ] = 1.
  return signal

emitters = [ [-2,2,"#ff0000"], [-2,0,"#00ff00"], [-2,-2,"#0033ff"] ]

if __name__ == "__main__":

  images = []

  pulses = [ 5 ]

  if False:
    i = nframes

    plt.close('all')
    fig = plt.figure(figsize=(12,8))
    ax = plt.subplot(1,2,1)
    ax1 = plt.subplot(2,2,2)
    ax2 = plt.subplot(2,2,4)

    # plot emitters
    for em in emitters: 
      ax.plot(em[0],em[1],'o',ms=10,c=em[2])

    # plot antenna
    ax.plot(1,1.5,'ow',ms=5)
    ax.plot(1,-0.5,'ow',ms=5)

    # plot rings
    r1 = np.sqrt(9.+1.5*1.5)
    r2 = np.sqrt(9.+0.5*0.5)
    plot_ring(ax, 1,1.5, r1, "#00ff00")
    plot_ring(ax, 1,-0.5, r2, "#00ff00")

    # format plot
    ax.set_facecolor('black')
    ax.set_aspect('equal')
    ax.set_xlim(-3,2)
    ax.set_ylim(-3,3)

    # plot signals for top receiver
    for em in emitters:
      signal = np.zeros(101)
      for pulse in pulses:
        signal = add_signal_for(em[0],em[1],pulse,1,1.5, signal)
      ax1.plot(signal[:i][::-1], ':', c=em[2])

    # plot signals for bottom receiver
    for em in emitters:
      signal = np.zeros(101)
      for pulse in pulses:
        signal = add_signal_for(em[0],em[1],pulse,1,-0.5, signal)
      ax2.plot(signal[:i][::-1], ':', c=em[2])

    # add some axes
    ax1.axhline(y=0,c='#ffffff',lw=1)
    #ax1.axvline(x=0,c='#ffffff',lw=1)
    ax2.axhline(y=0,c='#ffffff',lw=1)
    #ax2.axvline(x=0,c='#ffffff',lw=1)

    # format other plots
    ax1.set_ylim(-0.2,1.2)
    ax2.set_ylim(-0.2,1.2)
    ax1.set_xlim(0,40)
    ax2.set_xlim(0,40)
    ax1.set_facecolor('black')
    ax2.set_facecolor('black')
    ax1.set_ylabel('signal (top)',c='#ffffff',fontsize=20)
    ax2.set_xlabel('time',c='#ffffff',fontsize=20)
    ax2.set_ylabel('signal (bottom)',c='#ffffff',fontsize=20)

    ofname = "imgs/reversed.png".format(i)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(ofname, bbox_inches='tight', facecolor='black')


  for i in range(17,75):
  #for i in [ 72, 73, 74, 75 ]:

    print("processing frame {0:d} of {1:d}".format(i+1,nframes))

    plt.close('all')
    fig = plt.figure(figsize=(12,8))
    ax = plt.subplot(1,2,1)
    ax1 = plt.subplot(2,2,2)
    ax2 = plt.subplot(2,2,4)

    # plot emitters
    for em in emitters: 
      ax.plot(em[0],em[1],'o',ms=10,c=em[2])

    # plot antenna
    ax.plot(1,1.5,'ow',ms=5)
    ax.plot(1,-0.5,'ow',ms=5)

    # plot red rings
    em = emitters[0]
    signali = 43.6
    r = (i-signali)/nframes_per_second
    if r > 0:
      plot_ring(ax, 1.,1.5, r, em[2])
    signali = 35
    r = (i-signali)/nframes_per_second
    if r > 0:
      plot_ring(ax, 1.,-0.5, r, em[2])
    # plot green rings
    em = emitters[1]
    signali = 41
    r = (i-signali)/nframes_per_second
    if r > 0:
      plot_ring(ax, 1.,1.5, r, em[2])
    signali = 44
    r = (i-signali)/nframes_per_second
    if r > 0:
      plot_ring(ax, 1.,-0.5, r, em[2])
    # plot blue rings
    em = emitters[2]
    signali = 28.4
    r = (i-signali)/nframes_per_second
    if r > 0:
      plot_ring(ax, 1.,1.5, r, em[2])
    signali = 41
    r = (i-signali)/nframes_per_second
    if r > 0:
      plot_ring(ax, 1.,-0.5, r, em[2])


    # format plot
    ax.set_facecolor('black')
    ax.set_aspect('equal')
    ax.set_xlim(-3,2)
    ax.set_ylim(-3,3)

    # plot signals for top receiver
    for em in emitters:
      signal = np.zeros(101)
      for pulse in pulses:
        signal = add_signal_for(em[0],em[1],pulse,1,1.5, signal)
      ax1.plot(signal[:nframes-i][::-1], ':', c=em[2])

    # plot signals for bottom receiver
    for em in emitters:
      signal = np.zeros(101)
      for pulse in pulses:
        signal = add_signal_for(em[0],em[1],pulse,1,-0.5, signal)
      ax2.plot(signal[:nframes-i][::-1], ':', c=em[2])

    # add some axes
    ax1.axhline(y=0,c='#ffffff',lw=1)
    ax2.axhline(y=0,c='#ffffff',lw=1)

    # format other plots
    ax1.set_ylim(-0.2,1.2)
    ax2.set_ylim(-0.2,1.2)
    ax1.set_xlim(0,40)
    ax2.set_xlim(0,40)
    ax1.set_facecolor('black')
    ax2.set_facecolor('black')
    ax1.set_ylabel('signal (top)',c='#ffffff',fontsize=20)
    ax2.set_xlabel('time',c='#ffffff',fontsize=20)
    ax2.set_ylabel('signal (bottom)',c='#ffffff',fontsize=20)

    ofname = "imgs/reversed_{0:04d}.png".format(i)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(ofname, bbox_inches='tight', facecolor='black')

    images.append(imageio.imread(ofname))

  imageio.mimsave("movie_reversed.gif", images)

