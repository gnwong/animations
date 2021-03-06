"""

  Makes hero high resolution movie
  2020.11.22 gnw

$ python plot.py path/to/images/*h5

$ ffmpeg -framerate 20 -i image%*.png -s:v 1280x720 -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p ../out.mp4

"""

import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=10)
import numpy as np
import h5py
import sys
import os

from scipy.ndimage.filters import gaussian_filter

import ehtplot.color


# manually chosen
VMAX = 0.00085 
t0 = 4841.01856589023
t0 = 5357.02
show_time = 'M'  # 'days' or 'M'
overwrite = False
overwrite = "--overwrite" in sys.argv
blur = 0.5
blur = None


if __name__ == "__main__":

  if False:
    # ensure that the 1d plots align with how "extent" works in imshow
    ax = plt.subplot(1, 1, 1)
    ax.imshow(np.random.random((4,4)), extent=[-2,2,-2,2])
    txs = np.linspace(-2, 2, 5)
    txs = (txs[1:] + txs[:-1])/2.
    print(txs)
    plt.savefig('alignment.png')
    exit()

  for fname in sys.argv[1:]:

    if fname[-3:] != ".h5": continue
    print("plotting {0:s}".format(fname))

    ofname = os.path.basename(fname.replace(".h5", ".png"))

    if blur is not None and blur > 0:
        ofname = "imgs/blur_" + ofname
    else:
        ofname = "imgs/" + ofname

    if os.path.exists(ofname) and not overwrite:
      print(" - file exists: skipping")
      continue

    # load
    hfp = h5py.File(fname,'r')    
    dx = hfp['header']['camera']['dx'][()]
    t = hfp['header']['t'][()]
    tunit = hfp['header']['units']['T_unit'][()]
    dsource = hfp['header']['dsource'][()]
    lunit = hfp['header']['units']['L_unit'][()]
    fov_muas = dx / dsource * lunit * 2.06265e11
    scale = hfp['header']['scale'][()]
    evpa_0 = 'W'
    if 'evpa_0' in hfp['header']:
      evpa_0 = hfp['header']['evpa_0'][()]
    unpol = np.copy(hfp['unpol']).transpose((1,0))
    #imagep = np.copy(hfp['pol']).transpose((1,0,2))
    #I = imagep[:,:,0]
    #Q = imagep[:,:,1]
    #U = imagep[:,:,2]
    #V = imagep[:,:,3]
    hfp.close()

    I = unpol

    if blur is not None and blur > 0.:
        print("blurring to " + str(blur))
        I = gaussian_filter(I, blur, truncate=3)

    # get critical curve
    a = bhspin = 0.9375
    inc = 168./180. * np.pi
    r = np.linspace(2.208018820801882, 2.814799181479918, 100000)
    lamb = - r*r*r + 3.*r*r - a*a*(r+1.)
    lamb /= a * (r-1.)
    eta = r*r*r/a/a/(r-1.)/(r-1.)
    eta *= 4.*a*a - r*(r-3.)**2
    eta = np.sqrt(eta)
    alpha = - lamb / np.sin(inc)
    beta = np.sqrt(eta*eta + a*a*np.cos(inc)**2 - lamb*lamb / np.tan(inc)/np.tan(inc))
    xs = alpha[np.isfinite(beta)]
    ys = beta[np.isfinite(beta)]

    plt.close('all')

    facecolor = ehtplot.color.ctab.get_cmap('afmhot_10us')(0)

    if True:

        figside = 10
        fig = plt.figure(figsize=(figside, figside), facecolor=facecolor)
        ax = plt.subplot(1, 1, 1)

        nx, ny = I.shape
        dds = dx / nx

        xoff = 6.*dds
        yoff = dds/2.
        ax.imshow(I, cmap='afmhot_10us', vmin=0., vmax=VMAX, extent=[-dx/2+xoff, dx/2+xoff, -dx/2+yoff, dx/2+yoff], interpolation='none')

        # print(I.max(), t, tunit*t)

        time_in_M = r"{0:d} GM/c$^3$".format(int(t-t0))
        time_in_days = "{0:.1f} days".format((t-t0)*tunit / 3600/24.)

        if show_time == 'M':
            ax.text(13, 12, time_in_M, color='w', fontsize=32, horizontalalignment='right', bbox=dict(edgecolor='w', alpha=0.))
        elif show_time == 'days':
            ax.text(13, 12, time_in_days, color='w', fontsize=32, horizontalalignment='right', bbox=dict(edgecolor='w', alpha=0.))

        if False:
            ax.plot(xs,  ys, '-g', lw=0.5)
            ax.plot(xs, -ys, '-g', lw=0.5)

        if False:
            ax.set_xlim(-5, -4)
            ax.set_ylim(-0.5, 0.5)
        elif False:
            ax.set_xlim(0, 1)
            ax.set_ylim(4.25, 5.25)

        ax.set_xlim(-15, 15)
        ax.set_ylim(-15, 15)

        # make the axes white
        acolor = 'w'
        ax.spines['bottom'].set_color(acolor)
        ax.spines['left'].set_color(acolor)
        ax.spines['top'].set_visible(False) 
        ax.spines['right'].set_visible(False)

        ax.tick_params(axis='x', colors=acolor)
        ax.tick_params(axis='y', colors=acolor)
        ax.tick_params(axis='x', labelsize=22)
        ax.tick_params(axis='y', labelsize=22)

    else:
        # make some 1d plots to check the alignment

        fig = plt.figure(figsize=(8, 4))
        ax1 = plt.subplot(1, 2, 1)
        ax2 = plt.subplot(1, 2, 2)

        nx, ny = I.shape
        dds = dx / nx

        tx = np.linspace(-dx/2, dx/2, nx+1)
        tx = (tx[1:]+tx[:-1])/2. - dds/2.
        ty = np.linspace(-dx/2, dx/2, nx+1)
        ty = (ty[1:]+ty[:-1])/2. + 6.*dds

        ax1.plot(tx, I[:, ny//2], 'r')
        ax2.plot(tx, I[:, ny//2], 'r')
        ax1.plot(ty, I[nx//2, :], 'b')
        ax2.plot(ty, I[nx//2, :], 'b')

        ax1.axvline(xs[0], c='b', ls=':')
        ax1.axvline(xs[-1], c='b', ls=':')
        ax1.axvline(ys.max(), c='r', ls=':')
        ax1.axvline(-ys.max(), c='r', ls=':')

        ax2.axvline(xs[0], c='b', ls=':')
        ax2.axvline(xs[-1], c='b', ls=':')
        ax2.axvline(ys.max(), c='r', ls=':')
        ax2.axvline(-ys.max(), c='r', ls=':')

        Islice = I[nx//2, :]
        Isliceargmax = np.argmax(Islice == Islice.max())

        ax1.set_xlim(-5, -4)
        ax2.set_xlim(4.5, 5.5)

    plt.tight_layout(rect=[0, 0, 1, 1])
    fig.savefig(ofname)

