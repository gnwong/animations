"""

  Thoughts: run frames 3500 GM/c^3 -> 10k GM/c^3
    = 700 -> 2000

"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
import iharm
from scipy.interpolate import griddata, interp1d
import matplotlib as mpl

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import sys
import os

plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=10)

log = True
res = 512
mad_dir = "/data/bh-bd3/eht/GRMHD/MAD/a+0.94/384x192x192_IHARM/dumps/"
sane_dir = "/data/bh-bd3/eht/GRMHD/SANE/a+0.94/288x128x128_IHARM/dumps/"


def get_interpolated(X, Y, data, xres, yres, xlim, ylim, domax=False):
    Xflat = X.flatten()
    Yflat = Y[:, ::-1].flatten()
    points = np.zeros((len(Xflat),2))
    for i in range(len(Xflat)):
        points[i, 0] = Xflat[i]
        points[i, 1] = Yflat[i]
    if False:
        if domax:
            data[:, 0] = data.max()
            data[:, -1] = data.max()
            fval = data.max()
        else:
            data[:, 0] = data.min()
            data[:, -1] = data.min()
            fval = data.min()
    grid_x, grid_y = np.meshgrid(np.linspace(*xlim, xres), np.linspace(*ylim, yres))
    return griddata(points, data.flatten(), (grid_x, grid_y), method='linear') #, fill_value=fval)

def get_or_make_cache(fname, tlim, res):
    if "MAD" in fname:
        adj = "m"
    elif "SANE" in fname:
        adj = "s"
    else:
        adj = "q"
    adj += os.path.basename(fname).replace("dump", "")
    cfn = "cache/" + adj
    must_make = True
    if os.path.exists(cfn):
        hfp = h5py.File(cfn, 'r')
        if hfp['tlim'][()] == tlim and hfp['res'][()] == res:
            must_make = False
        hfp.close()
    if must_make:
        print(f" - making cache for {fname}")
        R,H,P = iharm.get_RHP(fname, centers=True)
        X = R * np.sin(H) * np.cos(P)
        Y = R * np.sin(H) * np.sin(P)
        Z = R * np.cos(H)
        hfp = h5py.File(fname, 'r')
        rho = np.array(hfp['prims'][:,:,:,0])
        hfp.close()
        N1,N2,N3 = rho.shape
        N2m = N2//2
        elim = [-tlim, tlim]
        interpd = get_interpolated(X[:,N2m], Y[:,N2m], rho[:,N2m], res, res, elim, elim)
        hfp = h5py.File(cfn, 'w')
        hfp['interpd'] = interpd
        hfp['tlim'] = tlim
        hfp['res'] = res
        hfp.close()
    hfp = h5py.File(cfn, 'r')
    interpd = np.array(hfp['interpd'])
    hfp.close()
    return interpd

def load_interp(fname):
    t, v = np.loadtxt(fname).T
    qskip = len(t) // 4 - 1
    return interp1d(t[::qskip], v[::qskip], bounds_error=False, fill_value='extrapolate')

mad_interp = load_interp('norms_mad.dat')
sane_interp = load_interp('norms_sane.dat')

def add_blackhole(ax, bhspin):
    reh = 1.+np.sqrt(1.-bhspin*bhspin)
    bhartist = plt.Circle((0, 0), reh, color='#fff', zorder=100)
    ax.add_artist(bhartist)
    bhartist2 = plt.Circle((0, 0), reh*0.95, color='#000', zorder=101)
    ax.add_artist(bhartist2)

def add_cbar(tax, vmin, vmax, cmap, label, ticks=None, ticklabels=None, labelpad=0):
    axins = inset_axes(tax,
                       width="5%", height="100%", loc='upper left',
                       bbox_to_anchor=(1.05, 0., 1, 1), bbox_transform=tax.transAxes, borderpad=0,
                       )
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cb1 = mpl.colorbar.ColorbarBase(axins, cmap=cmap, norm=norm, orientation='vertical')
    cb1.ax.tick_params(labelsize=11) 
    if ticks is not None:
        cb1.set_ticks(ticks)
        cb1.ax.set_yticklabels(ticklabels)
    cb1.set_label(label, fontsize=14, labelpad=labelpad)

def plot_on(ax, fname, time, tlim=15):
    interpd = get_or_make_cache(fname, tlim, res)

    if "MAD" in fname:
        vmax = mad_interp(time)
    elif "SANE" in fname:
        vmax = sane_interp(time)

    pdata = interpd
    vmin = 0
    if log:
        pdata = np.log10(interpd / vmax)
        vmax = 0
        vmin = -2

    ax.imshow(pdata, extent=[-tlim,tlim,-tlim,tlim], cmap='turbo', vmin=vmin, vmax=vmax) 

    add_blackhole(ax, 0.9375)

def make_frame(frame):

    time = frame * 5

    plt.close('all')
    fig = plt.figure(figsize=(10, 5))

    ax1 = plt.subplot(1, 2, 1)
    ax2 = plt.subplot(1, 2, 2)

    plot_on(ax1, mad_dir + "dump_{0:08d}.h5".format(frame), time)
    plot_on(ax2, sane_dir + "dump_{0:08d}.h5".format(frame), time)

    if log:
        add_cbar(ax2, -2, 0, plt.cm.turbo, r'log$_{10}$ density', labelpad=7, ticks=[-2,-1.5,-1,-.5,0], ticklabels=["-2.0", "-1.5", "-1.0", "-0.5", "0.0"])
    else:
        add_cbar(ax2, 0, 1, plt.cm.turbo, r'density (arbitrary units)', ticks=[0,0.2,0.4,0.6,0.8,1.], ticklabels=[], labelpad=7)

    ax2.set_yticklabels([])
    ax1.tick_params(axis='both', labelsize=14)
    ax2.tick_params(axis='both', labelsize=14)

    time = (frame-800) * 5
    fig.text(0.57, 0.94, f"${time}$ GM/c$^3$", ha='right', va='center', fontsize=20)
    ax1.set_title("MAD", fontsize=16)
    ax2.set_title("SANE", fontsize=18)

    plt.subplots_adjust(wspace=0, hspace=0.01)
    plt.tight_layout(rect=[0, 0.02, 0.93, 0.96])
    plt.savefig(f'imgs/frame_{frame:08d}.png')

if __name__ == "__main__":

    frame = int(sys.argv[1]) 

    print(frame)

    make_frame(frame)



