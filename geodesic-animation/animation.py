import numpy as np
import matplotlib.pyplot as plt
import geotracer.pytool.geotrace as geotrace
import h5py
import os
from scipy.optimize import brentq

plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=10)

color_hit = '#f00'
color_inf = '#00f'
bgcolor = '#000'

def get_geodesic_data(geot, rho, varphi):
    nx = 40
    ny = 40
    i = rho * np.cos(varphi) + nx//2 - 0.5
    j = rho * np.sin(varphi) + ny//2 - 0.5
    ii = int(i)
    jj = int(j)
    di = i - ii
    dj = j - jj
    X = np.zeros(4)
    Kcon = np.copy(X)
    Kcov = np.copy(X)
    geot.init_XK(ii,jj,X,Kcon,Kcov,nx,ny,di=di,dj=dj)
    geod_dl = geot.push_photon_full(X, Kcon, n=500, dl=True)
    geod = np.array(geod_dl[:,:4])
    geodBL = np.copy(geod)
    geot.bl_coord(geod,geodBL)
    geod_dl[:,:4] = geodBL
    return geod_dl[:,:4]

def find_zero(geot):
    rho_left = -0.1
    rho_right = 0.1
    phi_left = get_geodesic_data(geot, rho_left, 0)[:, 3]
    phi_right = get_geodesic_data(geot, rho_right, 0)[:, 3]
    root = brentq(lambda x: get_geodesic_data(geot, x, 0)[-1, 3], -0.1, 0.1)
    return root

def get_rays_for_spin_mass(spin, mass, count):
    rcam = 1.e8
    thetacam = 90.
    geot = geotrace.geotrace(bhspin=spin, mbh_msun=6.2e9, rmax_geo=100./mass)
    geot.set_camera_location(tcam=0., rcam=rcam, thetacam=thetacam)
    geot.set_camera_properties(DX=80./mass, Dsource_cm=16.9e6*3.085678e18)
    
    xs = []
    ys = []
    
    radius_eh = (1.+np.sqrt(1.-spin*spin)) / mass

    rho_offset = 0.009999999999998045
    if mass < 1.:
        rho_offset = find_zero(geot)
    
    for rho in np.linspace(-10, 10, count):

        # get geodesic
        geodBL = get_geodesic_data(geot, rho + rho_offset, 0.)
        t = geodBL[:, 0]
        r = geodBL[:, 1]
        theta = geodBL[:, 2]
        phi = geodBL[:, 3]
        
        rcut_index = np.argmax(r < radius_eh)
        if rcut_index > 0:
            t = t[:rcut_index]
            r = r[:rcut_index]
            theta = theta[:rcut_index]
            phi = phi[:rcut_index]

        # translate to x,y,z coordinates (units GM/c^2)
        x = np.sin(theta) * r * np.cos(phi) * mass
        y = np.sin(theta) * r * np.sin(phi) * mass
        z = np.cos(theta) * r * mass
        
        xs.append(x)
        ys.append(y)
        
    return xs, ys


def plot_frame(index, spin, massfrac, count, doplot=True):

    cachename = "cache/frame_{0:05d}.h5".format(index)
    loaded = False

    try:
        if os.path.exists(cachename):
            hfp = h5py.File(cachename, 'r')
            if hfp['spin'][()] == spin and hfp['massfrac'][()] == massfrac and hfp['count'][()] == count:
                loaded = True
            hfp.close()
    except:
        pass

    if not loaded:
        xs, ys = get_rays_for_spin_mass(spin, massfrac, count)
        hfp = h5py.File(cachename, 'w')
        hfp['spin'] = spin
        hfp['massfrac'] = massfrac
        hfp['count'] = count
        for i in range(count):
            hfp[f'xs_{i}'] = xs[i]
            hfp[f'ys_{i}'] = ys[i]
        hfp.close()

    hfp = h5py.File(cachename, 'r')
    xs = []
    ys = []
    for i in range(count):
        xs.append(np.array(hfp[f'xs_{i}']))
        ys.append(np.array(hfp[f'ys_{i}']))
    hfp.close()

    if doplot:

        # create figure
        plt.close('all')
        plt.figure(figsize=(8, 8), facecolor=bgcolor)
        ax1 = plt.subplot(1, 1, 1)

        for i, x in enumerate(xs):
            y = ys[i]
            r = x*x+y*y
            color = color_inf
            if r[-1] < 20:
                color = color_hit
            ax1.plot(y, -x, lw=1, c=color)
            
        # add black hole for style
        radius_eh = 1.+np.sqrt(1.-spin*spin) * 1.05
        bh = plt.Circle((0, 0), radius_eh, color='#aaa', zorder=10)
        ax1.add_artist(bh)

        # format plot
        tlim = 15
        ax1.set_xlim(-tlim, tlim)
        ax1.set_ylim(-tlim, tlim)
        ax1.set_aspect('equal')
        ax1.set_xlabel(r'$x\ c^2/GM$')
        ax1.set_ylabel(r'$y\ c^2/GM$')
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])

        ax1.set_facecolor(bgcolor)

        acolor = 'w'
        ax1.spines['bottom'].set_color(acolor)
        ax1.spines['left'].set_color(acolor)
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)

        ax1.tick_params(axis='x', colors=acolor)
        ax1.tick_params(axis='y', colors=acolor)
        ax1.tick_params(axis='x', labelsize=22)
        ax1.tick_params(axis='y', labelsize=22)

        ax1.set_title(str(spin) + " " + str(massfrac))

        plt.tight_layout(rect=[0, 0, 1, 1])
        plt.savefig("imgs/frame_{0:05d}.png".format(index))


if __name__ == "__main__":

    count = 71

    skip = 4
    nframes_part1 = 101
    nframes_part2 = 100

    index = 0

    if True:
        ## first go through massfracs. increase mass linearly
        massfracs = np.logspace(-4, 0, nframes_part1)
        massfracs = np.linspace(1.e-3, 1., nframes_part1)
        for massfrac in massfracs[::skip]:
            print(index)
            plot_frame(index, 0, massfrac, count, doplot=True)
            index += skip

    if True:
        ## then spin it up
        index = nframes_part1
        spins = np.linspace(0, 0.99, 100)
        for spin in spins[::skip]:
            print(index)
            plot_frame(index, spin, 1., count, doplot=True)
            index += skip

