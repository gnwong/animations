import numpy as np
import matplotlib.pyplot as plt
import geotracer.pytool.geotrace as geotrace
import h5py
import os
import sys
from scipy.optimize import brentq
from scipy.interpolate import interp1d
from colour import Color

def hex_range(c1, c2, count):
    def clamp(x):
        x = int(x * 256)
        return max(0, min(x, 255))
    def hex_from_rgb(r,g,b):
        return "#{0:02x}{1:02x}{2:02x}".format(clamp(r), clamp(g), clamp(b))
    c1tup = np.array(Color(c1).rgb)
    c2tup = np.array(Color(c2).rgb)
    vals = []
    for i in range(count):
        tv = 1. * i / count
        vals.append( hex_from_rgb( *(c1tup*(1.-tv) + c2tup*tv) ) )
    return vals

plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=10)

nobg = "--nobg" in sys.argv

if nobg:
    fade = int(sys.argv[3])
    ogfade = fade
    if fade > 100:
        fade = 200 - fade 

# v1
bgcolor = '#fff'
dark_bgcolor = '#111'
hit_bgcolor = '#900'
inf_bgcolor = '#00f'
color_hit = '#daa'
color_inf = '#aad'

# v2
hit_bgcolor = '#420000'
inf_bgcolor = '#000042'
color_hit = '#ecc'
color_inf = '#cce'

# v3
hit_bgcolor = '#a00'
inf_bgcolor = '#00a'

acolor = '#000'

runtype = "12"
nsamps = 20

face_bgcolor = inf_bgcolor
if nobg:
    c_inf_seq = hex_range('#00f', '#cce', 101)
    color_inf = c_inf_seq[fade]

    c_hit_seq = hex_range('#f00', '#cce', 101)
    color_hit = c_hit_seq[fade]

    c_dbg_seq = hex_range('#000', '#00a', 101)
    c_hbg_seq = hex_range('#000', '#a00', 101)
    dark_bgcolor = c_dbg_seq[fade]
    hit_bgcolor = c_hbg_seq[fade]

    face_bgcolor = dark_bgcolor


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

def get_pruned_r(r, radius_eh): 
    rcut_index = np.argmax(r < radius_eh)
    if rcut_index > 0:
        return 1
        r = r[:rcut_index]
    return r[-1]

def find_edge(geot, rho_left, rho_right, radius_eh, nsamps=10):
    r_left = get_pruned_r(get_geodesic_data(geot, rho_left, 0)[:, 1], radius_eh)
    r_right = get_pruned_r(get_geodesic_data(geot, rho_right, 0)[:, 1], radius_eh)
    rhos_to_return = []
    for i in range(nsamps):
        rho_trial = (rho_left + rho_right)/2.
        r_trial = get_pruned_r(get_geodesic_data(geot, rho_trial, 0)[:, 1], radius_eh)
        if r_left > 10.:
            if r_trial > 10.:
                rho_left = rho_trial
            else:
                rho_right = rho_trial
            rhos_to_return = [rho_right, rho_left]
        else:
            if r_trial > 10.:
                rho_right = rho_trial
            else:
                rho_left = rho_trial
            rhos_to_return = [rho_left, rho_right]
    return rhos_to_return

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

    # bookkeeping for red zone edges
    state = "out"
    last_rho = -100
    edges = []

    for rho in np.linspace(-10, 10, count):

        # get geodesic
        geodBL = get_geodesic_data(geot, rho + rho_offset, 0.)
        t = geodBL[:, 0]
        r = geodBL[:, 1]
        theta = geodBL[:, 2]
        phi = geodBL[:, 3]
     
        nextstate = "out"   

        rcut_index = np.argmax(r < radius_eh)
        if rcut_index > 0:
            t = t[:rcut_index]
            r = r[:rcut_index]
            theta = theta[:rcut_index]
            phi = phi[:rcut_index]
            nextstate = "in"

        if r[-1] < 10.:
            nextstate = "in"

        if nextstate == "in" and state == "out":
            edges.append(find_edge(geot, last_rho + rho_offset, rho + rho_offset, 
                                   radius_eh, nsamps=nsamps))
        elif nextstate == "out" and state == "in":
            edges.append(find_edge(geot, last_rho + rho_offset, rho + rho_offset, 
                                   radius_eh, nsamps=nsamps))
        state = nextstate

        # translate to x,y,z coordinates (units GM/c^2)
        x = np.sin(theta) * r * np.cos(phi) * mass
        y = np.sin(theta) * r * np.sin(phi) * mass
        z = np.cos(theta) * r * mass
        
        xs.append(x)
        ys.append(y)

        last_rho = rho
    
    xedges = []
    yedges = []
    for edge in edges[0] + edges[1]:
        geodBL = get_geodesic_data(geot, edge, 0.)
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
        xedges.append(np.sin(theta) * r * np.cos(phi) * mass)
        yedges.append(np.sin(theta) * r * np.sin(phi) * mass)

    return xs, ys, xedges, yedges

def plot_other_blob(xs, ys, color, ax):
    ax.fill_betweenx(ys, xs, 0, color=color, zorder=8)

def plot_reds_pre(xs1, ys1, xs2, ys2, color, ax):
    ax.fill_betweenx(ys1, xs1, xs2, color=color)

def plot_background_above(xs, ys, color, ax, offset=0):
    if xs[1] < 0.:
        xidx = np.argmax(xs[1:] > 0.) + 1 + offset
        xs = xs[:xidx]
        ys = ys[:xidx]
        ax.fill_between(xs, ys, -20*np.ones_like(ys), color=color)
    else:
        xidx = np.argmax(xs[1:] < 0.) + 1 + offset
        xs = xs[:xidx]
        ys = ys[:xidx]
        ax.fill_between(xs, ys, -20*np.ones_like(ys), color=color)

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

    if not loaded or "--overwritec" in sys.argv:
        xs, ys, xes, yes = get_rays_for_spin_mass(spin, massfrac, count)
        hfp = h5py.File(cachename, 'w')
        hfp['spin'] = spin
        hfp['massfrac'] = massfrac
        hfp['count'] = count
        for i in range(count):
            hfp[f'xs_{i}'] = xs[i]
            hfp[f'ys_{i}'] = ys[i]
        for i in range(4):
          hfp[f'xes_{i}'] = xes[i]
          hfp[f'yes_{i}'] = yes[i]
        hfp.close()

    hfp = h5py.File(cachename, 'r')
    xs = []
    ys = []
    for i in range(count):
        xs.append(np.array(hfp[f'xs_{i}']))
        ys.append(np.array(hfp[f'ys_{i}']))
    xes = [ np.array(hfp['xes_0']), np.array(hfp['xes_1']), 
            np.array(hfp['xes_2']), np.array(hfp['xes_3']) ]
    yes = [ np.array(hfp['yes_0']), np.array(hfp['yes_1']),  
            np.array(hfp['yes_2']), np.array(hfp['yes_3']) ]
    hfp.close()

    if doplot:

        # create figure
        plt.close('all')

        #if nobg:
        #    plt.figure(figsize=(8, 8), facecolor='w')
        #else:
        plt.figure(figsize=(8, 10), facecolor=bgcolor)

        ax1 = plt.subplot(1, 1, 1)

        for i, x in enumerate(xs):
            y = ys[i]
            r = x*x+y*y
            color = color_inf
            #color = '#ccc'
            zorder = 10
            if r[-1] < 20:
                color = color_hit
                zorder = 20
            ax1.plot(y, -x, lw=1, c=color, zorder=zorder)
   
        # left black boundary
        #ax1.plot(yes[1], -xes[1], c='w', lw=1, zorder=7)
        #if not nobg:
        plot_background_above(yes[1], -xes[1], dark_bgcolor, ax1)

        # right black boundary
        #ax1.plot(yes[3], -xes[3], c='w', lw=1, zorder=7)
        #if not nobg:
        plot_background_above(yes[3], -xes[3], dark_bgcolor, ax1)
        if index < 4:
            plot_reds_pre(yes[1], -xes[1], yes[3], -xes[3], '#000', ax1)
            #plot_reds_pre(yes[1], -xes[1], yes[3], -xes[3], dark_bgcolor, ax1)

        orient_up = False
        #ax1.plot(yes[0], -xes[0], c='w', lw=1, zorder=7)
        #ax1.plot(yes[2], -xes[2], c='w', lw=1, zorder=7)
        #if not nobg:
        if index < 230:
            #plot_reds_pre(yes[0], -xs[0], hit_bgcolor, ax1) 
            #plot_reds_pre(yes[2], -xs[2], hit_bgcolor, ax1) 
            #ax1.fill_betweenx(yes[0], -xs[0], -xs[2], color=hit_bgcolor)
            #print(-xs[0].shape, -xs[2].shape)
            plot_reds_pre(yes[0], -xes[0], yes[2], -xes[2], hit_bgcolor, ax1)
        else:
            plot_background_above(yes[0], -xes[0], hit_bgcolor, ax1, offset=3)
            plot_background_above(yes[2], -xes[2], hit_bgcolor, ax1, offset=3)

        if index > 300:
            plot_other_blob(yes[2], -xes[2], hit_bgcolor, ax1)

        # add black hole for style
        radius_eh = 1.+np.sqrt(1.-spin*spin) * 1.05
        bh = plt.Circle((0, 0), radius_eh, color='#aaa', zorder=100)
        ax1.add_artist(bh)

        # format plot
        tlim = 15
        ax1.set_xlim(-tlim, tlim)
        ax1.set_ylim(-tlim, tlim)
        ax1.set_aspect('equal')
        #ax1.set_xlabel(r'$x\ c^2/GM$')
        #ax1.set_ylabel(r'$y\ c^2/GM$')
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])

        ax1.set_facecolor(face_bgcolor)
        #ax1.set_facecolor('#000000')

        ax1.spines['bottom'].set_color(acolor)
        ax1.spines['left'].set_color(acolor)
        ax1.spines['left'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)

        ax1.tick_params(axis='x', colors=acolor)
        ax1.tick_params(axis='y', colors=acolor)
        ax1.tick_params(axis='x', labelsize=22)
        ax1.tick_params(axis='y', labelsize=22)

        ax1.set_yticklabels([])
        ax1.set_yticks([])

        #ax1.set_title(str(spin) + " " + str(massfrac))

        plt.tight_layout(rect=[0.01, 0.01, 1, 1])
        if nobg:
            plt.savefig("imgs/frame_{0:05d}_{1:04d}.png".format(index, ogfade))
        else:
            plt.savefig("imgs/frame_{0:05d}.png".format(index))


if __name__ == "__main__":

    count = 141
    #count = 281
    #count = 561

    skip = 1
    nframes_part1 = 401
    nframes_part2 = 400

    index = 0
    
    tindex = None
    if len(sys.argv) > 1:
        tindex = int(sys.argv[1])

    if "1" in runtype:
        ## first go through massfracs. increase mass linearly
        massfracs = np.logspace(-4, 0, nframes_part1)
        massfracs = np.linspace(1.e-3, 1., nframes_part1)
        for massfrac in massfracs[::skip]:
            if tindex is None or tindex == index:
                print(index)
                plot_frame(index, 0, massfrac, count, doplot=True)
            index += skip

    if "2" in runtype:
        ## then spin it up
        index = nframes_part1
        spins = np.linspace(0, 1, nframes_part2+1)[:-1]
        for spin in spins[::skip]:
            if tindex is None or tindex == index:
                print(index)
                plot_frame(index, spin, 1., count, doplot=True)
            index += skip

