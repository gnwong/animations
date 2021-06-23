import numpy as np
import matplotlib.pyplot as plt
import geotracer.pytool.geotrace as geotrace

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

def get_rays_for_spin_mass(spin, mass):
    rcam = 1.e8
    thetacam = 90.
    geot = geotrace.geotrace(bhspin=spin, mbh_msun=6.2e9, rmax_geo=100./mass)
    geot.set_camera_location(tcam=0., rcam=rcam, thetacam=thetacam)
    geot.set_camera_properties(DX=80./mass, Dsource_cm=16.9e6*3.085678e18)
    
    xs = []
    ys = []
    
    radius_eh = (1.+np.sqrt(1.-spin*spin)) / mass
    
    for rho in np.linspace(-10, 10, 51):

        # get geodesic
        geodBL = get_geodesic_data(geot, rho, 0.)
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


def plot_frame(index, spin, massfrac):

    xs, ys = get_rays_for_spin_mass(spin, massfrac)    

    # create figure
    plt.figure(figsize=(6, 6), facecolor='w')
    ax1 = plt.subplot(1, 1, 1)

    for i, x in enumerate(xs):
        y = ys[i]
        ax1.plot(y, -x)
        
    # add black hole for style
    radius_eh = 1.+np.sqrt(1.-spin*spin)
    bh = plt.Circle((0, 0), radius_eh, color='#000000', zorder=10)
    ax1.add_artist(bh)

    # format plot
    tlim = 20
    ax1.set_xlim(-tlim, tlim)
    ax1.set_ylim(-tlim, tlim)
    ax1.set_aspect('equal')
    ax1.set_xlabel('x c^2/GM')
    ax1.set_ylabel('y c^2/GM')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    ax1.set_title(str(spin)  + " " + str(massfrac))

    plt.savefig("imgs/frame_{0:05d}.png".format(index))


if __name__ == "__main__":

    index = 0

    ## first go through massfracs. increase mass linearly
    for massfrac in np.logspace(-4, 0):
        print(index)
        plot_frame(index, 0, massfrac)
        index += 1

    ## then spin it up



