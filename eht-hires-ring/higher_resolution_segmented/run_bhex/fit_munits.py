import os
import sys

from wongutils.grrt import ipole

if __name__ == "__main__":

    base_directory = '/data/gnwong/GRMHD/Ma+0.94_hicadence/dumps/'
    base_directory = '/data/gnwong/GRMHD/Ma+0.85_hicadence/dumps/'

    if len(sys.argv) > 1:
        spin = float(sys.argv[1])

    if spin == 0:
        base_directory = '/scratch/lustre/gnwong/GRMHD/eht_v1/MAD/a0/384x192x192_IHARM/dumps/'
        logname = 'bhex_iharm_00_40.log'
    else:
        base_directory = '/scratch/lustre/gnwong/GRMHD/eht_v1/MAD/a+0.94/384x192x192_IHARM/dumps/'
        logname = 'bhex_iharm_94_40.log'

    dumps = sorted([base_directory + f for f in os.listdir(base_directory) if f.endswith('.h5')])

    dumps = dumps[::len(dumps)//10]
    print(dumps)

    ipole.fit_munit(dumps, 0.65, 1.e20, 1.e30, rhigh=40, thetacam=60, logname=logname, verbose=True, res=80)

    #
