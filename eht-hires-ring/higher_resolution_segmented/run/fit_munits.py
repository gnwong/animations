import os

from wongutils.grrt import ipole

if __name__ == "__main__":

    base_directory = '/data/gnwong/GRMHD/Ma+0.94_hicadence/dumps/'
    base_directory = '/data/gnwong/GRMHD/Ma+0.85_hicadence/dumps/'

    dumps = sorted([base_directory + f for f in os.listdir(base_directory) if f.endswith('.h5')])

    dumps = dumps[::len(dumps)//5]
    print(dumps)

    ipole.fit_munit(dumps, 0.65, 1.e20, 1.e30, rhigh=40, logname='iharm_85_40.log', verbose=True, res=80)

    #
