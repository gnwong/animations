import numpy as np
import h5py
import matplotlib.pyplot as plt
import iharm

import os

mad_dir = "/data/bh-bd3/eht/GRMHD/MAD/a+0.94/384x192x192_IHARM/dumps/"
sane_dir = "/data/bh-bd3/eht/GRMHD/SANE/a+0.94/288x128x128_IHARM/dumps/"


def make_frame(frame):
    
    fig = plt.figure(figsize=(10, 5))


    plt.savefig('imgs/frame_{0:08d}.png')

if __name__ == "__main__":

    

    make_frame(0)



