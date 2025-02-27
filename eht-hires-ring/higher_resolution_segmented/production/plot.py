import numpy as np
import matplotlib.pyplot as plt
import ehtplot.color

import os

fnames = [
  "dump_00004784_segA_2.npy",
  "dump_00005712_segA_2.npy",
  "dump_00009424_segA_2.npy",
  "dump_00011280_segA_2.npy"
]

bfname = "/mnt/lustre/user/gnwong/images/highres/iharm_94_10/combined/"
bfname = "/Users/gnwong/r/codes/gnwong/animations/eht-hires-ring/higher_resolution_segmented/plot/remote/"


if __name__ == "__main__":

    if not os.path.exists('cache.npy'):
        
        imdata = None
        wsize = None
        hsize = None

        print(' - loading files to make cache')
        for fname in fnames:
            data = np.load(bfname + fname, allow_pickle=True).item()
            if imdata is None:
                imdata = data['data']
                wsize = data['wsize']
                hsize = data['hsize']
            else:
                imdata += data['data']

        print(' - saving cache')
        data_to_save = dict(data=imdata, wsize=wsize, hsize=hsize)
        np.save('cache.npy', np.array(data_to_save, dtype=object))

    print(' - loading data')
    ddict = np.load('cache.npy', allow_pickle=True).item()
    imdata = ddict['data']
    wsize = ddict['wsize']
    hsize = ddict['hsize']

    print(' - setting image weights')
    """
    pdata = np.power(imdata, 0.5) * 4.
    pdata += np.power(imdata, 0.3) * 10.
    pdata += np.power(imdata, 2) * 15
    vmax = np.max(pdata) * np.sqrt(1.2)"""
    pdata = np.power(imdata, 0.8) * 5.
    pdata += np.power(imdata, 0.32) * 12.
    pdata += np.power(imdata, 3) * 20
    vmax = np.max(pdata) * np.sqrt(1.1)

    print(' - plotting')
    fig = plt.figure(figsize=(wsize, hsize), facecolor='w')
    ax1 = plt.subplot(1, 1, 1)
    ax1.imshow(pdata[:, ::-1], origin='lower', cmap='afmhot', vmin=0, vmax=vmax)
    ax1.set_position([0, 0, 1, 1])

    print(' - saving')
    plt.savefig('cache_afmhot.png', dpi=100)


