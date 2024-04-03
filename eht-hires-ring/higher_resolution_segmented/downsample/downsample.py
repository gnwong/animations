import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys
import os


def combine_to_smaller(data):
    combined_data = data[::2, ::2]
    combined_data += data[1::2, ::2]
    combined_data += data[::2, 1::2]
    combined_data += data[1::2, 1::2]
    return combined_data / 4.


if __name__ == "__main__":

    ifname = sys.argv[1]

    n_downsamples = 2  # at least 1 to deal with the final stride of 2

    print(f' - loading {ifname}')
    ofname = ifname.replace('.h5', '') + f'_{n_downsamples}.npy'
    if os.path.exists(ofname):
        print(f' - {ofname} already exists. skipping.')
        exit()

    with h5py.File(ifname, 'r') as hfp:
        data = np.array(hfp['unpol'])

    wsize = 409.6
    hsize = wsize

    print(' - downsampling...')
    for _ in range(n_downsamples):
        data = combine_to_smaller(data)
        wsize /= 2.
        hsize /= 2.

    print(f' - saving to {ofname}')
    data_to_save = dict(data=data, wsize=wsize, hsize=hsize)
    np.save(ofname, np.array(data_to_save, dtype=object))

