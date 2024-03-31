import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys


def combine_to_smaller(data):
    combined_data = data[::2, ::2]
    combined_data += data[1::2, ::2]
    combined_data += data[::2, 1::2]
    combined_data += data[1::2, 1::2]
    return combined_data / 4.


if __name__ == "__main__":

    ifnames = sys.argv[1:]

    if len(ifnames) > 1:

        data = None

        for fn in ifnames:
            print(f' - processing {fn}')
            with h5py.File(fn, 'r') as hfp:
                unpol = np.array(hfp['unpol'])
            if data is None:
                data = unpol
            else:
                data += unpol

        fig = plt.figure(figsize=(10, 10), facecolor='w')
        ax1 = plt.subplot(1, 1, 1)
        ax1.imshow(data[::32, ::32].T, cmap='turbo')

        plt.savefig('image.png', dpi=100)

    else:

        n_downsamples = 1  # at least 1 to deal with the final stride of 2

        print(f' - loading {ifnames[0]}')
        ofname = ifnames[0].replace('.h5', f'_{n_downsamples}.png')

        with h5py.File(ifnames[0], 'r') as hfp:
            data = np.array(hfp['unpol'])

        wsize = 409.6
        hsize = wsize

        print(' - downsampling...')
        for _ in range(n_downsamples):
            data = combine_to_smaller(data)
            wsize /= 2.
            hsize /= 2.

        print(' - plotting')
        fig = plt.figure(figsize=(wsize, hsize), facecolor='w')
        ax1 = plt.subplot(1, 1, 1)
        ax1.imshow(data.T, cmap='turbo')
        ax1.set_position([0, 0, 1, 1])

        print(f' - saving {ofname}')
        plt.savefig(ofname, dpi=100)


