import h5py
import sys
import numpy as np
import tqdm

if __name__ == "__main__":

    flux = "mad"
    if "s_" in sys.argv[1]:
        flux = "sane"

    fp = open("norms_" + flux + ".dat", 'w')

    for fname in tqdm.tqdm(sorted(sys.argv[1:])):

        hfp = h5py.File(fname, 'r')
        interpd = np.array(hfp['interpd'])
        hfp.close()

        time = int(fname.replace(".h5", "").split('_')[-1]) * 5

        fp.write("{0:g} {1:g}\n".format(time, interpd.max()))

    fp.close()
          


