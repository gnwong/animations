"""



"""


import numpy as np
import h5py
from subprocess import call

import sys
import os

if __name__ == "__main__":

    # for 32 stride
    # /mnt/lustre/flash/gnwong/imgs/highres/iharm_94_10

    # python combine.py /mnt/lustre/user/gnwong/images/highres/subrings/dump_00005712_*_res40960_n000.h5 /mnt/lustre/user/gnwong/images/highres/subrings/combined/dump_00005712_res40960_n000.h5
    # python combine.py /mnt/lustre/user/gnwong/images/highres/subrings/dump_00005712_*_res40960_n001.h5 /mnt/lustre/user/gnwong/images/highres/subrings/combined/dump_00005712_res40960_n001.h5
    # python combine.py /mnt/lustre/user/gnwong/images/highres/subrings/dump_00005712_*_res40960_n002.h5 /mnt/lustre/user/gnwong/images/highres/subrings/combined/dump_00005712_res40960_n002.h5

    # python combine.py /mnt/lustre/user/gnwong/images/highres/subrings/dump_00008496_*_res40960_n000.h5 /mnt/lustre/user/gnwong/images/highres/subrings/combined/dump_00008496_res40960_n000.h5
    # python combine.py /mnt/lustre/user/gnwong/images/highres/subrings/dump_00008496_*_res40960_n001.h5 /mnt/lustre/user/gnwong/images/highres/subrings/combined/dump_00008496_res40960_n001.h5
    # python combine.py /mnt/lustre/user/gnwong/images/highres/subrings/dump_00008496_*_res40960_n002.h5 /mnt/lustre/user/gnwong/images/highres/subrings/combined/dump_00008496_res40960_n002.h5

    # python combine.py /mnt/lustre/user/gnwong/images/highres/subrings/dump_00010352_*_res40960_n000.h5 /mnt/lustre/user/gnwong/images/highres/subrings/combined/dump_00010352_res40960_n000.h5
    # python combine.py /mnt/lustre/user/gnwong/images/highres/subrings/dump_00010352_*_res40960_n001.h5 /mnt/lustre/user/gnwong/images/highres/subrings/combined/dump_00010352_res40960_n001.h5
    # python combine.py /mnt/lustre/user/gnwong/images/highres/subrings/dump_00010352_*_res40960_n002.h5 /mnt/lustre/user/gnwong/images/highres/subrings/combined/dump_00010352_res40960_n002.h5
    
    # python combine.py /mnt/lustre/user/gnwong/images/highres/subrings/dump_00003856_*_res40960_n000.h5 /mnt/lustre/user/gnwong/images/highres/subrings/combined/dump_00003856_res40960_n000.h5
    # python combine.py /mnt/lustre/user/gnwong/images/highres/subrings/dump_00003856_*_res40960_n001.h5 /mnt/lustre/user/gnwong/images/highres/subrings/combined/dump_00003856_res40960_n001.h5
    # python combine.py /mnt/lustre/user/gnwong/images/highres/subrings/dump_00003856_*_res40960_n002.h5 /mnt/lustre/user/gnwong/images/highres/subrings/combined/dump_00003856_res40960_n002.h5

    ifnames = sys.argv[1:-1]
    ofname = sys.argv[-1]

    if os.path.exists(ofname):
        print("File exists, skipping")
        sys.exit(0)

    unpol = None

    for idx, ifn in enumerate(ifnames):
        print(f' - loading {ifn} ({idx+1} of {len(ifnames)})')
        with h5py.File(ifn, 'r') as hfp:
            if unpol is None:
                unpol = np.array(hfp['unpol'])
            else:
                unpol += np.array(hfp['unpol'])

    print(f' - saving {ofname}')
    call(['cp', ifnames[0], ofname])
    with h5py.File(ofname, 'r+') as ohfp:
        del ohfp['unpol']
        ohfp['unpol'] = unpol
        



