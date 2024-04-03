import os
import sys

from wongutils.grrt import ipole

munits = {
    'iharm_85_10': 7.796696396950134e+24,
    'iharm_85_40': 1.365574591044775e+25,
    'iharm_94_10': 7.434391294578116e+24,
    'iharm_94_40': 1.1736963391247541e+25
}

res = 40960

if __name__ == "__main__":

    # rough estimate at 320^2 resolution is
    # 2560.0 pixels/second

    if len(sys.argv) > 1:
        model = sys.argv[1]
    else:
        model = 'iharm_94_10'

    rhigh = 40
    if model.endswith('_10'):
        rhigh = 10

    base_directory = '/data/gnwong/GRMHD/Ma+0.94_hicadence/dumps/'
    if '85' in model:
        base_directory = '/data/gnwong/GRMHD/Ma+0.85_hicadence/dumps/'

    munit = munits[model]

    dumps = sorted([base_directory + f for f in os.listdir(base_directory) if f.endswith('.h5')])

    dumps = dumps[::len(dumps)//10]

    strides = [0, 32, 0, 32]
    # for stride 32, 32, it took 606.932 seconds to complete a segment, which makes ~170 hours
    # for a single image. producing an image at only 16x16 resolution takes 1/4 of the time, so
    # we can expect a 16x16 image to take ~42 hours to produce

    # gnwong@thalassa is about 600
    # gnwong@nereid is about 550
    # gnwong@ganymede is about 790
    # gnwong@neso is about 600

    stride2 = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]
    stride4 = [0, 4, 8, 12, 16, 20, 24, 28]
    stride16 = [0, 16]

    strides = [0, 8, 0, 8]
    stride2e = [0, 2, 4, 6]
    stride2o = [1, 3, 5, 7]

    for dump in dumps[:-1]:

        if '00004784' not in dump and '00005712' not in dump and '00009424' not in dump and '00011280' not in dump:
            continue

        for si in stride2o:
            for sj in stride2o:

                strides[0] = si
                strides[2] = sj

                stride_string = f'_{strides[0]}_{strides[1]}_{strides[2]}_{strides[3]}_res{res}'

                outfile = '/mnt/lustre/user/gnwong/images/highres/' + model + '/' + os.path.basename(dump).replace('.h5', '') + stride_string + '.h5'
                if not os.path.exists(os.path.dirname(outfile)):
                    os.makedirs(os.path.dirname(outfile))

                if os.path.exists(outfile):
                    continue

                args = ipole.run_ipole(dump, outfile=outfile, rhigh=rhigh, munit=munit, res=res, unpol=True, onlyargs=True, strides=strides)
                print(' '.join(args))

    #ipole.fit_munit(dumps, 0.65, 1.e20, 1.e30, rhigh=40, logname='iharm_85_40.log', verbose=True, res=80)

    #
