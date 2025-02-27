
basedir = '/mnt/lustre/user/gnwong/images/highres/bhex_iharm_94_30/'

fnames = [1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900]

for f in fnames:

    globp = basedir + f'dump_{f:08d}_*_res20480.h5'
    opath = basedir + f'combined/dump_{f:08d}.h5'

    c = 'python combine.py ' + globp + ' '  + opath

    print(c)

