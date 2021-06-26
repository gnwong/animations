from subprocess import call
import glob
import os

fnames = glob.glob("/data/bh29-fs1/gnwong2/Ma+0.94_384_hc/dump*h5")
fnames = sorted(fnames)

if __name__ == "__main__":

    ofp = open('run', 'w')

    for fname in fnames[::4]:
        
        ifname = "/data/bh29-fs1/gnwong2/Ma+0.94_384_hc/imgs/"
        ifnameog = ifname + os.path.basename(fname)
        ifname = "imgs/"
        ifname += os.path.basename(fname)

        dfname = "../dumps/" + os.path.basename(fname)

        args  = "./ipole --freqcgs=230.e9 --MBH=6.2e9 --nx=1280 --thetacam=168 --dump=" + dfname
        args += " --M_unit=1.76917e+25 --trat_large=160 --ny=1280 -unpol --outfile2=" + ifname
        args += " --outfile=" + ifnameog
        ofp.write(args + "\n")

    ofp.close()

