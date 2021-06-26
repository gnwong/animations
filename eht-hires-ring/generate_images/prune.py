import os
import sys

if __name__ == "__main__":

    count = 0
    newlines = []

    fname = sys.argv[1]

    fp = open(fname, 'r')
    for line in fp.readlines():
        count += 1
        try:
            ofn = line.split()[-1].replace("--outfile=", "")
            if not os.path.exists(ofn):
                newlines.append(line)
        except:
            pass
    fp.close()

    ofp = open(fname, 'w')
    for newline in newlines:
        ofp.write(newline.strip() + "\n")
    ofp.close()

    print(len(newlines), count)

