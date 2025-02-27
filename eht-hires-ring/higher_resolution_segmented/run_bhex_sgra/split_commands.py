import sys

# 28 cores or 96 cores

preamble = """#!/bin/bash

# run on typhon

#SBATCH --output=imrun-%j.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH --time=72:00:00

pwd
module list
date

"""

if __name__ == "__main__":

    tagname = ""

    ifname = sys.argv[1]
    fp = open(ifname, 'r')
    lines = fp.readlines()
    fp.close()
    commands = [x.strip() for x in lines if len(x) > 0]

    if len(sys.argv) > 2:
        tagname = sys.argv[2]

    ncommands = len(commands)

    print(f"found {ncommands} different commands. splitting into as many submission scripts.")

    n_for_split = 15

    nsplit = int(ncommands / n_for_split)
    nsplit = max(nsplit, 1)

    vfr = ""

    subnames = []

    tag = "im_sgra94_"
    if len(tagname) > 0:
        tag = tagname + "_"


    command_groups = []
    for _ in range(n_for_split):
        command_groups.append([])

    idx = 0
    for cmd in commands:
        command_groups[idx % n_for_split].append(cmd)
        idx += 1

    for v in range(n_for_split):
        subname = '{1:s}{0:04d}.sb'.format(v, tag)
        commands_to_print = command_groups[v]
        if len(commands_to_print) < 1:
            continue
        ofp = open(subname, 'w')
        ofp.write(preamble)
        for cmd in commands_to_print:
            ofp.write(vfr + cmd + "\n")
        ofp.close()
        subnames.append(subname)

    for sn in subnames:
        print('sbatch ' + sn)


