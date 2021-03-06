#!/usr/bin/env python

# Ryan A. Melnyk
# Bioinformatics Scientist - Pivot Bio
# ryan@pivotbio.com

import os
import argparse
import random
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(description='''
    Script that removes redundancy and/or generates jackknife resampling
    (i.e. without replacement).
    ''')
    parser.add_argument('align_file', type=str, help='''path to the alignment
    file for sampling''')
    parser.add_argument('prefix', type=str, help='prefix for jackknife files')
    parser.add_argument('--size', type=int, help='''length of jackknife
    (default 10000)''')
    parser.add_argument('--num', type=int, help='''number of jackknifes
    (default 1)''')
    parser.add_argument('--remove_redundant', action='store_true', help='''use
    if you wish to remove sites that are non-informative (i.e. the same in all
    sequences)''')
    parser.add_argument('--remove_gapped', type=float, help='''use if you wish
    to remove gapped sites - enter max proportion sites that can be gapped.''')
    return parser.parse_args()


def parse(align_file):
    seqdata = {}
    print("Reading alignment into memory...")
    for seq in SeqIO.parse(open(align_file, 'r'), 'fasta'):
        seqdata[str(seq.id)] = str(seq.seq)
    print("Done!")
    return seqdata


def remove_redundant(seqdata):
    first = seqdata.keys()[0]
    newdata = {s: [] for s in seqdata.keys()}
    length = len(seqdata[first])
    print(length, "residues to scan...")
    count = 0
    uniqcount = 0
    for i in range(0, length):
        res = seqdata[first][i]
        matching = True
        for s in seqdata:
            if seqdata[s][i] != res:
                matching = False
        if not matching:
            for s in seqdata:
                newdata[s].append(seqdata[s][i])
            uniqcount += 1
        count += 1
        if count % 100000 == 0:
            if count < 1e6:
                print("{}00K residues parsed...".format(
                    str(int(count/100000)))
                )
            else:
                print("{:.1f}M residues parsed...".format(count/1e6))

    print("Finished!")
    print(uniqcount, "informative residues out of", length, "total positions.")
    return {s: "".join(newdata[s]) for s in newdata}


def select_sites(prefix, seqdata, size, num):
    first = list(seqdata.keys())[0]
    sites = list(range(0, len(seqdata[first])))
    print("Beginning jackknife replicates of size {}...".format(str(size)))
    for i in range(0, num):
        print("jackknife replicate {} of {}...".format(str(i + 1), str(num)))
        jackknife = {s: [] for s in seqdata}
        selected = random.sample(sites, size)
        [sites.remove(s) for s in selected]
        for s in selected:
            for seq in seqdata:
                jackknife[seq].append(seqdata[seq][s])
        o = open(os.path.join(prefix+"_{}.faa".format(str(i + 1))), 'w')
        for j in jackknife:
            o.write(">{}\n{}\n".format(j, "".join(jackknife[j])))
        o.close()
    return


def remove_gapped(seqdata, t):
    first = list(seqdata.keys())[0]
    newdata = {s: [] for s in seqdata.keys()}

    length = len(seqdata[first])
    print(length, "residues to scan...")
    count = 0
    gap_totalcount = 0
    for i in range(0, length):
        gap_rescount = 0
        for s in seqdata:
            if seqdata[s][i] == "-":
                gap_rescount += 1

        prop = gap_rescount/len(seqdata.keys())
        if prop < t:
            for s in seqdata:
                newdata[s].append(seqdata[s][i])
            gap_totalcount += 1

        count += 1
        if count % 100000 == 0:
            if count < 1e6:
                print("{}00K residues parsed...".format(
                    str(int(count/100000)))
                )
            else:
                print("{:.1f}M residues parsed...".format(count/1e6))

    print("Done!")
    print(gap_totalcount, "informative residues out of", length,
          "total positions.")
    return {s: "".join(newdata[s]) for s in newdata}


def main():
    args = parse_args()
    align_file = os.path.abspath(args.align_file)
    prefix = os.path.abspath(args.prefix)

    seqdata = parse(align_file)

    if args.remove_gapped:
        t = args.remove_gapped
        seqdata = remove_gapped(seqdata, t)

    if args.remove_redundant:
        seqdata = remove_redundant(seqdata)

    if args.size:
        size = args.size
    else:
        size = 10000

    if args.num:
        num = args.num
    else:
        num = 1

    select_sites(prefix, seqdata, size, num)


if __name__ == '__main__':
    main()
