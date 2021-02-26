#!/usr/bin/env python

from Bio import SeqIO
import sys

if len(sys.argv) > 2:
    seq_format = sys.argv[2]
else:
    seq_format = 'fasta'

lengths = {}
for seq in SeqIO.parse(open(sys.argv[1], 'r'), seq_format):
    seq_len = len(seq.seq)
    if seq_len not in lengths:
        lengths[seq_len] = 1
    else:
        lengths[seq_len] += 1

for seq_len in sorted(lengths.keys()):
    print(lengths[seq_len], "sequences of length", seq_len)
