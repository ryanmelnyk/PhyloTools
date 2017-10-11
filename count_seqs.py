from Bio import SeqIO
import sys

lengths = {}
for seq in SeqIO.parse(open(sys.argv[1],'r'),'fasta'):
	l = len(seq.seq)
	if l not in lengths:
		lengths[l] = 1
	else:
		lengths[l] += 1

for l in lengths:
	print lengths[l], "sequences of length", l
