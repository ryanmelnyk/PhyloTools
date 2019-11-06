#Ryan A. Melnyk
#schmelnyk@gmail.com

import argparse, os
import numpy as np


def parse_args():
	parser = argparse.ArgumentParser(description='''
Given a list of strains of interest, pull out information on unique genes.
	''')
	parser.add_argument('pypdir',type=str,help="path to PyParanoid directory")
	parser.add_argument('group',type=str,help="group ID")
	return parser.parse_args()


def extract_hits(pypdir,group):
	mat_file = open(os.path.join(pypdir,"homolog_matrix.txt"),'r')
	header = mat_file.readline().rstrip().split("\t")
	for line in mat_file:
		if line.startswith(group):
			vals = line.rstrip().split('\t')

	return dict(zip(header[1:],vals[1:]))


def dump_hits(hits):
	for x in sorted(hits.items(), key = lambda x: x[1],reverse=True):
		print("{}\t{}".format(x[0],x[1]))
	return

def main():
	args = parse_args()
	pypdir = os.path.abspath(args.pypdir)
	group = args.group

	hits = extract_hits(pypdir,group)
	dump_hits(hits)

if __name__ == '__main__':
	main()
