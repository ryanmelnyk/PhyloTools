#!/usr/bin/env python
#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import os, argparse
from Bio import SeqIO

def parse_args():
	parser = argparse.ArgumentParser(description='''
Filter an alignment based on a strainlist
	''')
	parser.add_argument('infile', type=str,help='alignment file')
	parser.add_argument('outfile', type=str,help='output file to write filtered alignment')
	parser.add_argument('strainlist',type=str,help='list of strains')
	return parser.parse_args()

def main():
	args = parse_args()
	infile = open(os.path.abspath(args.infile),'r')
	outfile = open(os.path.abspath(args.outfile),'w')
	strains = [line.rstrip() for line in open(os.path.abspath(args.strainlist),'r')]


	for seq in SeqIO.parse(infile,'fasta'):
		if str(seq.id) in strains:
			outfile.write(">{}\n{}\n".format(seq.id,str(seq.seq)))

if __name__ == '__main__':
	main()
