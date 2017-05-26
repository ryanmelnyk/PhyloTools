
#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import os, argparse
from Bio import SeqIO

def parse_args():
	parser = argparse.ArgumentParser(description='''
Filter an alignment based on highly gapped sequences
	''')
	parser.add_argument('infile', type=str,help='alignment file')
	parser.add_argument('outfile', type=str,help='output file to write filtered alignment')
	parser.add_argument('threshold',type=float,help='minimum proportion of gaps')
	return parser.parse_args()

def main():
	args = parse_args()
	infile = open(os.path.abspath(args.infile),'r')
	outfile = open(os.path.abspath(args.outfile),'w')
	t = args.threshold

	count = 0
	for seq in SeqIO.parse(infile,'fasta'):
		gapprop = float(str(seq.seq).count("-"))/float(len(str(seq.seq)))
		if gapprop > t:
			print "{} has too many gaps. ({})".format(seq.id,str(gapprop))
			count += 1
		else:
			outfile.write(">{}\n{}\n".format(seq.id,str(seq.seq)))

	print count, "total strains removed."

if __name__ == '__main__':
	main()
