#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import os, argparse, string


def parse_args():
	parser = argparse.ArgumentParser(description='''
Convert stockholm format to fasta.
	''')
	parser.add_argument('infile', type=str,help='alignment file')
	parser.add_argument('outfile', type=str,help='output file to write filtered alignment')
	return parser.parse_args()


def main():
	args = parse_args()
	align_data = {}
	str_tbl=str.maketrans("","",string.ascii_lowercase + ".")
	for line in open(args.infile,'r'):
		if line.startswith("#") or line.startswith("//"):
			continue
		else:
			vals = line.rstrip().split()
		if len(vals) < 1:
			continue
		elif vals[0] in align_data:
			align_data[vals[0]].append(vals[1].translate(str_tbl))
		else:
			align_data[vals[0]] = [vals[1].translate(str_tbl)]
	o = open(args.outfile,'w')
	for a in align_data:
		o.write(">{}\n{}\n".format(a,"".join(align_data[a])))
	o.close()


if __name__ == '__main__':
	main()
