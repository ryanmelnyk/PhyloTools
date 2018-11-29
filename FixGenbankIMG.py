#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import os, argparse,re
from Bio import SeqIO

def parse_args():
	parser = argparse.ArgumentParser(description='''
Replacing the locus tags in the IMG genbank file with accession numbers.
	''')
	parser.add_argument('strain', type=str,help='strain ID for genomedb files')
	parser.add_argument("genomedb",type=str,help='path to genomedb')
	parser.add_argument("output",type=str,help="")
	return parser.parse_args()


def get_tags(strain,genomedb):
	tags = {}
	for seq in SeqIO.parse(open(os.path.join(genomedb,"pep","{}.pep.fa".format(strain)),'r'),'fasta'):
		tags[seq.description.split()[1]] = seq.id
	return tags

def replace_names(strain,genomedb,tags,output):
	o = open(output,'w')
	for line in open(os.path.join(genomedb,"gbk","{}.gbk".format(strain)),'r'):
		vals = line.split()
		if vals[0].startswith("/locus_tag=\""):
			try:
				o.write(line.replace(vals[0].split("\"")[1],tags[vals[0].split("\"")[1]]))
			except KeyError:
				o.write(line)
		else:
			o.write(line)
	o.close()
	return

def main():
	args = parse_args()
	strain = args.strain
	genomedb = os.path.abspath(args.genomedb)
	output = os.path.abspath(args.output)

	tags = get_tags(strain,genomedb)

	replace_names(strain,genomedb,tags,output)


if __name__ == '__main__':
	main()
