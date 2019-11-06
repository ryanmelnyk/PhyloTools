#!/usr/bin/env python
#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import os,argparse

def parse_args():
	parser = argparse.ArgumentParser(description='''
Using a tab-separated list of human readable headers to replace tree taxon names.
	''')
	parser.add_argument('tree', type=str,help='newick tree file')
	parser.add_argument('names',type=str,help='tab-separated name list')
	parser.add_argument('outf',type=str,help='name of new tree file')
	return parser.parse_args()


def main():
	args = parse_args()

	renamed = {v[0] : v[1] for v in [line.rstrip().split("\t") for line in open(os.path.abspath(args.names),'r').readlines()]}

	tree = open(os.path.abspath(args.tree),'r').read()

	for r in renamed:
		tree = tree.replace(r+":",renamed[r]+":")

	o = open(os.path.abspath(args.outf),'w')
	o.write(tree+"\n")
	o.close()


if __name__ == '__main__':
	main()
