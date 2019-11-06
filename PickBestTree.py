#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import sys
import os, argparse

def parse_args():
	parser = argparse.ArgumentParser(description='''
From a set of independently calculated RAxML "Best Trees", select the one with
the highest max likelihood value.
	''')
	parser.add_argument('prefix', type=str,help='name of files i.e. RAxML_info.prefix_0')
	parser.add_argument('inferences',type=int,help='number of independent consecutively numbered inferences')
	return parser.parse_args()


def main():
	args = parse_args()
	prefix = args.prefix
	inferences = args.inferences
	scores = []
	for i in range(0,inferences):
		try:
			for line in open("RAxML_info.{}_{}".format(prefix,i),'r'):
				if line.startswith("Final GAMMA-based Score of best tree"):
					scores.append(float(line.rstrip().split()[-1]))
					print "Tree {}: {}".format(i,scores[-1])
		except IOError:
			print "Tree {}: File not found."

	print scores
	print max(scores)
	print scores.index(max(scores))


if __name__ == '__main__':
	main()
