#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import os, argparse
import ete3
from itertools import combinations

def parse_args():
	parser = argparse.ArgumentParser(description='''
Given a tree, extract a distance matrix and reduce the dataset below a certain threshold
by eliminating very similar strains.
	''')
	parser.add_argument('tree', type=str,help='path to tree file')
	parser.add_argument('mindist', type=float,help='average distance to compress')
	parser.add_argument('output',type=str,help="path to write output strainlist to")
	return parser.parse_args()

def main():
	args = parse_args()
	t = ete3.Tree(os.path.abspath(args.tree))
	mindist = args.mindist
	o = open(os.path.abspath(args.output),'w')
	good_strains = []
	to_skip = []
	for node in t.iter_descendants("preorder"):

		if node in to_skip:
			continue
		else:
			leafnodes = [x for x in node.get_leaves()]
			print(node.name,len(leafnodes))
			if len(leafnodes) > 500:
				continue
			elif len(leafnodes) > 1:
				pairs = [p for p in combinations(leafnodes,2)]
				dist = 0.0
				for p in pairs:
					dist += node.get_distance(p[0],p[1])
				if dist/len(pairs) < mindist:
					good_strains.append(leafnodes[0].name)
					[to_skip.append(desc) for desc in node.iter_descendants("preorder")]
			else:
				if node.is_leaf():
					good_strains.append(node.name)

	print(len(good_strains), "at threshold", mindist)
	o.write("\n".join(good_strains))

if __name__ == '__main__':
	main()
