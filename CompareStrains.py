#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import argparse, os
import numpy as np


def parse_args():
	parser = argparse.ArgumentParser(description='''
Given a list of strains of interest, pull out information on unique genes.
	''')
	parser.add_argument('strains', type=str,help='path to list of strains')
	parser.add_argument('pypdir',type=str,help="path to PyParanoid directory")
	parser.add_argument('prefix',type=str,help="prefix for output files")
	return parser.parse_args()

def get_corr_coef(strains,pypdir, prefix):
	infile = open(os.path.join(pypdir,"homolog_matrix.txt"),'r')
	header_line = infile.readline().rstrip().split("\t")
	indices = [header_line.index(s) for s in strains]
	# print dict(zip(indices,strains))

	lines = []
	groups = []
	for line in infile:
		vals = line.rstrip().split("\t")
		# convert input data into boolean presence/absence
		lines.append([int(bool(int(vals[i]))) for i in indices])
		groups.append(vals[0])
	a = np.stack(lines)
	cc = np.corrcoef(a.T)
	o = open(prefix+".corrcoef.txt",'w')
	o.write("\t{}\n".format("\t".join(strains)))
	for i in range(0,len(strains)):
		o.write("{}\t{}\n".format(strains[i],"\t".join([str(x) for x in cc[i]])))
	return a, groups

def find_unique(a, groups, strains):
	unique = {}
	missing = {}
	common = []
	for s in strains:
		unique[s] = []
		missing[s] = []
	for i in range(0,a.shape[0]):
		nz = np.nonzero(a[i])
		if len(nz[0]) == 1:
			unique[strains[nz[0][0]]].append(groups[i])
		elif len(nz[0]) == len(strains)-1:
			missing[strains[np.nonzero(a[i] < 1)[0][0]]].append(groups[i])
		elif len(nz[0]) == len(strains):
			common.append(groups[i])
		else:
			pass
	print strains
	for s in strains:
		print s
		print "\t", len(unique[s]), "unique"
		print "\t", len(missing[s]), "missing"
	print len(common), "common to all strains."
	return unique,missing, common

def get_gene_info(unique,missing,pypdir,strains,prefix):
	# get a locus tag and a description
	all_groups = {}
	for s in strains:
		for g in unique[s] + missing[s]:
			all_groups[g] = {}

	infile = open(os.path.join(pypdir,"locustag_matrix.txt"),'r')
	header_line = infile.readline().rstrip().split("\t")
	indices = [header_line.index(s) for s in strains]
	for line in infile:
		vals = line.rstrip().split("\t")
		if vals[0] in all_groups:
			all_groups[vals[0]]["locustags"] = [vals[i] for i in indices]
	infile.close()

	for line in open(os.path.join(pypdir,"group_descriptions.txt"),'r'):
		vals = line.rstrip().split("\t")
		if vals[0] in all_groups:
			counts = {s : vals[1:].count(s) for s in list(set(vals[1:]))}
			all_groups[vals[0]]["description"] = sorted(counts.items(), key = lambda x: x[1],reverse=True)[0][0]

	# for g in all_groups:
	# 	if "description" not in all_groups[g]:
	# 		all_groups[g]["description"] = "no description"

	o = open(prefix+".geneinfo.txt",'w')
	o.write("groupname\tdescription\t{}\n".format("\t".join(strains)))
	for s in strains:
		for g in unique[s]:
			o.write("{}\t{}\t{}\n".format(g,all_groups[g]["description"],"\t".join(all_groups[g]["locustags"])))
		for g in missing[s]:
			o.write("{}\t{}\t{}\n".format(g,all_groups[g]["description"],"\t".join(all_groups[g]["locustags"])))
	o.close()

	return



def main():
	args = parse_args()
	strains = [line.rstrip() for line in open(os.path.abspath(args.strains),'r').readlines()]
	pypdir = os.path.abspath(args.pypdir)
	prefix = args.prefix

	a, groups = get_corr_coef(strains, pypdir, prefix)
	unique, missing, common = find_unique(a, groups, strains)
	get_gene_info(unique,missing,pypdir,strains,prefix)

if __name__ == '__main__':
	main()
