#!/usr/bin/env python

import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(description='''
Extract protein sequences from a genbank file.
    ''')
    parser.add_argument('gbk', type=str, help='Path to Genbank file')
    parser.add_argument('faa', type=str, help='Path to Output file')
    parser.add_argument('--locus_tag_field', type=str,
                        help="alternate to locus_tag qualifier")
    return parser.parse_args()


def main():
    args = parse_args()

    if args.locus_tag_field:
        tag_field = args.locus_tag_field
    else:
        tag_field = "locus_tag"

    pep = open(args.faa, 'w')
    for seq in SeqIO.parse(open(args.gbk, 'r'), 'genbank'):
        for feat in seq.features:
            try:
                if feat.type == "CDS":
                    if feat.qualifiers["translation"][0]:
                        aa_seq = SeqRecord(
                            Seq(feat.qualifiers["translation"][0])
                        )
                        aa_seq.id = feat.qualifiers[tag_field][0]
                        aa_seq.description = feat.qualifiers["product"][0]
                        SeqIO.write(aa_seq, pep, 'fasta')
            except KeyError:
                pass

    pep.close()


if __name__ == '__main__':
    main()
