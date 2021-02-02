import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

pep = open(sys.argv[2], 'w')
for seq in SeqIO.parse(open(sys.argv[1], 'r'), 'genbank'):
    for feat in seq.features:
        try:
            if feat.type == "CDS":
                if feat.qualifiers["translation"][0]:
                    aa_seq = SeqRecord(Seq(feat.qualifiers["translation"][0]))
                    aa_seq.id = feat.qualifiers["locus_tag"][0]
                    aa_seq.description = feat.qualifiers["product"][0]
                    SeqIO.write(aa_seq, pep, 'fasta')
        except KeyError:
            print("Key missing")

pep.close()
