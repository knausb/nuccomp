#!/usr/bin/env python3

# https://docs.python.org/3/tutorial/
# https://docs.python.org/3/howto/argparse.html#id1
import argparse
import ntpath
from Bio import SeqIO


parser = argparse.ArgumentParser(description='Counts the occurrence of a motif in FASTA files. Requires python >=3.7.3 and Biopython >= 1.78.')
parser.add_argument("FASTA", help="FASTA file containing nucleotides.")
parser.add_argument("-m", "--motif", default="CTTAAG",
                    help="the motif to count, default = CTTAAG")
parser.add_argument("-v", "--verbose", help="increase output verbosity",
                    action="store_true")
                    
args = parser.parse_args()

##### ##### ##### ##### #####
#
# Main.
#
##### ##### ##### ##### #####

if args.verbose:
    print("verbosity turned on")

if args.verbose:
    print("Input FASTA file: " + args.FASTA)
    print("Motif: " + args.motif)

outfile = ntpath.basename(args.FASTA)
outfile = outfile + '_' + args.motif + '_motifcounts.csv'

if args.verbose:
    print("Output file: " + outfile)


with open(outfile, 'w') as f:
    print('Id', 'Length', "Count", sep = ',', file=f)
    for record in SeqIO.parse(args.FASTA, "fasta"):
        counts = record.seq.count(args.motif) # Get count of each nucl and store in list
        print(record.id, len(record), counts, sep = ',', file=f)

f.close()

# EOF