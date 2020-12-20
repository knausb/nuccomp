#!/usr/bin/env python3

# https://docs.python.org/3/tutorial/
# https://docs.python.org/3/howto/argparse.html#id1
import argparse
import ntpath
from Bio import SeqIO
#from Bio.Alphabet import IUPAC


parser = argparse.ArgumentParser(description='Determine nucleotide composition of FASTA files. Requires python >=3.7.3 and Biopython >= 1.78.')
parser.add_argument("FASTA", help="FASTA file containing nucleotides.")
parser.add_argument("-v", "--verbose", help="increase output verbosity",
                    action="store_true")
                    
args = parser.parse_args()

##### ##### ##### ##### #####
# Functions

#read_fasta

##### ##### ##### ##### #####
#
# Main.
#
##### ##### ##### ##### #####

if args.verbose:
    print("verbosity turned on")

if args.verbose:
    print("Input FASTA file: " + args.FASTA)

outfile = ntpath.basename(args.FASTA)
outfile = outfile + '_nuccomp.csv'

if args.verbose:
    print("Output file: " + outfile)

with open(outfile, 'w') as f:
    nucls = "aAcCgGtTwWsSmMkKrRyYnN" # Nucleotides to get counts of
    print('Id', 'Length', *nucls, sep = ',', file=f)
    for record in SeqIO.parse(args.FASTA, "fasta"):
        counts = [ record.seq.count(b) for b in nucls ] # Get count of each nucl and store in list
        print(record.id, len(record), *counts, sep = ',', file=f)


# EOF.