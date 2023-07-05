#!/usr/bin/env python3

# https://docs.python.org/3/tutorial/
# https://docs.python.org/3/howto/argparse.html#id1
# https://hifiasm.readthedocs.io/en/latest/faq.html#how-do-i-get-contigs-in-fasta

import argparse
import ntpath
import os
import re
import gzip

parser = argparse.ArgumentParser(description='Convert *.gfa file to FASTA.')
#parser.add_argument("FASTA", help="FASTA file containing nucleotides.")
parser.add_argument("INFILE", help="gfa file containing nucleotides.")
parser.add_argument("-v", "--verbose", help="increase output verbosity",
                    action="store_true")

args = parser.parse_args()


##### ##### ##### ##### #####

if args.verbose:
    print("verbosity turned on")

outfile = ntpath.basename(args.INFILE)
outfile = re.sub(".gfa$", "", outfile)
#outfile = outfile + '.fasta'
outfile = outfile + '.fasta.gz'

if args.verbose:
    print("outfile:", outfile)

f = open(args.INFILE, 'r', encoding="utf-8")
#f2 = open(outfile, 'w', encoding="utf-8")
f2 = gzip.open(outfile, 'rt')


for line in f:
    if line.startswith("S"):
#    print(line, end='')
        tmp = line.split()
        print(">", tmp[1], sep = '', file=f2)
        print(tmp[2], sep = '', file=f2)


f.close()
f2.close()


# EOF.
