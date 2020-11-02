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
#print("outfile: " + outfile)
#outfile = args.FASTA.split('.')[0]
#outfile = outfile.split('.')[0]
outfile = outfile + '_nuccomp.csv'

if args.verbose:
    print("Output file: " + outfile)

# os.linesep

with open(outfile, 'w') as f:

    #f.write('Id,Length,a,A,c,C,g,G,t,T,w,W,s,S,m,M,k,K,r,R,y,Y,n,N\n')
    nucls = "aAcCgGtTwWsSmMkKrRyYnN" # Nucleotides to get counts of
    print('Id', 'Length', *nucls, sep = ',', file=f)
    for record in SeqIO.parse(args.FASTA, "fasta"):
        counts = [ record.seq.count(b) for b in nucls ] # Get count of each nucl and store in list
        print(record.id, len(record), *counts, sep = ',', file=f)
        """
        f.write(record.id)
        f.write(',')
        f.write(str(len(record)))
        f.write(',')
        f.write(str(record.seq.count("a")))
        f.write(',')
        f.write(str(record.seq.count("A")))
        f.write(',')
        f.write(str(record.seq.count("c")))
        f.write(',')
        f.write(str(record.seq.count("C")))
        f.write(',')
        f.write(str(record.seq.count("g")))
        f.write(',')
        f.write(str(record.seq.count("G")))
        f.write(',')
        f.write(str(record.seq.count("t")))
        f.write(',')
        f.write(str(record.seq.count("T")))
        f.write(',')
        f.write(str(record.seq.count("w")))
        f.write(',')
        f.write(str(record.seq.count("W")))
        f.write(',')
        f.write(str(record.seq.count("s")))
        f.write(',')
        f.write(str(record.seq.count("S")))
        f.write(',')
        f.write(str(record.seq.count("m")))
        f.write(',')
        f.write(str(record.seq.count("M")))
        f.write(',')
        f.write(str(record.seq.count("k")))
        f.write(',')
        f.write(str(record.seq.count("K")))
        f.write(',')
        f.write(str(record.seq.count("r")))
        f.write(',')
        f.write(str(record.seq.count("R")))
        f.write(',')
        f.write(str(record.seq.count("y")))
        f.write(',')
        f.write(str(record.seq.count("Y")))
        f.write(',')
        f.write(str(record.seq.count("n")))
        f.write(',')
        f.write(str(record.seq.count("N")))
        f.write('\n')
        """
#    count = count + 1
#        print("%s %i" % (record.id, len(record)))
#        print("%s %i %i" % ("aA:", record.seq.count("a"), record.seq.count("A")))
#        print("%s %i %i" % ("cC:", record.seq.count("c"), record.seq.count("C")))
#        print("%s %i %i" % ("gG:", record.seq.count("g"), record.seq.count("G")))
#        print("%s %i %i" % ("tT:", record.seq.count("t"), record.seq.count("T")))
#        print("%s %i %i" % ("wW:", record.seq.count("w"), record.seq.count("W")))
#        print("%s %i %i" % ("sS:", record.seq.count("s"), record.seq.count("S")))
#        print("%s %i %i" % ("mM:", record.seq.count("m"), record.seq.count("M")))
#        print("%s %i %i" % ("kK:", record.seq.count("k"), record.seq.count("K")))
#        print("%s %i %i" % ("rR:", record.seq.count("r"), record.seq.count("R")))
#        print("%s %i %i" % ("yY:", record.seq.count("y"), record.seq.count("Y")))
#        print("%s %i %i" % ("nN:", record.seq.count("n"), record.seq.count("N")))
#    print(record.seq.count("a"), record.seq.count("A"))
#    print(record.seq.count("w"), record.seq.count("W"))
#    print(record.seq.count("n"), record.seq.count("N"))

#    my_seq = record.seq
#    print(my_seq.count("n"), my_seq.count("N"))
#    print(my_seq[0:10])
#    print(my_seq.count("A"))

# EOF.
