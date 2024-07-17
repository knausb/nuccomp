#!/usr/bin/env python3

# https://docs.python.org/3/tutorial/
# https://docs.python.org/3/howto/argparse.html#id1
import argparse
from Bio import SeqIO
import gzip
import ntpath
import os
import re


parser = argparse.ArgumentParser(description='Determine nucleotide composition of FASTA files. Requires python >=3.7.3 and Biopython >= 1.78.')
parser.add_argument("INFILE", help="FASTA file containing nucleotides.")
parser.add_argument("-v", "--verbose", help="increase output verbosity",
                    action="store_true")

args = parser.parse_args()

##### ##### ##### ##### #####
#
# Function definitions.
#
##### ##### ##### ##### #####


def is_gz_file(filepath):
    # https://stackoverflow.com/a/47080739
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'


def check_file(infile):
    file_properties = {'inname': infile, 'format': 'none', 'compression': 'none', 'outname': 'none'}
    # Remove PATH from infile.
    file_properties[ 'outname' ] = os.path.basename(infile)

    if is_gz_file(infile):
        file_properties[ 'compression' ] = 'GZ'
        file_properties[ 'outname' ] = re.sub(".gz$", "", file_properties[ 'outname' ])
        f = gzip.open(infile, 'rt')
    else:
        f = open(infile, encoding="utf-8")
    line = f.readline()
    f.close()
    line = line.rstrip()
    if re.match("^@", line):
        file_properties[ 'format' ] = 'FASTQ'
        file_properties[ 'outname' ] = re.sub(".fq$", "", file_properties[ 'outname' ])
        file_properties[ 'outname' ] = re.sub(".fastq$", "", file_properties[ 'outname' ])
        return( file_properties )
    elif re.match("^>", line):
        file_properties[ 'format' ] = 'FASTA'
        file_properties[ 'outname' ] = re.sub(".fa$", "", file_properties[ 'outname' ])
        file_properties[ 'outname' ] = re.sub(".faa$", "", file_properties[ 'outname' ])
        file_properties[ 'outname' ] = re.sub(".fna$", "", file_properties[ 'outname' ])
        file_properties[ 'outname' ] = re.sub(".fsa$", "", file_properties[ 'outname' ])
        file_properties[ 'outname' ] = re.sub(".fasta$", "", file_properties[ 'outname' ])
        return( file_properties )
    else:
        print("Unexpected line:", line)
        sys.exit("Expecting a FASTA or FASTQ format file.")


def count_nucleotides(fh_in, file_dict):
    nucls = "aAcCgGtTwWsSmMkKrRyYnN" # Nucleotides to get counts of
    fout = open(file_dict[ 'outname' ], 'w')
    print('Id', 'Length', *nucls, sep = ',', file=fout)
    if( file_dict[ 'format' ] == 'FASTA' ):
        for record in SeqIO.parse( file_dict[ 'inname' ], "fasta"):
            counts = [ record.seq.count(b) for b in nucls ] # Get count of each nucl and store in list
            print(record.id, len(record), *counts, sep = ',', file=fout)
    if( file_dict[ 'format' ] == 'FASTQ' ):
        for record in SeqIO.parse( file_dict[ 'inname' ], "fastq"):
            counts = [ record.seq.count(b) for b in nucls ] # Get count of each nucl and store in list
            print(record.id, len(record), *counts, sep = ',', file=fout)


##### ##### ##### ##### #####
#
# Main.
#
##### ##### ##### ##### #####

if args.verbose:
    print("verbosity turned on")

file_dict = check_file(args.INFILE)
file_dict[ 'outname' ] = file_dict[ 'outname' ] + '_nuccomp.csv'


if args.verbose:
    print("Input file: "  + file_dict[ 'inname' ])
    print("Output file: " + file_dict[ 'outname' ])
    print("File compression: " + file_dict[ 'compression' ])


if file_dict[ 'compression' ] == 'GZ':
    f = gzip.open( file_dict[ 'inname' ], 'rt' )
else:
    f = open( file_dict[ 'inname' ], encoding="utf-8" )
nucls = "aAcCgGtTwWsSmMkKrRyYnN" # Nucleotides to get counts of
fout = open(file_dict[ 'outname' ], 'w')
print('Id', 'Length', *nucls, sep = ',', file=fout)
if( file_dict[ 'format' ] == 'FASTA' ):
    for record in SeqIO.parse( f, "fasta"):
        counts = [ record.seq.count(b) for b in nucls ] # Get count of each nucl and store in list
        print(record.id, len(record), *counts, sep = ',', file=fout)
if( file_dict[ 'format' ] == 'FASTQ' ):
    for record in SeqIO.parse( f, "fastq"):
        counts = [ record.seq.count(b) for b in nucls ] # Get count of each nucl and store in list
        print(record.id, len(record), *counts, sep = ',', file=fout)


f.close()
fout.close()


# EOF.
