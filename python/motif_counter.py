#!/usr/bin/env python3

# https://docs.python.org/3/tutorial/
# https://docs.python.org/3/howto/argparse.html#id1
import argparse
from Bio import SeqIO
import gzip
import ntpath
import os
import re

parser = argparse.ArgumentParser(
    prog='motif_counter.py',
    description='Count motifs in windows of FAST[AQ] files. Requires python >=3.7.3 and Biopython >= 1.78.')

parser.add_argument("INFILE", help="FAST[AQ] file containing nucleotides.")
parser.add_argument('--motif', nargs='?', default="CG", const="CG", type=str, help="Motif to count in each window [CG].")
parser.add_argument('--invert', dest='invert', action='store_true', help="Invert the scaled counts [default].")
parser.add_argument('--no-invert', dest='invert', action='store_false', help="Do not invert the scaled counts.")
parser.set_defaults(invert=True)
parser.add_argument('--win_size', nargs='?', const=1000000, type=str, default=1000000, help="Window size [default = 1000000]")
parser.add_argument("-v", "--verbose", help="increase output verbosity.",
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


def chrom_to_win(chrom, win_size, motif, out_file):
    chrom.seq = chrom.seq.upper()
    # Initialize the number of windows in this chromosome.
    # Taking the int() is equivalent to floor(), so we add one.
    win_number = int( len(chrom) / win_size )
    win_number = win_number + 1

    # Initialize lists for chromosomal summaries.
    my_starts = [0] * win_number
    my_ends = [0] * win_number
    my_lens = [0] * win_number
    my_counts = [0] * win_number
    my_scores = [0] * win_number
    my_rgbs = [0] * win_number

    # Initialize a min and max for scaling.
    chrom_min = -999
    chrom_max = -999

    # R> dput(apply(col2rgb(viridisLite::magma( n=11, alpha = 1, begin = 0.2, end = 0.9, direction  =1)), MARGIN = 2, function(x){ paste(x, collapse = ",") }))
    viridis_magma = ["59,15,112", "89,21,126", "116,33,129",
                     "145,43,129", "173,52,124", "203,62,114",
                     "228,79,100", "245,107,92", "252,140,99",
                      "254,173,119", "254,206,145"]

    viridis_magma.reverse()

    # Loop over windows
    # to count regex matches
    # and collect chromosome max and min.
    for i in range( win_number ):
        # Window coordinates.
        my_starts[i] = i * win_size
        my_ends[i] = my_starts[i] + win_size - 1
        if my_ends[i] > len(chrom):
            my_ends[i] = len(chrom)
        my_seq = chrom.seq[my_starts[i]:my_ends[i] + 1]
        my_lens[i] = len(my_seq)
        my_counts[i] = len( re.findall( motif, str(my_seq) ) )
        my_scores[i] = my_counts[i] / my_lens[i]
        # Check chromosome minimum and maximum values, store for scaling.
        if i == 0:
            chrom_min = my_scores[i]
            chrom_max = my_scores[i]
        if my_scores[i] < chrom_min:
            chrom_min = my_scores[i]
        if my_scores[i] > chrom_max:
            chrom_max = my_scores[i]

    for i in range( win_number ):
            if (chrom_max - chrom_min) > 0:
                my_scores[i] = my_scores[i] - chrom_min
                my_scores[i] = my_scores[i] / ( chrom_max - chrom_min )
            if args.invert:
                my_scores[i] = 1 - my_scores[i]
            my_scores[i] = my_scores[i] * 1000
            my_scores[i] = int( my_scores[i] )
            my_index = int( my_scores[i] / 100 )
            my_rgbs[i] = viridis_magma[ my_index ]


    for i in range( win_number ):
        print( chrom.id,
               my_starts[i],
               my_ends[i],
               ".",
               my_scores[i],
               ".",
               ".",".",
               my_rgbs[i],
               my_counts[i],
               my_lens[i],
               sep = "\t", file = out_file)



##### ##### ##### ##### #####
#
# Main.
#
##### ##### ##### ##### #####

if args.verbose:
    print("verbosity turned on")

win_size = args.win_size
win_size = float(win_size)
win_size = int(win_size)

file_dict = check_file(args.INFILE)

file_dict[ 'outname' ] = file_dict[ 'outname' ] + "_" + str(args.motif) + '_wins.bed'

if args.verbose:
    print("Input file: "  + file_dict[ 'inname' ])
    print("Output file: " + file_dict[ 'outname' ])
    print("File compression: " + file_dict[ 'compression' ])
    print("Window size: " + str(win_size))
    print("Invert score: " + str(args.invert))


if file_dict[ 'compression' ] == 'GZ':
    f = gzip.open( file_dict[ 'inname' ], 'rt' )
else:
    f = open( file_dict[ 'inname' ], encoding="utf-8" )

fout = open(file_dict[ 'outname' ], 'w')
print("# ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRGB', 'blockCount', 'blockSizes']",
      sep = '\t', file=fout)
print("# c('chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRGB', 'blockCount', 'blockSizes')",
      sep = '\t', file=fout)
print("# motif_counter.py --motif ", args.motif, " --win_size ", args.win_size, " ", args.INFILE,
      sep = '', file=fout)

if( file_dict[ 'format' ] == 'FASTA' ):
    for record in SeqIO.parse( f, "fasta"):
        chrom_to_win(record, win_size, args.motif, out_file = fout )
if( file_dict[ 'format' ] == 'FASTQ' ):
    for record in SeqIO.parse( f, "fastq"):
        chrom_to_win(record, win_size, args.motif, out_file = fout )


f.close()
fout.close()



# EOF.
