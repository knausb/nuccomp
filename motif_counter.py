#!/usr/bin/env python3

# https://docs.python.org/3/tutorial/
# https://docs.python.org/3/howto/argparse.html#id1
import argparse
from Bio import SeqIO
import gzip
import ntpath
import os
import re

parser = argparse.ArgumentParser(description='Count motifs in windows of FAST[AQ] files. Requires python >=3.7.3 and Biopython >= 1.78.')
parser.add_argument("INFILE", help="FAST[AQ] file containing nucleotides.")
parser.add_argument('--motif', nargs='?', default="CG", const="CG", type=str, help="Motif to count in each window (CG).")
parser.add_argument('--invert', action='store_true', help="Invert the scaled counts [True]")
parser.add_argument('--no-invert', dest='invert', action='store_false')
parser.set_defaults(invert=True)
#parser.add_argument('--win_size', nargs='?', const=1000000, type=int, default=1000000)
parser.add_argument('--win_size', nargs='?', const=1000000, type=str, default=1000000)
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
#    print("file_properties[ 'outname' ]", file_properties[ 'outname' ])

    if is_gz_file(infile):
#        print("File is gzipped.")
        file_properties[ 'compression' ] = 'GZ'
        file_properties[ 'outname' ] = re.sub(".gz$", "", file_properties[ 'outname' ])
        f = gzip.open(infile, 'rt')
    else:
#        print("File is not gzipped.")
        f = open(infile, encoding="utf-8")
    line = f.readline()
    f.close()
    line = line.rstrip()
    if re.match("^@", line):
#        print("FASTQ line:", line)
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
#        print("FASTA line:", line)
    else:
        print("Unexpected line:", line)
#        raise ValueError('Expecting a FASTQ format file.')
        sys.exit("Expecting a FASTA or FASTQ format file.")


def chrom_to_win(chrom, win_size, motif, out_file):
    chrom.seq = chrom.seq.upper()
    #chrom.seq = chrom.seq.lower()
    win_number = int( len(chrom) / win_size )
    win_number = win_number + 1
    # print( "win_number: ", win_number )
    my_starts = [0] * win_number
    my_ends = [0] * win_number
    my_lens = [0] * win_number
    my_counts = [0] * win_number
    my_scores = [0] * win_number
    my_rgbs = [0] * win_number
    chrom_min = 1000000
    chrom_max = 0
    
#    viridis_magma = ["59,15,112", "96,24,128", "131,38,129",
#                     "168,50,125", "205,64,113", "235,87,96",
#                     "250,127,94", "254,170,116", "254,212,150",
#                     "252,253,191"]
                    
    viridis_magma = ['59,15,112', '92,22,127', '124,35,130',
                     '156,46,127','190,58,119', '222,73,104',
                     '244,102,92', '252,140,99', '254,178,122',
                     '254,216,154','252,253,191']

#    print("win_number: ", win_number)
    # Loop over windows.
    for i in range( win_number ):
        # Window coordinates.
        my_starts[i] = i * win_size
        my_ends[i] = my_starts[i] + win_size - 1
        if my_ends[i] > len(chrom):
            my_ends[i] = len(chrom)
        my_seq = chrom.seq[my_starts[i]:my_ends[i] + 1]
        my_lens[i] = len(my_seq)
        my_counts[i] = len( re.findall( motif, str(my_seq) ) )
        # Check chromosome minimum and maximum values, store for scaling.
        if my_counts[i] < chrom_min:
            chrom_min = my_counts[i]
        if my_counts[i] > chrom_max:
            chrom_max = my_counts[i]


#    print( "chrom_min: ", chrom_min)
#    print( "chrom_max: ", chrom_max)
    # Compute score and RGB columns.
    if win_number == 1:
        i = 0
        my_scores[i] = my_counts[i] / my_lens[i]    
        if args.invert:
            my_scores[i] = 1 - my_scores[i]
        my_scores[i] = my_scores[i] * 1000
        my_scores[i] = int( my_scores[i] )
        my_index = int( my_scores[i] / 100 )
#        print( my_index )
        my_rgbs[i] = viridis_magma[ my_index ]        
    else:
        for i in range( win_number ):
#            print("    my_counts[i]: ", my_counts[i])
            my_scores[i] = my_counts[i] / my_lens[i]
            my_scores[i] = (my_scores[i] - chrom_min/my_lens[i] )
            my_scores[i] = my_scores[i] / ( (chrom_max - chrom_min)/my_lens[i] )
            if args.invert:
                my_scores[i] = 1 - my_scores[i]
            my_scores[i] = my_scores[i] * 1000
            my_scores[i] = int( my_scores[i] )
            my_index = int( my_scores[i] / 100 )
#        print( my_index )
            my_rgbs[i] = viridis_magma[ my_index ]

    
    for i in range( win_number ):
        print( chrom.id,
               my_starts[i],
               my_ends[i],
               ".",
               my_scores[i],
               ".",
               ".",".",
               #".",
               my_rgbs[i],
               my_counts[i],
               my_lens[i], 
#               chrom_min,
#               chrom_max,
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
#    print("Input file: " + args.FASTA)
#    print("Output file: " + outfile)
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
