###############################################################################
## Script to set up a variant table for input to DesignPrimeEditor.py
## Jesse Engreitz
## October 3, 2020

import pybedtools
from pybedtools import BedTool
import pandas as pd
import numpy as np
import argparse
from pandas.io.parsers import read_table
from Bio.SeqFeature import FeatureLocation
from collections import OrderedDict

def parseargs():
    class formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass

    epilog = ("")
    parser = argparse.ArgumentParser(description='Input: BED file. Output: table with chr start end name ref alt describing Output list of potential prime editor gRNAs for desired variants',
                                     epilog=epilog,
                                     formatter_class=formatter)
    readable = argparse.FileType('r')
    
    parser.add_argument('--input', required=True, help="BED file")
    parser.add_argument('--genomeFasta', required=False, help="Genome FASTA file, if using a BED file")
    parser.add_argument('--output', required=True, help="Variants file")
    parser.add_argument('--window', type=int, help="Window size to mutate")
    parser.add_argument('--offset', type=int, help="Sliding resolution")
    parser.add_argument('--mutagenesis', default=False, help="Include mutagenesis of each subsequence")
    parser.add_argument('--scramble', default=False, help="Include scramble (reverse) of each subsequence")
    parser.add_argument('--includeRef', default=False, help="Include a variant to convert ref to ref nucleotide (no change)")

    args = parser.parse_args()
    return(args)


###############################################
def addVariant(edits, chr, currStart, currEnd, currSeq, newSeq):
    edit = pd.Series( OrderedDict(( 
            ('chr', chr), 
            ('start', str(currStart)), 
            ('end', str(currEnd)), 
            ('name', chr + ':' + str(currStart) + ':' + str(currSeq) + '>' + newSeq),
            ('ref', currSeq), 
            ('alt', newSeq))))
    edits = edits.append(edit, ignore_index=True)[edit.index.tolist()]
    return edits

def createVariants(chr, start, end, seq, window, offset, mutagenesis=True, scramble=False, includeRef=True):

    edits = pd.DataFrame()

    currStart = start
    currEnd = start + window
    while currEnd <= end:
        currSeq = (FeatureLocation(currStart, currEnd)+(-start)).extract(seq).upper()

        if mutagenesis:
            for newBase in ['A','T','C','G']:
                newSeq = newBase * (currEnd-currStart)
                if (newSeq != currSeq):
                    edits = addVariant(edits, chr, currStart, currEnd, currSeq, newSeq)

        if scramble:
            ## Scramble by reversing (not reverse complementing) the sequence
            newSeq = currSeq[::-1]
            edits = addVariant(edits, chr, currStart, currEnd, currSeq, newSeq)

        if includeRef:
            newSeq = currSeq
            edits = addVariant(edits, chr, currStart, currEnd, currSeq, newSeq)


        ## Increment the window
        currStart = currStart + offset
        currEnd = currStart + window

    return edits


def createVariantsMutagenesis(chr, start, end, seq, mutagenesis=True, scramble=False, includeRef=True):
    return createVariants(chr, start, end, seq, 1, 1, mutagenesis, scramble, includeRef)


##############################################
def main(args):

    ## Read input BED4 file
    bed = read_table(args.input, names=['chr','start','end','name'])

    ## Get FASTA file to get sequences
    fasta = pybedtools.example_filename(args.genomeFasta)

    results = pd.DataFrame()
    for index, row in bed.iterrows():
        seq = BedTool.seq((row['chr'], row['start'], row['end']), fasta)
        curr = createVariants(
            row['chr'], 
            row['start'], 
            row['end'], 
            seq, 
            args.window, 
            args.offset, 
            args.mutagenesis,
            args.scramble,
            args.includeRef)
        if (len(curr) > 0):
            results = results.append(curr)[curr.columns.tolist()]
    
    results.to_csv(args.output, sep='\t', header=True, index=False)


##############################################
if __name__ == '__main__':
    args = parseargs()
    main(args)
