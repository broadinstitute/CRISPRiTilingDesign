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
from Bio import SeqUtils
from itertools import product
from random import sample


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
    parser.add_argument('--randomBalanced', type=int, default=0, help="Include up to N variants that are random sequences (>=50%% GC content, with >=75%% bases different from original. Will select from list of up to N*3 shared sequences across all mutations, so the same replacements are frequently used")
    parser.add_argument('--reverse', default=False, help="Include variant that is reverse of each subsequence")
    parser.add_argument('--replaceWithN', default=False, help="For each window, include a pegRNA encode an N in the pegRNA rather than mutagenesis or reverse")
    parser.add_argument('--includeRef', default=False, help="Include a variant to convert ref to ref nucleotide (no change)")
    parser.add_argument('--substitute', type=str, default=None, help="Substitute the sequence for (each of) the following replacement sequence(s), comma-delimited")
    args = parser.parse_args()
    return(args)


###############################################
def addVariant(edits, chr, currStart, currEnd, currSeq, newSeq, seqName):
    edit = pd.Series( OrderedDict(( 
            ('chr', chr), 
            ('start', str(currStart)), 
            ('end', str(currEnd)), 
            ('name', chr + ':' + str(currStart) + ':' + str(currSeq) + '>' + newSeq),
            ('ref', currSeq), 
            ('alt', newSeq),
            ('region', seqName))))
    edits = edits.append(edit, ignore_index=True)[edit.index.tolist()]
    return edits


def getRandomBalancedSeq(seqLength, n=1, minGcPct=50):
    ## Dumb algorithm:  Create all possible sequences, then filter (will not scale well for large regions)
    if seqLength >= 10: 
        raise ValueError("randomBalanced algorithm will not scale well to scrambling sequences longer than 10 nucleotides")

    allSeqs = ["".join(x) for x in product('ACGT', repeat=seqLength)]
    allSeqTable = pd.DataFrame({
        'seq' : allSeqs,
        'gcPct' : [SeqUtils.GC(x) for x in allSeqs]
        })

    filteredTable = allSeqTable[allSeqTable['gcPct'] >= minGcPct]

    ## Select n evenly spaced across the table:
    idx = np.round(np.linspace(0, len(filteredTable['seq']) - 1, n)).astype(int)
    selected = [list(filteredTable['seq'])[i] for i in idx]

    return selected


def selectBalancedSeqs(origSeq, replacements, n=1, minNewBasesPct=75):

    ## Find replacement sequences different from original:
    replacementsDf = pd.DataFrame({
        'seq' : replacements,
        'diff' : [sum(c1!=c2 for c1,c2 in zip(s,origSeq)) / len(origSeq) * 100 for s in replacements]
        })

    replacementsDf = replacementsDf[replacementsDf['diff'] >= minNewBasesPct]
    
    ## Select random N from among the replacements:
    result = sample(list(replacementsDf['seq']), min(n,len(replacementsDf)))
    return result 


def createVariants(chr, start, end, seqName, seq, window, offset, args): 

    edits = pd.DataFrame()  ## To do:  Change to edits = [] and use pd.concat() instead of df.append() for better speed

    if (args.randomBalanced > 0):
        replacementSeqs = getRandomBalancedSeq(window, args.randomBalanced*4, minGcPct=50)

    if (args.substitute is not None):
        substitutionSeqs = args.substitute.split(',')

    currStart = start
    currEnd = start + window
    while currEnd <= end:
        currSeq = (FeatureLocation(currStart, currEnd)+(-start)).extract(seq).upper()

        if args.mutagenesis:
            for newBase in ['A','T','C','G']:
                newSeq = newBase * (currEnd-currStart)
                if (newSeq != currSeq):
                    edits = addVariant(edits, chr, currStart, currEnd, currSeq, newSeq, seqName)

        if args.randomBalanced > 0:
            ## Include up to N variants that are random sequences (>=50%% GC content)
            selectedSequences = selectBalancedSeqs(currSeq, replacementSeqs, args.randomBalanced, minNewBasesPct=75)
            if (len(selectedSequences) == 0):
                raise ValueError("Did not find any substitution sequences for random-balanced replacements")
            for newSeq in selectedSequences:
                edits = addVariant(edits, chr, currStart, currEnd, currSeq, newSeq, seqName)

        if args.reverse:
            ## "Scramble" by reversing (not reverse complementing) the sequence
            newSeq = currSeq[::-1]
            edits = addVariant(edits, chr, currStart, currEnd, currSeq, newSeq, seqName)

        if args.replaceWithN:
            ## Replace with Ns
            newSeq = "N" * len(currSeq)
            edits = addVariant(edits, chr, currStart, currEnd, currSeq, newSeq, seqName)

        if args.substitute is not None:
            for newSeq in substitutionSeqs:
                edits = addVariant(edits, chr, currStart, currEnd, currSeq, newSeq, seqName)

        if args.includeRef:
            newSeq = currSeq
            edits = addVariant(edits, chr, currStart, currEnd, currSeq, newSeq, seqName)


        ## Increment the window
        currStart = currStart + offset
        currEnd = currStart + window

    return edits


#def createVariantsMutagenesis(chr, start, end, seqName, seq, mutagenesis=True, reverse=False, replaceWithN=False, includeRef=True):
#    return createVariants(chr, start, end, seq, 1, 1, mutagenesis, reverse, replaceWithN, includeRef)


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
            row['name'], 
            seq,
            args.window,
            args.offset,
            args)
        if (len(curr) > 0):
            results = results.append(curr)[curr.columns.tolist()]
    
    results.to_csv(args.output, sep='\t', header=True, index=False)


##############################################
if __name__ == '__main__':
    args = parseargs()
    main(args)
