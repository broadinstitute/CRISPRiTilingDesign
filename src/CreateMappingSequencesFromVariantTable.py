###############################################################################
## Jesse Engreitz
## September 29, 2021

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
    parser = argparse.ArgumentParser(description='Input: Variant table from CreateVariantTable.py. Output: same table with MappingSequence column added.',
                                     epilog=epilog,
                                     formatter_class=formatter)
    readable = argparse.FileType('r')
    
    parser.add_argument('--input', required=True, help="Input variant table")
    parser.add_argument('--genomeFasta', required=False, help="Genome FASTA file, if using a BED file")
    parser.add_argument('--output', required=True, help="Output variant table with MappingSequence column added")
    parser.add_argument('--bufferLength', type=int, default=5, help="Specify number of nucleotides on each side to output for mapping sequence.")
    args = parser.parse_args()
    return(args)


###############################################

def addMappingSequences(variantTable, bufferLength, fasta):
    ## For each row of the variant table, add a column corresponding to a 'mapping sequence' for downstream VFF analysis,
    ##  which contains the edited sequence plus some number of additional nucleotides on each side

    variantTable['MappingSequence'] = ""
    for index, row in variantTable.iterrows():
        editedSeq = BedTool.seq((row['chr'], row['start']-bufferLength, row['start']), fasta) + \
                     row['alt'] + \
                     BedTool.seq((row['chr'], row['end'], row['end']+bufferLength), fasta) 
        variantTable.loc[index,'MappingSequence'] = editedSeq

    return(variantTable)



##############################################
def main(args):

    ## Read input BED4 file
    variants = read_table(args.input)

    ## Get FASTA file to get sequences
    fasta = pybedtools.example_filename(args.genomeFasta)

    results = addMappingSequences(variants, args.bufferLength, fasta)

    results.to_csv(args.output, sep='\t', header=True, index=False)



##############################################
if __name__ == '__main__':
    args = parseargs()
    main(args)
