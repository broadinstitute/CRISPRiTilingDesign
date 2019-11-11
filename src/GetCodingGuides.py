###############################################################################
## Select sets of coding guides from Broad GPP Brunello or Brie
## Jesse Engreitz
## November 10, 2019
## Tested with:  "use .python-3.5.1; source /seq/lincRNA/Ben/VENV_MIP/bin/activate"


from __future__ import division
import pandas as pd
from pandas.io.parsers import read_table
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
import numpy as np
import sys
from scipy import stats
import csv

import warnings
warnings.filterwarnings("ignore")

def parseargs(required_args=True):
    class formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass

    epilog = ("")
    parser = argparse.ArgumentParser(description='Select coding gRNAs from Brie or Brunello',
                                     epilog=epilog,
                                     formatter_class=formatter)
    readable = argparse.FileType('r')
    
    parser.add_argument('--genes', required=required_args, help="One gene symbol per row")
    parser.add_argument('--outfile', required=required_args, help="Output design file")
    parser.add_argument('--guideLibrary', required=required_args, help="Path to gzipped Brie or Brunello library")
    #parser.add_argument('--nGuidesPerGene', default=4, help="Number of guides to choose per element (in the guideSet column)")
    args = parser.parse_args()
    return(args)


def main(args):
    allguides = read_table(args.guideLibrary, compression='gzip')
    genes = read_table(args.genes, header=None, names=['Target Gene Symbol'])


    ####### To do: Blat Brie and Brunello to genome and find genome coordinates to fit with other datasets    
    guides = genes.merge(allguides, on='Target Gene Symbol')
    guides['chr'] = 'NA'
    guides['start'] = 'NA'
    guides['end'] = 'NA'
    guides['score'] = 0
    guides['locus'] = 'NA'
    guides['strand'] = ''
    guides['GuideSequenceWithPAM'] = guides['sgRNA Target Sequence'] + guides['PAM Sequence']
    guides['guideSet'] = guides['Target Gene Symbol']
    guides['SSC'] = 0

    ## TO DO: implement the nGuidesPerGene option
    ## TO DO: Warn if there are genes that we don't find genes for

    guides = guides[['chr','start','end','locus','score','strand','GuideSequenceWithPAM','guideSet','SSC']]
    guides.to_csv(args.outfile, sep='\t', header=True, index=False)



if __name__ == '__main__':
    args = parseargs()
    main(args)
