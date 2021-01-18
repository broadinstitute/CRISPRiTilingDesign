###############################################################################
## Combining pegRNA pools for oligo orders
## Jesse Engreitz
## January 18, 2021
## Tested on Stanford Sherlock cluster with conda environment

import CombineGuidePools

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
import os.path
from Bio.Seq import Seq

import warnings
warnings.filterwarnings("ignore")

from JuicerCRISPRiDesign import *


def parseargs(required_args=True):
    class formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass

    epilog = ('''
Key columns in output file:
GuideSequenceMinusG     Sequence of the gRNA minus leading G - this should be completed by all guides, including negative controls
CoreOligo               Sequence of gRNA + flanking promoter + scaffold sequences
OligoSequence           Full oligo sequence to order, including subpool primers
subpool                 Name of the subpool of gRNAs with unique primer handles to amplify 
FwdPrimer               Primer 1 for amplifying the subpool (order this sequence)
RevPrimer               Primer 2 for amplifying the subpool (order this sequence)
guideSet                Gene / enhancer / target of this gRNA

Current behavior is, if design file has identical oligo sequences, to collapse these before re-duplicating to fill the space.
''')
    parser = argparse.ArgumentParser(description='''
Combine different oligo subpools output by MakeGuidePool.py into a final pool for ordering.  
Terminology: "subpool" refers to a set of oligos with the same PCR handles on the outside .. i.e. a set of oligos that will be amplified together.
''',
                                     epilog=epilog,
                                     formatter_class=formatter)
    readable = argparse.FileType('r')
    
    parser.add_argument('--config', required=required_args, help="Config file with columns: poolgroup pool DesignFile Multiply FwdPrimer RevPrimer")
    parser.add_argument('--outbase', required=required_args)
    parser.add_argument('--fillToOligoPoolSize', default=-1, type=int, help="Total size of oligo pool; if total number of guides is less than this number, will output oligo order file with this total number of gRNAs.  If total guides is more than this number, will truncate.")
    parser.add_argument('--includeReverseComplements', default=False, action='store_true', help="Include this flag to adding reverse complement oligos. Probably a bad idea if this array includes tiling sequences (or HyPR barcodes)")

    args = parser.parse_args()
    return(args)



def addHandles(df, fwdPrimer, revPrimer):
    ''' Adds PCR handles for PCR1 '''
    if np.isnan(fwdPrimer):
        fwdPrimer = ''
    if np.isnan(revPrimer):
        revPrimer = ''

    df["FwdPrimer"] = fwdPrimer
    df["RevPrimer"] = revPrimer
    df["RevPrimerOligo"]=df["RevPrimer"].apply(lambda x: str(Seq(x).reverse_complement()))
    df["OligoSequence"]=df["FwdPrimer"]+df["CoreOligo"]+df["RevPrimerOligo"]
    return df


def loadSubpool(poolconfig):
    currPool = read_table(poolconfig['DesignFile'])
    currPool = addHandles(currPool, poolconfig['FwdPrimer'], poolconfig['RevPrimer'])
    for col in poolconfig.index:
        if not col in currPool:
            currPool[col] = poolconfig[col]
    return currPool


def writeDesignFile(merged, outfile):
    DESIGNCOLS=["chr", "start", "end","name", "score", "strand", "GuideSequence", "GuideSequenceMinusG",
            "MappingSequence", "OffTargetScore", "target", "subpool", "OligoID"]
    for col in DESIGNCOLS:
        if not col in merged:
            merged[col] = ''
    design = merged[DESIGNCOLS]
    design.to_csv(outfile, sep='\t', header=True, index=False)


def writePcrPrimers(merged, outfile):
    primerpairs = merged[['subpool','FwdPrimer','RevPrimer']].drop_duplicates()
    primers = pd.DataFrame( {
        'PrimerName': "FWD-"+primerpairs['subpool'],
        'Sequence': primerpairs['FwdPrimer']
        })
    primers = primers.append(pd.DataFrame( {
        'PrimerName': "REV-"+primerpairs['subpool'],
        'Sequence': primerpairs['RevPrimer']
        }))
    primers.to_csv(outfile, sep='\t', header=True, index=False)


def writeSubpoolSummary(merged, outfile):
    summ = merged.copy()
    summ['nGuides'] = summ.groupby('subpool')['subpool'].transform('count')
    summ = summ[['subpool','FwdPrimer','RevPrimer','nGuides']].drop_duplicates()
    summ.to_csv(outfile, sep='\t', header=True, index=False)


def writeSequencesToOrder(oligos, outfile, poolMax, includeReverseComplements):
    if poolMax > 0 and len(oligos) > poolMax:
        print("Truncating oligo list from " + str(len(oligos)) + " to " + str(poolMax))
        towrite = oligos[0:poolMax]

    elif poolMax > 0 and len(oligos) < poolMax:
        print("Adding copies of oligos (and reverse complements) to bring pool from " + str(len(oligos)) + " to " + str(poolMax))
        strand="-"
        towrite = oligos.copy()
        oligosRevComp = pd.Series([str(Seq(oligo).reverse_complement()) for key,oligo in oligos.iteritems()])
        while len(towrite) < poolMax:
            nToAdd = min(len(oligos), poolMax-len(towrite))
            if strand == "-" and includeReverseComplements:
                towrite = towrite.append(oligosRevComp[0:nToAdd])
                strand = "+"
            elif strand == "+":
                towrite = towrite.append(oligos[0:nToAdd])
                strand = "-"

    else:
        towrite = oligos

    towrite.to_csv(outfile, header=False, index=False)


def main(args):
    config = read_table(args.config)

    oligos = pd.Series()
    merged = pd.DataFrame()
    for index, row in config.iterrows():
        currPool = loadSubpool(row)
        currOligos = currPool['OligoSequence'].drop_duplicates()

        merged = merged.append(currPool)
        ## Print oligos multiple times if desired, to balance the subpools
        if 'Multiply' in row:
            multiples = row['Multiply']
        else:
            multiples = 1
        for i in range(multiples):
            oligos = oligos.append(currOligos)

    merged.to_csv(args.outbase + ".full.txt", sep='\t', header=True, index=False)
    
    ## Write final files
    writeDesignFile(merged, args.outbase + ".design.txt")
    writePcrPrimers(merged, args.outbase + ".primers.txt") 
    writeSubpoolSummary(merged, args.outbase + '.subpools.txt')
    writeSequencesToOrder(oligos, args.outbase + ".SequencesToOrder.txt", args.fillToOligoPoolSize, args.includeReverseComplements)

    print("Total unique oligos: ", len(oligos.drop_duplicates()))
    print("Total oligos for order: ", len(oligos))

if __name__ == '__main__':
    args = parseargs()
    main(args)
