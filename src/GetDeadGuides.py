###############################################################################
## Create sets of dead guides (15bp) based on a previously designed pool
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


from JuicerCRISPRiDesign import *

def parseargs(required_args=True):
    class formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass

    epilog = ("")
    parser = argparse.ArgumentParser(description='Create sets of dead gRNAs',
                                     epilog=epilog,
                                     formatter_class=formatter)
    readable = argparse.FileType('r')
    
    parser.add_argument('--design', required=required_args, help="Input design file")
    parser.add_argument('--outfile', required=required_args, help="Output design file")
    parser.add_argument('--PoolID', required=required_args, help="Name of subpool")
    args = parser.parse_args()
    return(args)



def convertDesignToDeadGuide(design, PoolID, SEQCOL='seq', MapLength=21, finalLength=18):
    lengths = design.GuideSequenceWithPAM.str.len()
    if not all(lengths == 23):
        raise ValueError("Currently only supports 20bp guides (23bp with PAM)")

    design['strand'] = design['strand'].fillna('')
    design.loc[design['strand'] == "+",'start'] = design.loc[design['strand'] == "+",'start'] + 5
    design.loc[design['strand'] == "-",'end'] = design.loc[design['strand'] == "-",'end'] - 5

    design['GuideSequenceWithPAM'] = design['GuideSequenceWithPAM'].str.slice(5,23)
    design[SEQCOL] = design[SEQCOL].str.slice(5,20)
    design["GuideSequenceMinusG"]=design[SEQCOL].apply(lambda x: trimG(x))
    design["MappingSequence"]=(design["GuideSequenceMinusG"]+design["RightGA"]).apply(lambda x: x[0:MapLength])
    design["CoreOligo"]=design["LeftGA"]+design["GuideSequenceMinusG"]+design["RightGA"]

    design["pool"]=PoolID
    design["OligoID"]=[PoolID+"_"+str(x+1) for x in range(len(design))]
    design["name"]=design["OligoID"]

    return design



def main(args):
    design = read_table(args.design)
    design = convertDesignToDeadGuide(design, args.PoolID)
    design.to_csv(args.outfile, sep='\t', header=True, index=False)
    try:
        writeBed(design, args.outfile + ".bed")
    except:
        print("Design completed but failed to write BED file. Possible that some entries are missing chr:start-end coordinates", sys.exc_info()[0])

if __name__ == '__main__':
    args = parseargs()
    main(args)
