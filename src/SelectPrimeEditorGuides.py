###############################################################################
## Library code for designing CRISPR screens
## Jesse Engreitz
## November 10, 2019


from __future__ import division
import pandas as pd
from pandas.io.parsers import read_table
import argparse
import numpy as np
import sys
import os
import os.path
import random

import warnings
warnings.filterwarnings("ignore")

from JuicerCRISPRiDesign import *
import DesignPrimeEditor


def parseargs():
    class formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass

    epilog = ('''
Output:
- Integrated design file
- Design files separated by subpool

''')
    parser = argparse.ArgumentParser(description='''
Select prime editor gRNAs for a set of variants.
Logic:
- Select N pegRNAs per variant per subpool (default 1)
- Optional: In cases where there is more than one variant at a location, select a majority of pegRNAs from closest spacer
- Can specify different region names of variants through the 'edits' table than was used for the design process (this is useful if you want to design for large region, then select pegRNAs within several subregions or PCR amplicons)
''',
                                     epilog=epilog,
                                     formatter_class=formatter)
    readable = argparse.FileType('r')
    
    parser.add_argument('-i', '--input', action='append', help="Input design file (output by DesignPrimeEditor.py). Specify multiple times, in order, to select pegRNAs from these design files in order of priority.")
    parser.add_argument('--guidesPerVariant', type=int, default=1, help="Maximum number of pegRNAs to select per variant.")
    parser.add_argument('--spacerChoice', default='mixed', choices=['mixed','closest'], help="closest: Always choose the closest spacer. mixed: Choose the closest spacer a majority of the time, but include more distant spacers if available.")
    parser.add_argument('--edits', help="Input file with columns: chr start end name ref alt region")

    ## Output files and options
    parser.add_argument('-o','--outfile', required=True, help="Filebase for output files")
    parser.add_argument('--splitOutputByRegion', action='store_true', default=False, help="Set to true to output separate design files for each 'region' in the input edits file")

    args = parser.parse_args()
    return(args)


def selectByPriority(pegs, n):
    selected = []
    nSelected = 0
    for priority in np.sort(pegs['priority'].unique()):
        if nSelected < n:
            curr = pegs[pegs['priority'] == priority]
            toSelect = min(n-nSelected,len(curr))
            selected.append(curr.sample(n=toSelect))
            nSelected = nSelected + toSelect
    return selected


def selectPegRNAsForVariant(pegs, guidesPerVariant=1, selectClosest=True):

    if len(pegs) == 0:   ## is this the right check?
        return []

    ## if selectClosest is True, select the pegRNA with the closest spacer among the top-priority group   
    if (pegs['priority'] == 1).sum() > 0:
        idxClosest = pegs.loc[pegs['priority'] == 1,'editPositionRelativeToNick0Based'].idxmin()
        if selectClosest:
            pegs.at[idxClosest,'priority'] = pegs.at[idxClosest,'priority'] - 0.5
        else:
            pegs.at[idxClosest,'priority'] = pegs.at[idxClosest,'priority'] + 0.5
    
    selected = selectByPriority(pegs, guidesPerVariant)

    return selected


def selectPegRNAsForPosition(pegs, guidesPerVariant=1, minFractionClosest=0.5):
    selected = []

    nVariants = len(pegs['variantName'].unique())
    nClosest = int(np.ceil(nVariants*minFractionClosest))
    selectClosest = [True]*nClosest + [False]*(nVariants-nClosest)
    random.shuffle(selectClosest)

    i = 0
    for name, group in pegs.groupby('variantName'):
        currPegs = selectPegRNAsForVariant(group, guidesPerVariant, selectClosest[i])
        selected = selected + currPegs
        i = i + 1

    return selected


def mergePegRNAsByPriority(pegs):
    '''pegs : list of pegRNA data frames, in order of priority
return : merged data frame with new column "priority"'''
    result = pegs[0]
    result['priority'] = 1

    if len(pegs) > 1:
        for df in pegs[1:]:
            df['priority'] = max(result['priority']) + 1
            result = result.append(df)
            result = result.drop_duplicates(subset=['pegRNASequence','variantName'], keep='first', ignore_index=True)

    return(result)



def selectPegRNAs(pegs, variants, guidesPerVariant=1, minFractionClosest=0.5):

    pegDf = mergePegRNAsByPriority(pegs)

    selected = []
    for name, group in variants.groupby(['chr','start','end','region']):
        group.rename(columns={"start":"variantStart","end":"variantEnd","name":"variantName","region":"newRegion"}, inplace=True)
        currPegs = group.merge(pegDf, on=['chr','variantStart','variantEnd','variantName','ref','alt'])
        currPegs['region'] = currPegs['newRegion']
        currPegs.drop('newRegion', inplace=True, axis=1)
        selected = selected + selectPegRNAsForPosition(currPegs, guidesPerVariant, minFractionClosest)

    result = pd.concat(selected)
    return result



def main(args):
    variants = read_table(args.edits)

    pegs = []
    for f in args.input:
        currPegs = read_table(f)
        currPegs['source'] = f
        pegs.append(currPegs)

    ## Warn if the region names in 'variants' are different than the region names in 'pegs'

    if args.spacerChoice == 'mixed':
        minFractionClosest = 0.5
    else:
        minFractionClosest = 1

    results = selectPegRNAs(pegs, variants, args.guidesPerVariant, minFractionClosest)

    ## Write results
    DesignPrimeEditor.writePegRNAs(results, args.outfile, args.splitOutputByRegion)

    ## Plot statistics
    DesignPrimeEditor.plotPrimeEditorStats(results, args.outfile, args.edits, args.splitOutputByRegion)


if __name__ == '__main__':
    args = parseargs()
    main(args)
