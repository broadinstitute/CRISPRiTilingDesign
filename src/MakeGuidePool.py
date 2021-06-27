###############################################################################
## Library code for designing CRISPRi screens
## Jesse Engreitz
## November 10, 2019
## Based on Charlie's gRNA design
## Tested with:  "use .python-3.5.1; source /seq/lincRNA/Ben/VENV_MIP/bin/activate"

## TODO: Factor out helper code


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

    epilog = ("")
    parser = argparse.ArgumentParser(description='Create oligos for a CRISPRi pool or subpool. This is based on Charlie\'s design script',
                                     epilog=epilog,
                                     formatter_class=formatter)
    readable = argparse.FileType('r')
    
    parser.add_argument('--input', required=required_args, help="preDesign.bed file from GetTileGuides.py")
    parser.add_argument('--outdir', required=required_args)
    parser.add_argument('--nGuidesPerElement', required=required_args, type=int,    help="Number of guides to choose per element (in the guideSet column). Set to 0 to select all.")
    parser.add_argument('--negCtrlList', default="data/Weissman1000.negative_control.20bp.design.txt", help="Formatted guide file of negative controls.")
    parser.add_argument('--nCtrls', default=0, type=int, help="Number of negative control gRNAs to include")
    parser.add_argument('--safeCtrlList', default="data/Tycko2019SafeTargeting.txt", help="Formatted guide file for safe-targeting controls.")
    parser.add_argument('--nSafeCtrls', default=0, type=int, help="Number of safe-targeting gRNAs to include")
    parser.add_argument('--vectorDesigns', default="data/CloningDesigns.txt", help="Master file with gibson arm sequences for various plasmid designs")
    parser.add_argument('--vector', default='sgOpti', help="Name of vector to index into vectorDesigns")
    parser.add_argument('--PoolID', default="MyPool", help="Unique name of pool or subpool for naming oligos - e.g. 191110_GATA1")
    parser.add_argument('--seqCol', default='GuideSequenceWithPAM', help="Name of column with guide sequence")
    parser.add_argument('--noPAM', action="store_true", default=False, help="Whether to trim PAM from guide sequence in seqCol")
    parser.add_argument('--PAM', default="NGG", help="Only NGG PAM is currently supported")
    parser.add_argument('--selectMethod', default='score', choices=['score','even'], help="Method to select guides within elements")
    parser.add_argument('--forceGuides', help="Pass a file containing one guide name per line to force selection of these gRNAs (useful for including gRNAs that we know work)")
    parser.add_argument('--trimElements', default=0, type=int, help="Remove guides corresponding to elements that have fewer than this many guides (e.g., because not enough guides existed to choose from, repeats, etc.).")
    parser.add_argument('--excludeRestrictionSites', default="", help="Comma-delimited list of subsequences (e.g. restriction enzyme sites) to exclude gRNAs.")
    parser.add_argument('--barcodes', help="File with one unique barcode sequence per line, e.g. for HyPR screens. Will substitute these sequences in place of '[NNNNN]' sequence from vector design file.  Errors out if there are more unique guides than barcodes.")

    args = parser.parse_args()
    return(args)


def fillInDesign(df, SEQCOL="seq", sgOligos=False, TRIMPAM=True): 
    if TRIMPAM:
        df["seq"]=df[SEQCOL].apply(lambda x: trimPAM(x))
        SEQCOL="seq"
    df["GuideSequence"]=df[SEQCOL].apply(lambda x: prependGifNecessary(x))
    df["GuideSequenceMinusG"]=df["GuideSequence"].apply(lambda x: trimG(x))
    
    # make oligos for single guide cloning
    if sgOligos:
        df["Top"]=df["GuideSequenceMinusG"].apply(lambda x: makeOligos(x, "top"))
        df["Bot"]=df["GuideSequenceMinusG"].apply(lambda x: makeOligos(x, "bot"))

    cols=['chr','start','end','locus','strand','guideSet']
    for col in cols:
        if not col in df.columns:
            df[col] = ''

    num=['score','SSC']
    for col in num:
        if not col in df.columns:
            df[col] = 0

    return df


def addHandles(df, handleDict, SUBPOOLCOL="subpool"):
    ''' Adds PCR handles for PCR1 '''
    ## TODO:  Not used by this script.  Move to "CombinePools.py" script

    FWDDICT=readTDFtodict(handleDict, KEYCOL="subpool", VALCOL="FwdPrimer")
    REVDICT=readTDFtodict(handleDict, KEYCOL="subpool", VALCOL="RevPrimer")    

    df["FwdPrimer"]=df["subpool"].apply(lambda x: FWDDICT[x])
    df["RevPrimer"]=df["subpool"].apply(lambda x: REVDICT[x])
    
    # RC Rev primer to put onto oligo
    df["RevPrimerOligo"]=df["RevPrimer"].apply(lambda x: str(Seq(x).reverse_complement()))

    # add the PCR handles to the CoreOligo
    '''
    oligoSeq=[""]*len(df)
    for row in range(len(df)):
        tRow=df.irow(row)
        oligoSeq[row]=tRow["FwdPrimer"]+tRow["CoreOligo"]+tRow["RevPrimerOligo"]
    df["OligoSequence"]=oligoSeq
    '''
    df["OligoSequence"]=df["FwdPrimer"]+df["CoreOligo"]+df["RevPrimerOligo"]

    return df


def addGibsonArms(df, vectorDesignFile, vector, MapLength=21):
    vectorDesigns = read_table(vectorDesignFile)
    if not vector in vectorDesigns['CloningDesign'].values:
        raise ValueError("Chosen vector design '"+vector+"' is not in the Vector Design file.")
    leftGA = vectorDesigns.loc[vectorDesigns['CloningDesign'] == vector].LeftGibsonArm.item()
    rightGA = vectorDesigns.loc[vectorDesigns['CloningDesign'] == vector].RightGibsonArm.item()
    df["LeftGA"]=leftGA
    df["RightGA"]=rightGA
    df["CoreOligo"]=leftGA+df["GuideSequenceMinusG"]+rightGA
    df["MappingSequence"]=(df["GuideSequenceMinusG"]+df["RightGA"]).apply(lambda x: x[0:MapLength])
    return df
    
    
# select guides per element
def selectNguidesPerElement(df, NGuides, columnName="peakName", minGuides=0, method="score", includeGuides=None, trimElements=0):

    ''' Annotates with OffTargetScore (int(score) and Quality Score. Selects N guides per element'''
    df["OffTargetScore"]=df["score"].apply(lambda x: int(x)) 
    df["QualityScore"]=df.SSC*100+df.OffTargetScore # for now just get em close and add em
    
    # For each element, select the top N elements based on Quality score
    chosen=df.head(0)
    minGuidesInAnElement=len(df)

    # Force the inclusion of certain guides, regardless of any other considerations
    if includeGuides is not None:
        toAdd = pd.merge(df, includeGuides)
        chosen=pd.concat([chosen, toAdd])
        print("Forcing inclusion of " + str(len(toAdd)) + " gRNAs (out of " + str(len(includeGuides)) + " in the --forceGuides file)")

    for i in sorted(set(df[columnName]), reverse=True):
        tSet=df[df[columnName]==i]
        tSet=tSet.sort_values(by='start')

        if len(tSet)>=minGuides:
            
            if len(tSet)<minGuidesInAnElement:
                minGuidesInAnElement=len(tSet)

            toAdd = tSet.head(0)

            if includeGuides is not None:
                toAdd = pd.merge(tSet, includeGuides)

            if (len(tSet)<trimElements) and (len(toAdd) == 0):
                ## Skip elements with too few guides and where we weren't try to force include the guides
                print("Skipping " + str(i) + " because there are too few guides")
                continue

            # If guide has been already chosen for a different (overlapping) element, pre-select these guides
            #import pdb; pdb.set_trace()
            toAdd = pd.concat([toAdd, tSet[tSet['locus'].isin(chosen['locus'])]])

            nToAdd = max(0,NGuides-len(toAdd))
            tSet = tSet[~tSet['locus'].isin(toAdd['locus'])]

            if method == 'score':
                toAdd = pd.concat([toAdd, tSet.sort_values("QualityScore", ascending=False).head(nToAdd)])

            elif method == 'even':
                toAdd = pd.concat([toAdd, tSet.iloc[getEvenlySpacedIndices(tSet,nToAdd),]])
            else:
                raise ValueError("Guide selection method " + method + " is not supported.")
            
            chosen=pd.concat([chosen, toAdd])
    
    ## Now that we've made selections per guide, go back and pull in info from all guide-target pairs where
    ##  that guide has been selected for at least one target
    final = pd.merge(df, chosen[['GuideSequenceWithPAM']].drop_duplicates())

    print(minGuidesInAnElement, "Minimum guides per element")
    print("Selected " + str(len(chosen[['GuideSequenceWithPAM']].drop_duplicates())) + " unique guides in selectNguidesPerElement.")
    return final



def makeDesignFile(pool, PoolID):
    DESIGNCOLS=["chr", "start", "end","name", "score", "strand", "GuideSequence", "GuideSequenceMinusG",
                "MappingSequence", "OligoID"]
    # Add in an oligo ID based on name of oligo pool
    design = pool.copy()
    design["pool"]=PoolID

    uniqOligos=design['GuideSequenceMinusG'].unique()
    ids = pd.DataFrame({
        "OligoID": [PoolID+"_"+str(x+1) for x in range(len(uniqOligos))],
        "GuideSequenceMinusG": uniqOligos
        })
    design = pd.merge(design, ids, on='GuideSequenceMinusG') 
#    design["OligoID"]=[PoolID+"_"+str(x+1) for x in range(len(design))]
    design["name"]=design["OligoID"]
    design = design.astype({'start': pd.Int64Dtype(), 'end': pd.Int64Dtype()})  ## Int64Dtype required to allow NA values for negative_control guides

    ## Reorder columns
    cols = DESIGNCOLS + design.columns[~design.columns.isin(DESIGNCOLS)].tolist()
    design = design[cols]
    return design


def addBarcodeToOligo(oligo, barcode, toReplace):
    if (oligo.find(toReplace) == -1):
        raise ValueError("Failed adding barcode sequence to oligo: Expected to find '" + toReplace +"' in the oligo sequence: " + oligo)
    oligo = oligo.replace(toReplace, str(barcode))
    return oligo


def addBarcodes(design, barcodeFile, toReplace="[NNNNN]"):
    '''
    Add HyPR barcodes to the final CoreOligo sequence.
    '''
    barcodes = read_table(barcodeFile, header=None)
    barcodes.columns = ['BarcodeSequence']

    if len(barcodes) != len(barcodes.BarcodeSequence.unique()):
        raise ValueError("Not allowed: duplicate barcodes in " + barcodeFile)

    uniqOligos = design['CoreOligo'].unique()

    if len(uniqOligos) > len(barcodes):
        raise ValueError("Number of barcodes provided in " + barcodeFile + " (" + str(len(barcodes)) + ") is less than number of unique oligos (" + str(len(uniqOligos)) + ").")

    barcodes = barcodes.head(len(uniqOligos))
    barcodes['CoreOligo'] = uniqOligos

    cols = design.columns.to_list() + ['BarcodeSequence']
    design = pd.merge(design, barcodes, on='CoreOligo')
    design = design[cols]

    for idx,row in design.iterrows():
        design.loc[idx,'CoreOligo'] = addBarcodeToOligo(row['CoreOligo'], row['BarcodeSequence'], toReplace)

    return design


def getNegativeControlGuides(negCtrlList, nCtrls, reSites, vectorDesigns, vector):
    negCtrls = read_table(negCtrlList)
    negCtrls["OffTargetScore"]=200
    negCtrls["QualityScore"]=200
    negCtrls = addGibsonArms(negCtrls, vectorDesigns, vector)
    negCtrls = excludeRestrictionSites(negCtrls, reSites)

    if (nCtrls > len(negCtrls)):
        raise ValueError("Asking for more negative controls (" + str(nCtrls) + ") than are provided in the negCtrlList (" + str(len(negCtrls)) + ")")
    negCtrls = negCtrls.head(nCtrls)
    return negCtrls


def loadForceGuides(file):
    if file is not None:
        includeGuides = read_table(args.forceGuides, header=None) 
        includeGuides.columns = ['locus']
    else:
        includeGuides = None    
    return includeGuides


def excludeRestrictionSites(guides, excludeCommaList):
    #import pdb; pdb.set_trace()
    if excludeCommaList != '':
        for curr in excludeCommaList.split(','):
            ## Only look at guide sequence +/- 10bp on either side ... because the GA might have the RE sequence by design
            middle = pd.Series([row["LeftGA"][-10:] + row["GuideSequenceMinusG"] + row["RightGA"][10:] for key,row in guides.iterrows()])
            guides = guides[~middle.str.contains(curr, case=False)]
            rc = str(Seq(curr).reverse_complement())
            guides = guides[~middle.str.contains(rc, case=False)]
    return guides


def main(args):
    guides = read_table(args.input)
    guides = fillInDesign(guides, SEQCOL=args.seqCol, TRIMPAM=(not args.noPAM))
    guides = addGibsonArms(guides, args.vectorDesigns, args.vector)
    guides = excludeRestrictionSites(guides, args.excludeRestrictionSites)

    includeGuides = loadForceGuides(args.forceGuides)

    if args.nGuidesPerElement > 0:
        selected = selectNguidesPerElement(guides, args.nGuidesPerElement, columnName="guideSet", minGuides=0, method=args.selectMethod, includeGuides=includeGuides, trimElements=args.trimElements)
    else:
        print("Keeping all guides because --nGuidesPerElement is 0.")
        selected = guides

    if (args.nCtrls > 0):
        negCtrls = getNegativeControlGuides(args.negCtrlList, args.nCtrls, args.excludeRestrictionSites, args.vectorDesigns, args.vector)
        combined = pd.concat([selected, negCtrls])
    else:
        combined = selected

    if (args.nSafeCtrls > 0):
        safeCtrls = getNegativeControlGuides(args.safeCtrlList, args.nSafeCtrls, args.excludeRestrictionSites, args.vectorDesigns, args.vector)
        combined = pd.concat([combined, safeCtrls])

    design = makeDesignFile(combined, args.PoolID)

    ## Add HyPR barcodes if applicable
    if args.barcodes is not None:
        design = addBarcodes(design, args.barcodes)

    design.to_csv(os.path.join(args.outdir, args.PoolID + ".design.txt"), sep='\t', header=True, index=False)
    try:
        writeBed(design[['chr','start','end','name']].drop_duplicates(), os.path.join(args.outdir, args.PoolID + ".design.bed"))
    except:
        print("Failed to write BED file.", sys.exc_info()[0])


if __name__ == '__main__':
    args = parseargs()
    main(args)
