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
    parser.add_argument('--nGuidesPerElement', required=required_args, type=int,    help="Number of guides to choose per element (in the guideSet column)")
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

    args = parser.parse_args()
    return(args)



def makeOligos(STRING, side, vector="sgopti"):
    '''makes the single guide cloning oligos for a single guide.
    specify side (top or bottom) and vector (sgOpti or pZB)'''
    if vector.lower()=="sgopti":
        STRING=prependGifNecessary(STRING)
        if side.lower()=="top":
            return "CACC"+STRING
        elif side.lower()=="bot":
            SEQ=Seq(STRING)
            RC=SEQ.reverse_complement()
            return "AAAC"+str(RC)
        

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
def selectNguidesPerElement(df, NGuides, columnName="peakName", minGuides=0, method="score", includeGuides=None):

    ''' Annotates with OffTargetScore (int(score) and Quality Score. Selects N guides per element'''
    df["OffTargetScore"]=df["score"].apply(lambda x: int(x)) 
    df["QualityScore"]=df.SSC*100+df.OffTargetScore # for now just get em close and add em
    
    # For each element, select the top N elements based on Quality score
    chosen=pd.DataFrame()
    minGuidesInAnElement=len(df)
    
    for i in set(df[columnName]):
        tSet=df[df[columnName]==i]
        
        if len(tSet)>=minGuides:
            
            if len(tSet)<minGuidesInAnElement:
                minGuidesInAnElement=len(tSet)

            if includeGuides is not None:
                currInclude = pd.merge(tSet, includeGuides)
            else:
                currInclude = tSet.head(0)
            
            nToAdd = max(0,NGuides-len(currInclude))
            if method == 'score':
                chosen=pd.concat([chosen, currInclude, tSet.sort_values("QualityScore", ascending=False).head(nToAdd)])
            elif method == 'even':
                chosen=pd.concat([chosen, currInclude, tSet.iloc[getEvenlySpacedIndices(tSet,nToAdd),]])
            else:
                raise ValueError("Guide selection method " + method + " is not supported.")

    print(minGuidesInAnElement, "Minimum guides per element")
    return chosen



def makeDesignFile(pool, PoolID):
    DESIGNCOLS=["chr", "start", "end","name", "score", "strand", "GuideSequence", "GuideSequenceMinusG",
                "MappingSequence", "OffTargetScore", "target", "subpool", "OligoID"]
    # Add in an oligo ID based on name of oligo pool
    design = pool.copy()
    design["pool"]=PoolID
    design["OligoID"]=[PoolID+"_"+str(x+1) for x in range(len(design))]
    design["name"]=design["OligoID"]
    design["start"]=design["start"].apply(lambda x: smartInt(x))
    design["end"]=design["end"].apply(lambda x: smartInt(x))    
    return design


def getNegativeControlGuides(negCtrlList, nCtrls):
    negCtrls = read_table(negCtrlList)
    negCtrls["OffTargetScore"]=200
    negCtrls["QualityScore"]=200

    if (nCtrls > len(negCtrls)):
        raise ValueError("Asking for more negative controls than are provided in the negCtrlList")
    negCtrls = negCtrls.head(nCtrls)
    return negCtrls


def loadForceGuides(file):
    if file is not None:
        includeGuides = read_table(args.forceGuides, header=None) 
        includeGuides.columns = ['locus']
    else:
        includeGuides = None    
    return includeGuides



def main(args):
    guides = read_table(args.input)
    guides = fillInDesign(guides, SEQCOL=args.seqCol, TRIMPAM=(not args.noPAM))
    includeGuides = loadForceGuides(args.forceGuides)
    selected = selectNguidesPerElement(guides, int(args.nGuidesPerElement), columnName="guideSet", minGuides=0, method=args.selectMethod, includeGuides=includeGuides)

    if (args.nCtrls > 0):
        negCtrls = getNegativeControlGuides(args.negCtrlList, args.nCtrls)
        combined = pd.concat([selected, negCtrls])
    else:
        combined = selected

    if (args.nSafeCtrls > 0):
        safeCtrls = getNegativeControlGuides(args.safeCtrlList, args.nSafeCtrls)
        combined = pd.concat([combined, safeCtrls])

    design = addGibsonArms(combined, args.vectorDesigns, args.vector)
    design = makeDesignFile(design, args.PoolID)
    design.to_csv(os.path.join(args.outdir, args.PoolID + ".design.txt"), sep='\t', header=True, index=False)
    try:
        writeBed(design, os.path.join(args.outdir, args.PoolID + ".design.bed"))
    except:
        print("Failed to write BED file.", sys.exc_info()[0])


if __name__ == '__main__':
    args = parseargs()
    main(args)
