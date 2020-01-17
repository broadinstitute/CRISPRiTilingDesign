###############################################################################
## Library code for designing CRISPRi screens
## Jesse Engreitz
## January 16, 2020 - Annotate noncoding guides for gRNA design
## Tested with:  "use .python-3.5.1; source /seq/lincRNA/Ben/VENV_MIP/bin/activate"

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
from collections import OrderedDict
import itertools
import warnings
warnings.filterwarnings("ignore")


from JuicerCRISPRiDesign import *
import ReadUtils



def parseargs(required_args=True):
    class formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass

    epilog = ("")
    parser = argparse.ArgumentParser(description='Annotate noncoding guides with information about GC content, DNase counts, distance from peak summits, etc.  tabix and bedtools must be on path',
                                     epilog=epilog,
                                     formatter_class=formatter)
    readable = argparse.FileType('r')
    
    parser.add_argument('--input', required=required_args, help="preDesign.bed file from GetTileGuides.py")
    parser.add_argument('--outfile', required=required_args)
    parser.add_argument('--dhsList', help="File with columns: CellType BAM MACS2Summits, one row per cell type. If providing this, do not need --summits --bam --cellType")
    parser.add_argument('--summits', help="Name of column with guide sequence")
    parser.add_argument('--enhancers', required=required_args, help="EnhancerList.txt from ABC neighborhoods directory")
    parser.add_argument('--cellType', default="CellType1", help="Name of cell type (for annotation purposes only")
    parser.add_argument('--bam', help="BAM file for ATAC or DNase. .tagAlign.gz also okay")
    parser.add_argument('--PAM', default="NGG", help="PAM sequence for gRNAs")
    parser.add_argument('--genomeSizes', default='/seq/lincRNA/data/hg19/sizes', help="File with genome sizes")
    args = parser.parse_args()
    return(args)



def setupCellTypeParams(dhsList, summits, cellType, bam, enhancers):
    if dhsList is None:
        params = OrderedDict([ 
            ('CellType', [cellType]),
            ('BAM', [bam]),
            ('MACS2Summits',  [summits]),
            ('Enhancers', [enhancers]) ] )
        params = pd.DataFrame.from_dict(params)
    else:
        params = read_table(dhsList)
    return params


def readGuides(INFILE, PAM):
    # count the number of columns
    with open(INFILE) as f:
        reader=csv.reader(f, delimiter="\t", skipinitialspace=True)
        first_row = next(reader)
        num_cols = len(first_row)

        print(str(num_cols) + " columns in " + INFILE)

    # if there is a name column, include.
    if num_cols == 18:
        header=["chr", "start", "end", "locus", "score", "strand", "start.", "end.", "zero", "1", "20","zero.", "GuideSequenceWithPAM", "guideSet","peakChr", "peakStart", "peakEnd", "peakName"]
        guides=read_table(INFILE,names=header)
    
    elif num_cols == 17:
        header=["chr", "start", "end", "locus", "score", "strand", "start.", "end.", "zero", "1", "20","zero.", "GuideSequenceWithPAM", "guideSet","peakChr", "peakStart", "peakEnd"]
        guides=read_table(INFILE,names=header)

    elif num_cols == 14:
        header=["chr", "start", "end", "locus", "score", "strand", "start.", "end.", "zero", "1", "20","zero.", "GuideSequenceWithPAM", "guideSet"]
        guides=read_table(INFILE,names=header)
    
    else:
        print("Required 14, 17, or 18 columns")
        assert False

    guides["GuideSequence"]=guides["GuideSequenceWithPAM"].apply(lambda x: trimPAM(x, PAM))
    guides["GuideSequenceMinusG"]=guides["GuideSequence"].apply(lambda x: trimG(x))
    guides["GuideSequencePlusG"]=guides["GuideSequenceMinusG"].apply(lambda x: "G"+x)
    return guides


def getMaxRepetition(string, ch):
    if not ch in string:
        return 0
    z = [(x[0], len(list(x[1]))) for x in itertools.groupby(list(string))]
    return max([y for y in z if y[0] == ch], key=lambda x:x[1])[1]


def getMaxEndRepetition(string, ch):
    if not string[-1] == ch:
        return 0
    z = [(x[0], len(list(x[1]))) for x in itertools.groupby(list(string))]
    return z[-1][1]


def addSequenceFeatures(guides):
    guides["polyA"] = guides["GuideSequencePlusG"].apply(lambda x: getMaxRepetition(x,"A"))
    guides["polyC"] = guides["GuideSequencePlusG"].apply(lambda x: getMaxRepetition(x,"C"))
    guides["polyG"] = guides["GuideSequencePlusG"].apply(lambda x: getMaxRepetition(x,"G"))
    guides["polyT"] = guides["GuideSequencePlusG"].apply(lambda x: getMaxRepetition(x,"T"))
    guides["polyV"] = guides[['polyA','polyC','polyG']].max(axis=1)
    guides["GC"] = guides["GuideSequencePlusG"].apply(GCContent)
    guides["endTCount"] = guides["GuideSequencePlusG"].apply(lambda x: getMaxEndRepetition(x,"T"))
    guides["countA"] = guides["GuideSequencePlusG"].apply(lambda x: list(x).count('A'))
    guides["countC"] = guides["GuideSequencePlusG"].apply(lambda x: list(x).count('C'))
    guides["countG"] = guides["GuideSequencePlusG"].apply(lambda x: list(x).count('G'))
    guides["countT"] = guides["GuideSequencePlusG"].apply(lambda x: list(x).count('T'))   
    return guides


def addDHSCounts(guides, guideBed, cellType, bam, genome_sizes):
    bamTotal = ReadUtils.count_total(bam)

    guides = getDhsRpm(guides, guideBed, bam, guideBed + '.tmp.count20.bed', genome_sizes, cellType + '.dhsRpm20', bamTotal, cellType)

    bed100 = guideBed + '.bed100.bed'
    command = "bedtools slop -i {guideBed} -b 40 -g {genome_sizes} > {bed100}".format(**locals())
    (stdout, err) = ReadUtils.runCommand(command)
    guides = getDhsRpm(guides, bed100, bam, guideBed + '.tmp.count100.bed', genome_sizes, cellType + '.dhsRpm100', bamTotal, cellType)
    ReadUtils.runCommand("rm " + bed100)
    
    bed500 = guideBed + '.bed500.bed'
    command = "bedtools slop -i {guideBed} -b 240 -g {genome_sizes} > {bed500}".format(**locals())
    (stdout, err) = ReadUtils.runCommand(command)
    guides = getDhsRpm(guides, bed500, bam, guideBed + '.tmp.count500.bed', genome_sizes, cellType + '.dhsRpm500', bamTotal, cellType)
    ReadUtils.runCommand("rm " + bed500)

    return guides


def getDhsRpm(guides, guideBed, bam, tmpout, genome_sizes, colname, bamTotal, cellType):
    ## run count reads and extract the answer
    ReadUtils.run_count_reads(bam, tmpout, guideBed, genome_sizes, use_fast_count=True, count5pEnds=True)
    data = ReadUtils.read_bedgraph(tmpout)
    guides[colname] = data['score'] / (bamTotal / 1000000)
    col = cellType + '.nhbd.DHS.RPM'
    #guides.loc[guides[col] == '.', col] = 0
    guides[col] = [0 if x == '.' else float(x) for x in guides[col]]
    guides[colname + '.fractionOfRegion'] = guides[colname] / guides[col]
    #data = data.rename({'score':colname})
    #guides = guides.merge(data, on=['chr','start','end'])
    ReadUtils.runCommand("rm " + tmpout)
    return guides


def annotateWithEnhancers(guides, guideFile, eFile, cellType):
    enhancers = read_table(eFile)
    selectedData = enhancers[['chr','start','end','name','class','DHS.RPM','DHS.RPKM']]
    enhancerFile = guideFile + ".tmp.enhancers.txt"
    selectedData.to_csv(enhancerFile, sep='\t', header=False, index=False)

    tmpFile = guideFile + ".tmp.guides.txt"
    guides.to_csv(tmpFile, sep='\t', header=False, index=False)

    overlapFile = guideFile + ".tmp.enhancerOverlap.txt"
    command = "bedtools intersect -f 1.0 -a {tmpFile} -b {enhancerFile} -wao > {overlapFile}".format(**locals())
    ReadUtils.runCommand(command)
    result = read_table(overlapFile)
    result.columns = guides.columns.to_list() + [cellType + '.nhbd.' + x for x in selectedData.columns.to_list()] + [cellType + '.nhbd.overlapBases']

    ReadUtils.runCommand("rm " + tmpFile + " " + enhancerFile + " " + overlapFile)
    return(result)


def main(args):
    cellTypes = setupCellTypeParams(args.dhsList, args.summits, args.cellType, args.bam, args.enhancers)

    guides = readGuides(args.input, args.PAM)
    guides = addSequenceFeatures(guides)

    for index,row in cellTypes.iterrows():
        guides = annotateWithEnhancers(guides, args.input, row.Enhancers, row.CellType)
        guides = addDHSCounts(guides, args.input, row.CellType, row.BAM, args.genomeSizes)
        #guides = addSummitDistance(guides, row.MACS2Summits)  ## Not yet implemented

    guides.to_csv(args.outfile, sep='\t', header=True, index=False, float_format='%.5f')


if __name__ == '__main__':
    args = parseargs()
    main(args)
