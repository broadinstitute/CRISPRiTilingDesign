#! /usr/bin/env python

# USAGE: python ./GetTileGuides.py
# Must use Python-2.7

# ------------------------------------------- #
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

def limitToLocus(df, (CHR, START, END)):
    df=df[df["chr"]==CHR]
    df=df[df["start"]>=int(START)]
    df=df[df["end"]<=int(END)]
    return a

def GCContent(x):
    x=x.upper()
    return (x.count("G")+x.count("C"))/len(x)

def trimPAM(seq):
    assert seq[-2:]=="GG"
    return seq[0:-3]
def trimG(seq):
    seq=seq.upper()
    if seq[0]=="G":
        return seq[1:]
    else:
        return seq

def writeBed(df, OUTFILE, score=None,negative=False):
    df=df.copy()
    df["start"]=np.vectorize(str)(np.vectorize(int)(df["start"]))
    df["end"]=np.vectorize(str)(np.vectorize(int)(df["end"]))
    
    cols=["chr", "start", "end"]
    if score:
        cols=["chr", "start", "end", score]
    if negative:
        df[score]=-df[score]
    df[cols].to_csv(OUTFILE, sep='\t', header=False, index=False)

def drawDensity(vec, label=None, fill=False):
    x=np.arange(min(vec), max(vec), (max(vec)-min(vec))/1000)
    density=stats.kde.gaussian_kde(vec)
    if label:
        plt.plot(x, density(x), label=label)
    else:
        plt.plot(x, density(x))
    if fill:
        plt.fill_between(x, 0, density(x))

#----------------------------------------------------------------#
# SSC 
SSC19Loc="/seq/lincRNA/cfulco/bin/SSC/SSC0.1/matrix/human_CRISPRi_19bp.matrix"
SSC20Loc="/seq/lincRNA/cfulco/bin/SSC/SSC0.1/matrix/human_CRISPRi_20bp.matrix"
SSC21Loc="/seq/lincRNA/cfulco/bin/SSC/SSC0.1/matrix/human_CRISPRi_21bp.matrix"

def scoreByMatrix(STRING, MATRIX):
    STRING=STRING.upper()
    assert set(STRING)<=set(MATRIX[0]) # require same alphabets
    score=0
    for i in range(0,len(STRING)):
        score=score+MATRIX[i][STRING[i]]
    return score

def readSSCmatrix(MATRIXFILE):
    '''
    Read the .matrix file from SSC:
    intercept 
    A C G T
    0 -1 1 0.5 ... '''

    matrix=read_table(MATRIXFILE, header=1)
    LOD=[]
    for i in range(0, len(matrix)):
        LOD.append({matrix.columns[0]:matrix.irow(i)[matrix.columns[0]],
            matrix.columns[1]:matrix.irow(i)[matrix.columns[1]],
            matrix.columns[2]:matrix.irow(i)[matrix.columns[2]],
            matrix.columns[3]:matrix.irow(i)[matrix.columns[3]]})
    return LOD

def scoreByCorrectMatrix(STRING, LISTofMATRIX):
    for MATRIX in LISTofMATRIX:
        if len(STRING)==len(MATRIX):
            return scoreByMatrix(STRING, MATRIX)
    return scoreByMatrix(STRING[-len(MATRIX):], MATRIX) # score based on right of guide

# SSC data
[SSC19, SSC20, SSC21]=[readSSCmatrix(SSC19Loc),
                       readSSCmatrix(SSC20Loc),
                       readSSCmatrix(SSC21Loc)]
LISTofMATRIX=[SSC19, SSC20, SSC21]

# --------------------------------------------------------------------------------------- #
# Filtering Function
# --------------------------------------------------------------------------------------- #
def GetTileGuides(INFILE, OUTPREFIX, LOCUS, SSCMIN, OTMIN, POLYT, POLYV, ENDTCount, MAXG, MINGC, MAXGC, MINSTARTDIST, PLOT):
    
    # count the number of columns
    with open(INFILE) as f:
        reader=csv.reader(f, delimiter="\t", skipinitialspace=True)
        first_row = next(reader)
        num_cols = len(first_row)

        print num_cols, "number of columns in ", INFILE

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
        print "Required 14, 17, or 18 columns"
        assert False

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Prepare outfiles
    OUTFILE=OUTPREFIX+".getGuides.txt" # data about guide counts
    outfile=open(OUTFILE, 'w')

    OUTBED=OUTPREFIX+".preDesign.bed" # the file with guides.
    
    OUTPDF=OUTPREFIX+".getGuides.pdf" # plots

    # write the command with all the options
    outfile.write(" ".join([str(x) for x in sys.argv[:]])+ " \n")
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

    # will edit dataframe a
    df=guides.copy()

    if LOCUS[0] != "All" and LOCUS[0] is not None:
        df=limitToLocus(df, LOCUS)
    
    outfile.write("within region of interest: "+str(len(df.drop_duplicates("GuideSequenceWithPAM")))+"\n")
    
    # fill in other guide sequences       
    df["GuideSequence"]=df["GuideSequenceWithPAM"].apply(lambda x: trimPAM(x))
    df["GuideSequenceMinusG"]=df["GuideSequence"].apply(lambda x: trimG(x))

    #import pdb; pdb.set_trace()
    ## TODO JME: Factor out these guide filters so they can be easily applied to the negative control guides
    df=df[df["score"]>=OTMIN]
    if len(df) == 0:
        raise RuntimeException("No guides passed score filter")

    outfile.write("OTScore >= "+str(OTMIN)+" "+ str(len(df["GuideSequenceMinusG"].drop_duplicates()))+"\n")

    df=df[df["GuideSequenceMinusG"].apply(lambda x: "T"*POLYT not in "G"+x)]
    print str(1) + " " + str(len(df[df['guideSet'] == "K562-Roadmap-30"]))
    df=df[df["GuideSequenceMinusG"].apply(lambda x: "A"*POLYV not in "G"+x)]
    print str(2) + " " + str(len(df[df['guideSet'] == "K562-Roadmap-30"]))
    df=df[df["GuideSequenceMinusG"].apply(lambda x: "G"*POLYV not in "G"+x)]
    print str(3) + " " + str(len(df[df['guideSet'] == "K562-Roadmap-30"]))
    df=df[df["GuideSequenceMinusG"].apply(lambda x: "C"*POLYV not in "G"+x)]
    print str(4) + " " + str(len(df[df['guideSet'] == "K562-Roadmap-30"]))
    outfile.write("polyBaseFilters "+str(len(df["GuideSequenceMinusG"].drop_duplicates()))+"\n")

    df=df[df["GuideSequenceMinusG"].apply(lambda x: x[-ENDTCount:]!="T"*ENDTCount)]
    outfile.write("No guides ending with "+str(ENDTCount)+" Ts: "+str(len(df["GuideSequenceMinusG"].drop_duplicates()))+"\n")
    print str(5) + " " + str(len(df[df['guideSet'] == "K562-Roadmap-30"]))

    df=df[df["GuideSequenceMinusG"].apply(lambda x: x.count("G"))<MAXG]
    outfile.write("10 G Filter "+str(len(df["GuideSequenceMinusG"].drop_duplicates()))+"\n")
    print str(6) + " " + str(len(df[df['guideSet'] == "K562-Roadmap-30"]))

    df=df[df["GuideSequenceMinusG"].apply(lambda x: GCContent("G"+x))>MINGC]
    df=df[df["GuideSequenceMinusG"].apply(lambda x: GCContent("G"+x))<MAXGC]
    outfile.write("GCFilter "+str(len(df["GuideSequenceMinusG"].drop_duplicates()))+"\n")
    print str(7) + " " + str(len(df[df['guideSet'] == "K562-Roadmap-30"]))

    df["SSC"]=df["GuideSequence"].apply(lambda x: scoreByCorrectMatrix(x, LISTofMATRIX))
    print str(8) + " " + str(len(df[df['guideSet'] == "K562-Roadmap-30"]))

    # filter for very low SSC
    preF=len(df["GuideSequenceMinusG"].drop_duplicates())
    df=df[df["SSC"]>SSCMIN]
    postF=len(df["GuideSequenceMinusG"].drop_duplicates())
    outfile.write("SSC filter "+str(postF)+". fraction: "+str(1-postF/preF)+"\n")
    print str(9) + " " + str(len(df[df['guideSet'] == "K562-Roadmap-30"]))

    # - - - - - - - - - - - - - - - - - - - - - - - - #
    # Shrink array by removing very close guides, keeping the guide with higher SSC
    # for now ignore which peak it comes from
    
    b=df.sort(["chr", "start"])
    
    # if dont need to shrink (i.e. MINSTARTDIST<1), skip this step
    
    if MINSTARTDIST>0:
        b=b.drop_duplicates("GuideSequenceWithPAM")
        legal=[]

        legal.append(b.irow(0))
        for row in range(1,len(b)):
    
                # check if current guide overlaps last legal guide. Keep if it doesn't
                if abs(b.irow(row)["start"]-legal[-1]["start"])>=MINSTARTDIST:
                    legal.append(b.irow(row))
    
                # if there is an overlap,
                # if this guide is better REPLACE the current recent legal
                elif b.irow(row)["SSC"]>legal[-1]["SSC"]:
                    legal[-1]=b.irow(row)
    
                # if this guide is worse, do nothing
                else:
                    pass
    
        l=pd.DataFrame(legal)

        #only write thin bedgraph if doing the thinning
        writeBed(l, OUTPREFIX+".getGuides.thin.bedgraph", score="SSC")

    elif MINSTARTDIST<1:
        l=b.copy()
    else:
        print "ERROR IN MINSTARTDIST"

    outfile.write("Guides Passing Filters: "+str(len(b))+"\n")
    outfile.write("Guides With Thinning: "+str(len(l))+"\n")
    outfile.close()
    
    # - - - - - - - - - - - - - - - - - - - - - - - - #
    # write all the guides
    
    if num_cols == 18:
        cols=["chr", "start", "end", "locus", "score", "strand", "GuideSequenceWithPAM", "guideSet", "SSC", "peakName"]
    elif (num_cols == 17) or (num_cols == 14):
        cols=["chr", "start", "end", "locus", "score", "strand", "GuideSequenceWithPAM", "guideSet", "SSC"]
    

    l[cols].to_csv(OUTBED, sep="\t", index=False, header=True, quote=False)

    # write bedgraphs of thin and not thinned libraries
    #writeBed(l, OUTPREFIX+".getGuides.thin.bedgraph", score="SSC")
    writeBed(b, OUTPREFIX+".getGuides.all.bedgraph", score="SSC")

    # - - - - - - - - - - - - - - - - - - - - - - - - #
    # plot
    if PLOT:
        plt.subplot(2,1,1)
        drawDensity(l["SSC"], "Legal")
        drawDensity(b["SSC"], "all")
        plt.legend(loc="best")
        plt.xlabel("SSC")
    
        # adjacent guide distance
        plt.subplot(2,1,2)
        bDiff=[(b["start"].irow(x+1)-b["start"].irow(x)) for x in range(len(b)-1)]
        lDiff=[(l["start"].irow(x+1)-l["start"].irow(x)) for x in range(len(l)-1)]
    
        plt.hist(lDiff, normed=True, cumulative=True, histtype='step', label="legal", bins=700000)
        plt.hist(bDiff, normed=True, cumulative=True, histtype='step', label="all", bins=700000)
        plt.legend(loc="best")
        plt.xlim([-1,70])
        plt.xlabel("AdjacentGuideDistance")

        # save fig
        plt.savefig(OUTPDF)
        plt.clf()

# ------------------------------------------- #

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter for guides made by intersecting FilteredGuides with peaks')
    parser.add_argument("-i", "--infile", dest="infile", type=str, help="File output by bedIntersect -wa -wb -a FilteredGuides.bed -b peaks.bed")
    parser.add_argument("-o", "--outprefix", dest="outprefix", type=str, help="Outputs files outprefix.getGuidePlots.pdf, outprefix.getGuides.txt, outprefix.getGuides.all.bedgraph, outprefix.getGuides.thin.bedgraph, outprefix.getGuides..thin.bed")
    parser.add_argument("-c", "--chrom", dest="chrom", type=str, help="chromosome of region to tile")
    parser.add_argument("-s", "--start", dest="start", type=int, help="start position")
    parser.add_argument("-e", "--end", dest="end", type=int, help="end position")
    parser.add_argument("-S", "--SSCMin", dest="SSCMIN", type=float, help="hard minimum of SSC score",default=-1.0)
    parser.add_argument("-O", "--OTMin", dest="OTMIN", type=int, help="hard minimum of Off Target Score in score column of bedfile",default=50)
    parser.add_argument("-T", "--PolyT", dest="POLYT", type=int, help="Max consecutive Ts",default=5)
    parser.add_argument("-V", "--PolyV", dest="POLYV", type=int, help="Max consecutive As, Gs, or Cs",default=7)    
    parser.add_argument("-E", "--EndT", dest="ENDT", type=int, help="Max consecutive Ts at end of guide",default=4)
    parser.add_argument("-G", "--MaxG", dest="MAXG", type=int, help="Max Total G", default=10)
    parser.add_argument("-g", "--MinGC", dest="MINGC", type=float, help="Min GC content", default=0)
    parser.add_argument("-x", "--MaxGC", dest="MAXGC", type=float, help="Max GC content", default=1)
    parser.add_argument("-D", "--MinStartDistance", dest="MINSTARTDIST", type=int, help="Minimum distance between starts of consecutive guides. If closer, take higher SSC", default=0)
    parser.add_argument("-P", "--Plot", dest="PLOT", action='store_true')
    parser.set_defaults(PLOT=False)

    args = parser.parse_args() 
    # (INFILE, LOCUS, SSCMIN, OTMIN, POLYT, POLYV, ENDTCount, MAXG, MINGC, MAXGC, MINSTARTDIST)
    LOCUS=(args.chrom, args.start, args.end)
    GetTileGuides(args.infile, args.outprefix, LOCUS, args.SSCMIN, args.OTMIN, args.POLYT, args.POLYV, args.ENDT, args.MAXG, args.MINGC, args.MAXGC, args.MINSTARTDIST, args.PLOT)


