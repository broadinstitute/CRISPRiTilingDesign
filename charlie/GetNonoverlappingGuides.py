#! /usr/bin/env python

# USAGE: python ./GetTileGuides.py

# ------------------------------------------- #
import pandas as pd
from pandas.io.parsers import read_table
import argparse
import numpy as np
import sys
from scipy import stats

import warnings
warnings.filterwarnings("ignore")
# ------------------------------------------- #
# helpers
def limitToLocus(df, (CHR, START, END)):
    a=df[df["chr"]==CHR]
    a=a[a["start"]>=int(START)]
    a=a[a["end"]<=int(END)]
    return a

def rankswitch(RANK):
	#returns True if the best guides have score low, False if the best guides have score high
	if RANK.upper()=="LOW":
		return True
	elif RANK.upper()=="HIGH":
		return False
	else:
		raise ValueError('-r/--Rank accepts only Low or High as input')

def checkBEDOverlap(LINE, BED):

    for row in range(len(BED)):
        
        if checkLINEOverlap(LINE, BED.irow(row)):
            return True
    return False
    
def checkLINEOverlap(LINE1, LINE2):
    if LINE1["chr"]!=LINE2["chr"]:
        return False
    elif(LINE1["end"]>LINE2["start"]) & (LINE1["start"]<LINE2["end"]):
        return True
    return False

def writeOutput(df, OUTFILE):
    a=df.copy()
    a["start"]=np.vectorize(str)(np.vectorize(int)(a["start"]))
    a["end"]=np.vectorize(str)(np.vectorize(int)(a["end"]))
    a.to_csv(OUTFILE, sep='\t', header=True, index=False)

# ------------------------------------------- #
# Main method

def grabGuides(INFILE, LOCUS, SCORECOL, RANK, OUTFILE, N, NOOVERLAP):
	print str(N)
	print str(NOOVERLAP)
	
	df=read_table(INFILE)
	
	# limit to locus!
	df=limitToLocus(df, LOCUS)
	
	df=df.sort(SCORECOL, ascending=rankswitch(RANK))
	
	print LOCUS

	if NOOVERLAP:
		# automatically take first, then the rest
   		winners=pd.DataFrame(df.irow([0]))
		for row in range(1,len(df)):
		
			# check if overlapping and not done yet
			if (len(winners)<N) & (checkBEDOverlap(df.irow(row), winners)==False):
				winners=pd.concat([winners, df.irow([row])])
	else:
		winners=df.copy()
	writeOutput(winners.irow(range(0,N)), OUTFILE)

# ------------------------------------------- #

if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='Outputs rows from INFILE with best N non-overlaping guides')
        parser.add_argument("-i", "--infile", dest="INFILE", type=str, help="bed-like file with header of guides with a score column")
	parser.add_argument("-c", "--chrom", dest="chrom", type=str, help="chromosome of region to tile")
	parser.add_argument("-s", "--start", dest="start", type=int, help="start position")
	parser.add_argument("-e", "--end", dest="end", type=int, help="end position")
	parser.add_argument("-o", "--outfile", dest="OUTFILE", type=str, help="OUTFILE")
	parser.add_argument("-S", "--scoreCol", dest="SCORECOL", type=str, help="score column")
	parser.add_argument("-r", "--rank", dest="rank", type=str, default="Low" ,help="Best guides have 'High' or 'Low' scores")
	parser.add_argument("-N", "--NGuides", dest="N", type=int, help="hard minimum of Off Target Score in score column of bedfile",default=50)
	parser.add_argument("-O", "--AllowOverlap", dest="NoOverlap", action="store_false", help="Column with OffTargetScore",default=True)
	#parser.add_argument("-C", "--OTCOL", dest="OTCOL", type=int, help="Column with OffTargetScore",default)

	args = parser.parse_args() 
	LOCUS=(args.chrom, args.start, args.end)
	
	grabGuides(args.INFILE, LOCUS, args.SCORECOL, args.rank, args.OUTFILE, args.N, args.NoOverlap)
