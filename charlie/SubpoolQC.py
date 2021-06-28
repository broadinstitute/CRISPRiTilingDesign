#! /usr/bin/env python

# USAGE: python ./SubpoolQC.py COUNTFILE DESIGNFILE OUTFILE SUBPOOLCOLUMN IDCOLUMN

from __future__ import division
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import sys
import numpy as np
import pandas as pd
from pandas.io.parsers import read_table
from scipy import stats
import argparse


# ------------------------------------------- #
def skewRatio(df, col="count"):
    # skew ratio of top 10% to bottom 10% of guide counts
    [top_10, bottom_10] = [np.percentile(df[col], 90), np.percentile(df[col], 10)]
    if top_10 >0:
		return top_10/bottom_10
    else:
		return "Insufficient matches for skew ratio"


def drawDensity(vec, label=None, fill=False):
    x=np.arange(min(vec), max(vec), (max(vec)-min(vec))/1000)
    density=stats.kde.gaussian_kde(vec)
    if label:
        plt.plot(x, density(x), label=label)
    else:
        plt.plot(x, density(x))
    if fill:
        plt.fill_between(x, 0, density(x))

# ------------------------------------------- #
## Read down-sampling functions
def get_regression_stats(X,Y): #This function performs linear regression and returns slope,r-squared,and p value
    slope, intercept, r_value, p_value, std_err = stats.linregress(X,Y)
    dictionary = {'slope':slope,'R-squared':r_value**2,'P-value': p_value, 'R':r_value}
    return dictionary

def sampleReads(v, FractionReads):
    '''
    Generates a new vector of with downsampled reads for each item.
    Does not exactly hit number of reads
    '''
    assert FractionReads<1
    assert FractionReads>0
    
    sample=[0]*len(v)
    for row in range(len(v)):
        sample[row]=sum(np.random.rand(v[row])<FractionReads)
    return sample

def downsamplePearson(v, FractionReads, log):
    '''Finds pearson r of original vector and downsampled one. (optional log10 of both)'''
    sample=sampleReads(v, FractionReads)
    if log:
        return get_regression_stats(np.log10(v+1), np.log10(np.array(sample)+1))["R"]
    return get_regression_stats(v, sample)["R"]

def downsample(v, fracList, log=False):
    '''downsample several fractions and return df with fraction and pearsonr'''
    
    LOD=[]
    for frac in fracList:
        data=downsamplePearson(v, frac, log)
        LOD.append({"Fraction":frac, "PearsonR":data})
    return pd.DataFrame(LOD)

# ------------------------------------------- #
# QC Function
# ------------------------------------------- #

def subpoolQC(COUNTFILE, DESIGNFILE, SUBPOOLCOL, OUTPREFIX, MINCOUNTS, IDCOL, LOG):
	OUTFILE=OUTPREFIX+".txt"
	OUTPDF=OUTPREFIX+".pdf"

	count=read_table(COUNTFILE, names=["count", IDCOL])
    	design=read_table(DESIGNFILE)
    
	# --------------------------------------------- #
	# Plot to see if read coverage is saturating
	plt.subplot(1,2,1)
	v=count["count"]
        data=downsample(v, [0.01,0.05,0.25,0.5,0.9], log=LOG)
        plt.scatter(data["Fraction"], data["PearsonR"])
    	plt.xlabel("Fraction of Reads"+" log10: "+str(LOG))
	plt.ylabel("Pearson Correlation With All Reads")

	# --------------------------------------------- #	

    	data=pd.merge(design, count[count["count"]>=MINCOUNTS], on=IDCOL)
    	dataAll=pd.merge(design, count, on=IDCOL)
    
    	outfile=open(OUTFILE, 'w')
	
	# write the number of guides in each subpool
    	for i in zip(data[SUBPOOLCOL].value_counts().index, data[SUBPOOLCOL].value_counts()):
        	print >> outfile, i
    	print >>outfile, "---------------------------------"
    
    	print SUBPOOLCOL
    	print data[SUBPOOLCOL].drop_duplicates()
    	
	plt.subplot(1,2,2)
	for row in range(min(len(data[SUBPOOLCOL].drop_duplicates()), 7)):
        	tGene=dataAll[SUBPOOLCOL].value_counts().index[row]
		print tGene
        	ontarget=dataAll[dataAll[SUBPOOLCOL]==tGene]
        	drawDensity(np.log10(ontarget["count"]), label=tGene)
    
        	designed=len(design[design[SUBPOOLCOL]==tGene])
        	made=len(data[data[SUBPOOLCOL]==tGene])
        	madeAll=len(ontarget)
        
        	undetected=designed-madeAll
        	lowReads=madeAll-made
        	detected=made
        
        	outfile.write(str(tGene)+"\n")
        	outfile.write("Guides in subpool: "+str(designed)+"\n")
        	outfile.write("Guides detected: "+str(detected)+" "+ str(detected/designed)+"\n")
        	outfile.write("Guides low reads: "+str(lowReads)+" "+ str(lowReads/designed)+"\n")
        	outfile.write("Guides undetected "+str(undetected)+" "+ str(undetected/designed)+"\n")
        	outfile.write("Skew Ratio: "+str(skewRatio(ontarget))+"\n")
        	outfile.write("---------------------------------\n\n")
    
    	outfile.close()
    
    	plt.xlabel("log10 reads")
    	plt.legend(loc="best")
	plt.tight_layout()
	plt.savefig(OUTPDF)
    	plt.clf()

# ------------------------------------------- #

if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='Plot and count distributions of a subpool cloning')
        parser.add_argument("-c", "--countFile", dest="countFile", type=str,  help="File with unnamed columns for counts per guide and OligoID")
        parser.add_argument("-d", "--designFile", dest="designFile", type=str, help="File with named columns with at least OligoID and subpoolColumn")
        parser.add_argument("-s", "--subpoolColumn", dest="subpoolColumn", type=str, help="The name of the column in designFile designating subpool", default="target")
        parser.add_argument("-o", "--out", dest="out", type=str, help="Where to save output: [out].pdf, [out].txt", default="QC.out")
        parser.add_argument("-m", "--minCounts", dest="minCounts" , type=int, help="Minimum counts per guide for listing subpools represented", default=5)
	parser.add_argument("-i", "--idColumn", dest="idColumn" , type=str, help="The name of the column in designFile designating the guide ID used in the bowtie index", default="OligoID")
	parser.add_argument("-L", "--Linear", dest="LOG", action="store_false", help="Compute Pearson R using linear scale",default=True)

        args = parser.parse_args()

	subpoolQC(args.countFile, args.designFile, args.subpoolColumn, args.out, args.minCounts, args.idColumn, args.LOG)
