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



def writeBed(df, OUTFILE, score=None, negative=False):
    a=df.copy()
    a = a[a.chr.notnull()]
    a["start"]=np.vectorize(str)(np.vectorize(int)(a["start"]))
    a["end"]=np.vectorize(str)(np.vectorize(int)(a["end"]))
    cols=["chr", "start", "end"]
    if score:
        cols=["chr", "start", "end", score]
    if negative:
        a[score]=-a[score]
    a[cols].to_csv(OUTFILE, sep='\t', header=False, index=False)


def prependGifNecessary(thisGuide):
    '''Checks if the given guide starts with a G.
    If not, prepends one.'''
    if thisGuide[0]!="G":
        return "G"+thisGuide
    else:
        return thisGuide


def GCContent(x):
    x=x.upper()
    return (x.count("G")+x.count("C"))/len(x)


def trimPAM(seq, PAM="NGG"):
    if not PAM == "NGG":
        raise ValueError("Only NGG PAM supported at this time")
    assert seq[-2:]=="GG"
    return seq[0:-3]


def trimG(seq):
    seq=seq.upper()
    if seq[0]=="G":
        return seq[1:]
    else:
        return seq


def smartInt(x):
    try:
        return int(x)
    except ValueError:
        return x


def getEvenlySpacedIndices(arr, numElems):
    if numElems >= len(arr):
        return [i for i in range(len(arr))]
    else:
        return np.round(np.linspace(0, len(arr) - 1, numElems)).astype(int)

def getEvenlySpaced(arr, numElems):
    return arr[getEvenlySpacedIndices(arr, numElems)]

