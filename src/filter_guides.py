import sys
import pandas as pd
import argparse
import numpy as np

def parseargs():
    parser = argparse.ArgumentParser(description='filter guides by polyT and spacing')
    parser.add_argument('--polyTmax', help="guides with this or longer polyT sites will be removed", default=4)
    parser.add_argument('--minStartDist', help="minimum space between guide starts", default=5)
    parser.add_argument('--infile', help="input bed file", required=True)
    parser.add_argument('--outfile', help="output bed file", required=True)
    args = parser.parse_args()
    return(args)


if __name__ == '__main__':
    args = parseargs()

    guides = pd.read_csv(args.infile, names="chr,start,end,gRNA,cutting efficiency score,cutting specificity score,strand,offtargets sum,offtargets summary,gRNA group,gRNA label".split(','))
    guides = guides.dropna(axis=0)
    guides = guides[guides.chr != 'chromosome']
    guides['start'] = guides.start.astype(int)

    guides = guides[~ guides.gRNA.str.contains('T' * args.polyTmax)]


    def remove_too_close(grp):
        starts = grp.start.values
        dists = abs(starts.reshape((-1, 1)) - starts.reshape((1, -1)))
        dists = np.maximum(dists, np.tri(dists.shape[0]).T * (1 + args.minStartDist))
        bad = (dists < args.minStartDist).sum(axis=1) > 0
        return grp.loc[~bad, :]

    guides = guides.groupby("gRNA group").apply(remove_too_close)
    guides.to_csv(args.outfile, sep="\t", index=None)
