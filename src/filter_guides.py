import sys
import pandas as pd
import argparse
import numpy as np
from collections import Counter

def parseargs():
    parser = argparse.ArgumentParser(description='filter guides by polyT and spacing')
    parser.add_argument('--polyTmax', help="guides with this or longer polyT sites will be removed", default=4, type=int)
    parser.add_argument('--minStartDistance', help="minimum space between guide starts", default=5, type=int)
    parser.add_argument('--minCFD', help="minimum CFD", default=0.2, type=float)
    parser.add_argument('--startWithG', help="require guide sequences to start with G", action='store_true')
    parser.add_argument('--minScore', help="minimum Score", default=50, type=float)
    parser.add_argument('--infile', help="input bed file", required=True)
    parser.add_argument('--outfile', help="output bed file", required=True)
    args = parser.parse_args()
    return(args)

def filter_out(df, mask, reason):
    mask = mask & (df['filtered out reason'].str.len() == 0)
    print("removing {} of {} for {}".format(mask.sum(), (df['filtered out reason'].str.len() == 0).sum(), reason))
    df.loc[mask, 'filtered out reason'] = reason


if __name__ == '__main__':
    args = parseargs()

    guides = pd.read_csv(args.infile, names="chr,start,end,gRNA,cutting efficiency score,cutting specificity score,strand,offtargets sum,offtargets summary,gRNA group,gRNA label".split(','))
    guides = guides.dropna(axis=0)
    guides = guides[guides.chr != 'chromosome']
    guides['start'] = guides.start.astype(int)
    guides = guides[guides["cutting efficiency score"] != '*']  # ???
    guides["cutting efficiency score"] = guides["cutting efficiency score"].astype(int)
    guides["cutting specificity score"] = guides["cutting specificity score"].astype(float)
    guides['filtered out reason'] = ""

    print(guides.columns)
    # optionally filter out guides that don't start with G
    if args.startWithG:
        filter_out(guides, guides.gRNA.str.startswith('G'), 'does not start with G')

    # filter out any Ns
    filter_out(guides, guides.gRNA.str.contains('N'), 'guide sequence contains N')

    # filter out U repeats
    filter_out(guides,  guides.gRNA.str.contains('T' * args.polyTmax), "too many consecutive Ts")

    # filter out anything with a 5-mer repeat
    filter_out(guides,  guides.gRNA.str.contains('A' * 5), "too many consecutive As")
    filter_out(guides,  guides.gRNA.str.contains('C' * 5), "too many consecutive Cs")
    filter_out(guides,  guides.gRNA.str.contains('G' * 5), "too many consecutive Gs")
    filter_out(guides,  guides.gRNA.str.contains('T' * 5), "too many consecutive Ts")

    # UContentPredicate - filter out more than 1 T in the last 4 bases, or more than 40% T overall
    filter_out(guides,  guides.gRNA.str[-4:].str.count('T') > 1, "more than 2 Ts in final 4 bases")
    guidelens = guides.gRNA.str.len()
    filter_out(guides,  guides.gRNA.str.count('T') > 0.4 * guidelens, "more than 40% T")

    # filter by GC content
    GC_content = guides.gRNA.str.count('C') + guides.gRNA.str.count('G')
    filter_out(guides, (GC_content < 0.2 * guidelens), "less than 20% GC content")
    filter_out(guides, (GC_content > 0.9 * guidelens), "more than 90% GC content")

    # filter by CFD and score
    filter_out(guides, guides["cutting efficiency score"] < args.minScore, "score below {}".format(args.minScore))
    filter_out(guides, guides["cutting specificity score"] < args.minCFD, "CFD below {}".format(args.minCFD))

    # filter low complexity - long stretches of repeated sequence
    def low_complexity(seq):
        # for a given k-mer size, how long a repeat do we require to be considered low-complexity?
        subseq_lens = {2: 10,
                       3: 12,
                       4: 12,
                       5: 18,
                       6: 18}

        for kmer, subseq_len in subseq_lens.items():
            # check each substring of subseq_len
            for offset in range(len(seq) - subseq_len + 1):
                subseq = seq[offset:offset + subseq_len]
                # check if it's repetitive - (uses transitivity of ==)
                if subseq[:-kmer] == subseq[kmer:]:
                    return True
        return False
    filter_out(guides, guides.gRNA.apply(low_complexity), "low complexity")


    # keep guides spaced out
    def remove_too_close(starts):
        starts = starts.values
        dists = abs(starts.reshape((-1, 1)) - starts.reshape((1, -1)))
        dists = np.maximum(dists, np.tri(dists.shape[0]).T * (1 + args.minStartDistance))
        bad = (dists < args.minStartDistance).sum(axis=1) > 0
        return ~ bad
    spaced_out = guides[guides["filtered out reason"].str.len() == 0].groupby("gRNA group").start.transform(remove_too_close)
    too_close = ~ guides.index.isin(spaced_out.index)
    filter_out(guides, too_close, "too close to another guide")

    guides = guides.rename(columns={"gRNA group": "guideSet"})
    guides.to_csv(args.outfile, sep="\t", index=None)
