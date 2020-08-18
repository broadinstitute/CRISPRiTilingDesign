import sys
import pandas as pd
import argparse

def parseargs():
    parser = argparse.ArgumentParser(description='make name column of bed unique')
    parser.add_argument('suffix', help="extra suffix to include (e.g. 'TSS' will insert -TSS1, -TSS2...)")
    parser.add_argument('infile', help="input bed file")
    parser.add_argument('outfile', help="output bed file")
    args = parser.parse_args()
    return(args)


if __name__ == '__main__':
    args = parseargs()

    def add_suffix(grp):
        if len(grp) == 1:
            return [grp.name]
        else:
            return ['{}-{}{}'.format(grp.name, args.suffix, n) for n in range(1, len(grp) + 1)]

    genes = pd.read_table(args.infile, header=None)
    cols = list(genes.columns)
    cols[3] = 'symbol'
    genes.columns = cols
    tmp = genes.groupby("symbol").apply(add_suffix).reset_index(drop=True)
    # flatten list of lists
    genes['symbol'] = [item for sublist in tmp for item in sublist]
    genes.to_csv(args.outfile, header=None, index=None, sep="\t")
