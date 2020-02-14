###############################################################################
## Library code for designing CRISPRi screens
## Jesse Engreitz
## November 10, 2019
## Based on Charlie's gRNA design
## Tested with:  "use .python-3.5.1; source /seq/lincRNA/Ben/VENV_MIP/bin/activate"

## TODO: Factor out helper code


from PrimeEditor import *
import pybedtools
from pybedtools import BedTool
import pandas as pd
import numpy as np
import argparse
from pandas.io.parsers import read_table


def parseargs(required_args=True):
    class formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass

    epilog = ("")
    parser = argparse.ArgumentParser(description='Output list of potential prime editor gRNAs for desired variants',
                                     epilog=epilog,
                                     formatter_class=formatter)
    readable = argparse.FileType('r')
    
    parser.add_argument('--chr', help="Chromosome")
    parser.add_argument('--start', type=int, help="Start position for edit")
    parser.add_argument('--end', type=int, help="End position for edit")
    parser.add_argument('--ref', type=str, help="Confirm the genomic sequence intended to edit")
    parser.add_argument('--alt', type=str, help="Desired edited sequence")
    parser.add_argument('--edits', help="Input file with columns: chr start end name ref alt")
    parser.add_argument('--outfile', required=required_args, help="Output file")
    parser.add_argument('--guides', required=required_args, help="Guide file (e.g., .preDesign.bed file output by GetTileGuides.py) [cols: chr     start   end   strand  GuideSequenceWithPAM    guideSet]. Extra columns okay; will be ignored and output unchanged in final table")
    parser.add_argument('--fasta', required=required_args, default="/seq/lincRNA/data/hg19/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa", help="Indexed FASTA file")

    parser.add_argument('--minPbsLength', default=8, type=int, help="Minimum length of the primer binding sequence")
    parser.add_argument('--maxPbsLength', default=31, type=int, help="Maximum length of the primer binding sequence")
    parser.add_argument('--minRTPastEdit', default=8, type=int, help="Minimum length of the RT template")
    parser.add_argument('--maxRTPastEdit', default=20, type=int, help="Minimum length of the RT template")
    parser.add_argument('--maxRTTemplateLength', default=78, type=int, help="Maximum length of the RT template (78 = what was demonstrated with LoxP insertion)")
    parser.add_argument('--minPE3NickDistance', default=50, type=int, help="Minimum distance between the first nick and second PE3 nick")
    parser.add_argument('--maxPE3NickDistance', default=0, type=int, help="Maximum distance between the first nick and second PE3 nick")

    args = parser.parse_args()
    return(args)


def designPegRNAsForVariant(edit, guides, args):

    ## Expand bounds to acccount for max distance to PE3, RT template, etc.
    EXTEND = 200
    seqStart = edit['start'] - EXTEND
    seqEnd = edit['end'] + EXTEND

    ## Range intersect with guides
    currGuides = guides.loc[(guides['chr'] == edit['chr']) & (guides['start'] >= edit['start']) & (guides['end'] <= edit['end'])]

    ## Assert that base in genomic sequence = ref, as a sanity check
    ## TODO: Need to rewrite the guide designer so it can create guides for editing alt to ref
    fasta = pybedtools.example_filename(args.fasta)
    toreplace = BedTool.seq((edit['chr'], edit['start'], edit['end']), fasta)
    if not toreplace.upper() == edit['ref'].upper():
        raise ValueError("Reference sequence does not match proposed ref sequence: " + toreplace + " vs " + edit['ref'])

    ## Get genomic sequence
    seq = BedTool.seq((edit['chr'], seqStart, seqEnd), fasta)

    ## For each potential editing gRNA:
    pegs = pd.DataFrame()
    for index, guide in guides.iterrows():
        curr = getAllPegRNAs(edit['chr'], seqStart, seqEnd, seq, guide['start'], guide['end'], guide['strand'], edit['start'], edit['end'], edit['alt'], args.minPbsLength, args.maxPbsLength, args.minRTPastEdit, args.maxRTPastEdit, args.maxRTTemplateLength)
        if (len(curr) > 0):

            ## Get nicking sites (PE3)
            if args.maxPE3NickDistance > 0:
                nickGuides = getNickingGuides(guides, edit['chr'], seqStart, seqEnd, seq, guide['start'], guide['end'], guide['strand'], args.minPE3NickDistance, args.maxPE3NickDistance)
                if (len(nickGuides) > 0):
                    combos = [ row1.append(row2) for i1,row1 in curr.iterrows() for i2,row2 in nickGuides.iterrows() ]
                    curr = pd.DataFrame(combos)

            pegs = pegs.append(curr)[curr.columns.tolist()]

    ## Format output
    pegs['variantName'] = edit['name']
    pegs['ref'] = edit['ref']
    pegs['alt'] = edit['alt']

    ## Return
    return pegs



def main(args):
    guides = read_table(args.guides)

    required_cols = ['chr','start','end','strand','GuideSequenceWithPAM','guideSet']
    for col in required_cols:
        if not len(guides[col])>0:
            raise ValueError("Guide file must contain column: " + col)

    if args.edits is None:
        edits = pd.DataFrame({ 
            'chr' : [args.chr],
            'start' : [args.start],
            'end' : [args.end],
            'ref' : [args.ref],
            'alt' : [args.alt],
            'name' : ["Custom"]
            })
    else:
        edits = read_table(args.edits)

    results = pd.DataFrame()
    for index, edit in edits.iterrows():
        curr = designPegRNAsForVariant(edit, guides, args)
        if (len(curr) > 0):
            results = results.append(curr)[curr.columns.tolist()]
    
    results.to_csv(args.outfile, sep='\t', header=True, index=False)



if __name__ == '__main__':
    args = parseargs()
    main(args)
