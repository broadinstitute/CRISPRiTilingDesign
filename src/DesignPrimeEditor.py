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
    
    ## Provide either the variant on the command line:
    parser.add_argument('--chr', help="Chromosome")
    parser.add_argument('--start', type=int, help="Start position for edit")
    parser.add_argument('--end', type=int, help="End position for edit")
    parser.add_argument('--ref', type=str, help="Confirm the genomic sequence intended to edit")
    parser.add_argument('--alt', type=str, help="Desired edited sequence")

    ## Or provide a table of edits:
    parser.add_argument('--edits', help="Input file with columns: chr start end name ref alt region")

    ## Other inputs:
    parser.add_argument('--guides', required=required_args, help="Input guide file (e.g., .preDesign.bed file output by GetTileGuides.py) [cols: chr     start   end   strand  GuideSequenceWithPAM    guideSet]. Extra columns okay; will be ignored and output unchanged in final table")
    parser.add_argument('--fasta', required=required_args, default="/seq/lincRNA/data/hg19/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa", help="Input indexed FASTA file of genome")

    ## Output files and options
    parser.add_argument('--outfile', required=required_args, help="Filebase for output files")
    parser.add_argument('--splitOutputByRegion', action='store_true', default=False, help="Set to true to output separate design files for each 'region' in the input edits file")

    ## pegRNA design options:
    parser.add_argument('--minPbsLength', default=8, type=int, help="Minimum length of the primer binding sequence")
    parser.add_argument('--maxPbsLength', default=31, type=int, help="Maximum length of the primer binding sequence")
    parser.add_argument('--minRTPastEdit', default=8, type=int, help="Minimum length of the RT template")
    parser.add_argument('--maxRTPastEdit', default=20, type=int, help="Minimum length of the RT template")
    parser.add_argument('--maxRTTemplateLength', default=78, type=int, help="Maximum length of the RT template (78 = what was demonstrated with LoxP insertion)")
    parser.add_argument('--minPE3NickDistance', default=50, type=int, help="Minimum distance between the first nick and second PE3 nick")
    parser.add_argument('--maxPE3NickDistance', default=0, type=int, help="Maximum distance between the first nick and second PE3 nick. Default value (0) ignores PE3 nicking")
    parser.add_argument('--minPbsGcContent', default=30, type=float, help="Minimum GC content (%%) of the primer binding sequence")
    parser.add_argument('--minRTGcContent', default=30, type=float, help="Minimum GC content (%%) of the RT template")
    parser.add_argument('--maxUUcount', default=3, type=float, help="Maximum number of UU dinucleotides in the PBS + RT template (to avoid sequences likely to terminate Pol III Transcription)")

    ## Options for including sequences for pooled cloning
    parser.add_argument('--vectorDesigns', default="data/CloningDesigns.txt", help="Master file with gibson arm sequences for various plasmid designs")
    parser.add_argument('--vector', type=str, help="Name of vector to index into vectorDesigns")

    args = parser.parse_args()
    return(args)



def designPegRNAsForVariant(edit, guides, args):

    ## Expand bounds to account for max distance to PE3, RT template, etc.
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

    ## If this information is provided, filter guides to those corresponding to the edit set
    if ('region' in edit) and ('guideSet' in guides):
        guides = guides[guides['guideSet'] == edit['region']]

    ## For each potential editing gRNA:
    pegs = []
    for index, guide in guides.iterrows():
        curr = getAllPegRNAs(edit['chr'], seqStart, seqEnd, seq, guide['start'], guide['end'], guide['strand'], edit['start'], edit['end'], edit['alt'], args.minPbsLength, args.maxPbsLength, args.minRTPastEdit, args.maxRTPastEdit, args.maxRTTemplateLength)
        if (len(curr) > 0):
            ## Get nicking sites (PE3)
            if args.maxPE3NickDistance > 0:
                nickGuides = getNickingGuides(guides, edit['chr'], seqStart, seqEnd, seq, guide['start'], guide['end'], guide['strand'], args.minPE3NickDistance, args.maxPE3NickDistance)
                if (len(nickGuides) > 0):
                    combos = [ row1.append(row2) for i1,row1 in curr.iterrows() for i2,row2 in nickGuides.iterrows() ]
                    curr = pd.DataFrame(combos)
            pegs.append(curr)

    if (len(pegs) > 0):
        pegs = pd.concat(pegs)
        ## Format output
        pegs['variantName'] = edit['name']
        pegs['ref'] = edit['ref']
        pegs['alt'] = edit['alt']
        pegs['region'] = edit['region']
    else:
        pegs = pd.DataFrame()

    ## Return
    return pegs



def filterPegs(pegs, maxUUcount=4, minPbsGcContent=0, minRTGcContent=0):
    fp = pegs

    ## Filter based on GC content:
    fp = fp[fp["pbsGCpct"] >= minPbsGcContent] 
    fp = fp[fp["rtGCpct"] >= minRTGcContent]

    ## Filter based on UU count:
    fp = fp[fp["uuCount"] <= maxUUcount]

    ## ? Remove any pegRNAs that have 7-bp windows with 5 or more T/Us (which would reduce Pol III transcription)

    ## To do: Add filter to choose only the closest spacer? 

    return fp


def addPoolCloningOligos(df, vectorDesignFile, vector):  #, MapLength=21):
    vectorDesigns = read_table(vectorDesignFile)
    if not vector in vectorDesigns['CloningDesign'].values:
        raise ValueError("Chosen vector design '"+vector+"' is not in the Vector Design file.")
    leftGA = vectorDesigns.loc[vectorDesigns['CloningDesign'] == vector].LeftGibsonArm.item()
    rightGA = vectorDesigns.loc[vectorDesigns['CloningDesign'] == vector].RightGibsonArm.item()
    scaffold = vectorDesigns.loc[vectorDesigns['CloningDesign'] == vector].pegRNAScaffold.item()
    df["GuideSequence"]=df['spacer'].apply(lambda x: prependGifNecessary(x))
    df["GuideSequenceMinusG"]=df["GuideSequence"].apply(lambda x: trimG(x))
    df["LeftGA"]=leftGA
    df["RightGA"]=rightGA
    df["CoreOligo"]=[leftGA+getPegRNAFromPandas(row).getPegRNASequence(scaffold=scaffold)+leftGA for idx,row in df.iterrows()]
    #df["MappingSequence"]=(df["GuideSequenceMinusG"]+df["RightGA"]).apply(lambda x: x[0:MapLength])
    return df
    

def pegsToBED(pegTable):
    bed = []
    for idx,row in pegTable.iterrows():
        bed.append(getPegRNAFromPandas(row).toBED())
    bed = pd.concat(bed, axis=1).transpose()
    bed['name'] = pegTable['pegID'].values
    return bed


def getVariantCoverageBedgraph(pegTable):
    bedgraph = pegTable.groupby(['chr','variantStart','variantEnd']).size().reset_index()
    return bedgraph


def writePegRNAs(results, outfile, splitOutputByRegion):
    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    results.to_csv(outfile+".full.tsv", sep='\t', header=True, index=False)
    resultsBED = pegsToBED(results)
    resultsBED.to_csv(outfile + ".full.bed", sep='\t', header=False, index=False)
    resultsBedgraph = getVariantCoverageBedgraph(results)
    resultsBedgraph.to_csv(outfile + ".full.bedgraph", sep='\t', header=False, index=False)

    if splitOutputByRegion:
        for name, resultsGroup in results.groupby('region'):
            resultsGroup.to_csv(outfile+"."+name+".tsv", sep='\t', header=True, index=False)
            pegsToBED(resultsGroup).to_csv(outfile+"."+name+".bed", sep='\t', header=False, index=False)
            getVariantCoverageBedgraph(resultsGroup).to_csv(outfile+"."+name+".bedgraph", sep='\t', header=False, index=False)


def main(args):
    guides = read_table(args.guides)

    required_cols = ['chr','start','end','strand','GuideSequenceWithPAM','guideSet']
    for col in required_cols:
        if not len(guides[col])>0:
            raise ValueError("Guide file must contain column: " + col)

    ## Parse the variants/edits to create
    if args.edits is None:
        edits = pd.DataFrame({ 
            'chr' : [args.chr],
            'start' : [args.start],
            'end' : [args.end],
            'ref' : [args.ref],
            'alt' : [args.alt],
            'name' : ["Custom"],
            'region' : ["CustomRegion"]
            })
    else:
        edits = read_table(args.edits)

    ## Design pegRNAs
    results = []
    for index, edit in edits.iterrows():
        curr = designPegRNAsForVariant(edit, guides, args)
        if (len(curr) > 0):
            results.append(curr)
    
    ## Format and filter pegRNA list
    results = pd.concat(results)
    results = results.fillna(0)
    results = filterPegs(results, args.minPbsGcContent, args.minRTGcContent)
    results['pegID'] = ["peg"+str(n)+"-"+str(region)+"-"+vname for region,n,vname in zip(results['region'],range(1,len(results)+1),results['variantName'])]
    results = results.astype(str)

    ## Add pool cloning sequences, if requested
    if args.vector:
        results = addPoolCloningOligos(results, args.vectorDesigns, args.vector)

    ## Write results
    writePegRNAs(results, args.outfile, args.splitOutputByRegion)


if __name__ == '__main__':
    args = parseargs()
    main(args)
