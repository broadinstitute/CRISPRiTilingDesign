###############################################################################
## Library code for designing CRISPRi screens
## Jesse Engreitz
## November 14, 2019
## Design Prime Editor gRNAs


import Bio
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from Bio import SeqUtils
import pandas as pd
import numpy as np
from collections import OrderedDict
from JuicerCRISPRiDesign import *

NICK_POSITION=-3
SCAFFOLD     ="GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
OPTI_SCAFFOLD="GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"

class PrimeEditor:
    ## TODO:  Change to extend pd.DataFrame or pd.Series

    def __init__(
        self, 
        spacer, 
        primerBindingSite, 
        rtTemplate, 
        guideRegion=None, 
        variantRegion=None, 
        pbsRegion=None, 
        rtTemplateRegion=None,
        name=""):
        '''
        spacer, primerBindingSite, rtTemplate are strings representing sequence in pegRNA.
        Regions are in genomic coordinates
        '''
        self.guideRegion = guideRegion
        self.variantRegion = variantRegion
        self.pbsRegion = pbsRegion
        self.rtTemplateRegion = rtTemplateRegion
        self.spacer = str(spacer)
        self.primerBindingSite = str(primerBindingSite)
        self.rtTemplate = str(rtTemplate)
        #self.scaffold = str(scaffold)
        self.spacerPlusG = prependGifNecessary(self.spacer)
        self.name = name

    def __str__(self):
        return "PrimeEditor: "+ str(self.getPegRNASequence())
        
    def getPegRNASequence(self, prependGifNecessary=True, scaffold=SCAFFOLD):
        if prependGifNecessary is True:
            return self.spacerPlusG + scaffold + self.getExtensionSequence()
        else:
            return self.spacer + scaffold + self.getExtensionSequence()
    
    def getExtensionSequence(self):
        return self.rtTemplate + self.primerBindingSite
    
    '''
    def getLiuLabCloningOligos(self):
        return pd.DataFrame({
            "SpacerTop" : ["CACC" + self.spacer + "GTTTT"],
            "SpacerBot" : [str("CTCTAAAAC" + Seq(self.spacer).reverse_complement())],
            "ExtensionTop" : ["GTGC" + self.getExtensionSequence()],
            "ExtensionBot" : [str("AAAA" + Seq(self.getExtensionSequence()).reverse_complement())]
        })
    '''

    def getRtLengthPastEdit(self):
        if self.guideRegion is None:
            raise ValueError("PrimeEditor was not initialized with a guide")
        return abs(self.rtTemplateRegion.three_prime().distanceFrom(self.variantRegion.three_prime()))

    def editOverlapsPam(self):
        pamGG = self.guideRegion.slop(-(self.guideRegion.end-self.guideRegion.start+1), 3, considerStrand=True)
        return ((pamGG.end > self.variantRegion.start) & (pamGG.start < self.variantRegion.end))

    def toPandas(self):
        return pd.Series( OrderedDict((
            ("chr", self.guideRegion.chromosome),
            ("guideStart", self.guideRegion.start),
            ("guideEnd", self.guideRegion.end),
            ("guideStrand", self.guideRegion.strand),
            ("variantStart", self.variantRegion.start),
            ("variantEnd" , self.variantRegion.end),
            ("pbsStart" , self.pbsRegion.start),
            ("pbsEnd" , self.pbsRegion.end),
            ("rtTemplateStart" , self.rtTemplateRegion.start),
            ("rtTemplateEnd" , self.rtTemplateRegion.end),
            ("rtTemplateLength" , len(self.rtTemplate)),
            ("rtLengthPastEdit" , self.getRtLengthPastEdit()),
            ("warningFlapEndsInG" , self.rtTemplate[0] == "C"),
            ("editPositionRelativeToNick0Based", self.getPositionRelativeToNick()),
            ("editOverlapsPam" , self.editOverlapsPam()),
            ("pbsLength" , len(self.primerBindingSite)),
            ("pbsGCpct" , SeqUtils.GC(self.primerBindingSite)),
            ("spacer" , self.spacer),
            ("primerBindingSite" , self.primerBindingSite),
            ("rtTemplate" , self.rtTemplate),
            ("rtGCpct" , SeqUtils.GC(self.rtTemplate)),
            ("uuCount", self.getExtensionSequence().count("TT")),
            ("pegRNASequence" , self.getPegRNASequence()),
            ("SpacerOligoTop", "CACC" + self.spacerPlusG + "GTTTT"),
            ("SpacerOligoBot", str("CTCTAAAAC" + Seq(self.spacerPlusG).reverse_complement())),
            ("ExtensionOligoTop", "GTGC" + self.getExtensionSequence()),
            ("ExtensionOligoBot", str("AAAA" + Seq(self.getExtensionSequence()).reverse_complement())),
            ("GibsonOligoTop", "aaaggacgaaacacc" + self.spacerPlusG + "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAAC"),
            ("GibsonOligoBot", str("gcggcccaagcttaaaaaaa" + Seq(self.getExtensionSequence()).reverse_complement() + "GCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCT")),
            ("sgOpti.pegRNASequence" , self.getPegRNASequence(scaffold=OPTI_SCAFFOLD)),
            ("sgOpti.SpacerOligoTop", "CACC" + self.spacerPlusG + "GTTTA"),
            ("sgOpti.SpacerOligoBot", str("CTCTTAAAC" + Seq(self.spacerPlusG).reverse_complement())),
            ("sgOpti.ExtensionOligoTop", "GTGC" + self.getExtensionSequence()),
            ("sgOpti.ExtensionOligoBot", str("AAAA" + Seq(self.getExtensionSequence()).reverse_complement())),
            ("sgOpti.GibsonOligoTop", "aaaggacgaaacacc" + self.spacerPlusG + "GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAAC"),
            ("sgOpti.GibsonOligoBot", str("gcggcccaagcttaaaaaaa" + Seq(self.getExtensionSequence()).reverse_complement() + "GCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCT"))
            )))

    def getPositionRelativeToNick(self):
       return PrimeEditor.GetPositionRelativeToNick(self.variantRegion.get5pEnd(), self.guideRegion)

    def GetNickCoordinate(guideRegion):
        if guideRegion.strand == "-":
            return guideRegion.start - NICK_POSITION
        else:
            return guideRegion.end + NICK_POSITION

    def GetPositionRelativeToNick(position, guideRegion):
        ''' Note: 0-based '''
        distance = position - PrimeEditor.GetNickCoordinate(guideRegion)
        if guideRegion.strand == "-":
            return -1*distance
        else:
            return distance

    def toBED(self):
        if (self.guideRegion.strand == "+"):
            blockSizes = [self.guideRegion.length()-3, self.rtTemplateRegion.end-self.variantRegion.start]
            blockStarts = [0,self.variantRegion.start-self.guideRegion.start]
        else:
            blockSizes = [self.variantRegion.end-self.rtTemplateRegion.start, self.guideRegion.length()-3]
            blockStarts = [0,self.guideRegion.start+3-self.rtTemplateRegion.start]
        return(pd.Series( OrderedDict((
            ("chr", self.guideRegion.chromosome),
            ("start", min([self.guideRegion.start,self.rtTemplateRegion.start])),
            ("end", max([self.guideRegion.end,self.rtTemplateRegion.end])),
            ("name", self.name),
            ("score", 0),
            ("strand", self.guideRegion.strand),
            ("thickStart", self.variantRegion.start),
            ("thickEnd", self.variantRegion.end),
            ("itemRgb", "0,0,0"),
            ("blockCount", 2),
            ("blockSizes", ",".join([str(i) for i in blockSizes])),
            ("blockStarts", ",".join([str(i) for i in blockStarts]))
            ))))

    def getBounds(self):
        if self.guideRegion.strand == "-":
            return GenomicRange(self.guideRegion.chromosome, self.rtTemplateRegion.start, self.guideRegion.end, self.guideRegion.strand)
        else:
            return GenomicRange(self.guideRegion.chromosome, self.guideRegion.start, self.rtTemplateRegion.end, self.guideRegion.strand)



class GenomicRange:
    ## TODO:  Rewrite to extend pyranges

    def __init__(self, chromosome, start, end, strand="*", name=None):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.strand = strand
        self.name = name

        if not (strand in set(['+','-','*',''])):
            raise ValueError("Genomic Range: strand must be + - * or blank")

        if end < start:
            self.start = end
            self.end = start
            #raise ValueError("GenomicRange: end must be greater than start")

    def __str__(self):
        return self.chromosome + ":" + str(self.start) + "-" + str(self.end)

    def shift(self, distance, considerStrand=False):
        return self.slop(-distance, distance, considerStrand)

    def length(self):
        return self.end - self.start

    def slop(self, left, right, considerStrand=False):
        x = self.copy()
        if x.strand == "-" and considerStrand is True:
            x.start = x.start - right
            x.end = x.end + left
        else:
            x.start = x.start - left
            x.end = x.end + right
        if x.end < x.start:
            tmp = x.end
            x.end = x.start
            x.start = tmp
        return x

    def get3pEnd(self):
        if self.strand == "-":
            return self.start
        else:
            return self.end

    def get5pEnd(self):
        if self.strand == "-":
            return self.end
        else:
            return self.start

    def five_prime(self):
        x = self.copy()
        if x.strand == "-":
            x.start = x.end
        else:
            x.end = x.start
        return x

    def three_prime(self):
        x = self.copy()
        if x.strand == "-":
            x.end = x.start
        else:
            x.start = x.end
        return x

    def distanceFrom(self, other, considerStrand=True):
        if not self.chromosome == other.chromosome:
            return np.NaN

        distance1 = self.start - other.end
        distance2 = other.start - self.end
        
        if ((distance1 < 0) & (distance2 < 0)):
            return 0

        elif abs(distance1) <= abs(distance2):
            distance = distance1

        else:
            distance = -distance2

        if considerStrand & (other.strand == "-"):
            distance = -1 * distance
        
        return distance

    def contains(self, other, considerStrand=False):
        if (self.chromosome == other.chromosome) and (self.start <= other.start) and (self.end >= other.end) and ( (not considerStrand) or (self.strand == other.strand) ):
            return True
        else:
            return False

    def copy(self):
        return GenomicRange(self.chromosome, self.start, self.end, self.strand, self.name)



def getFeatureStrand(strand="+"):
    '''Converts +/- strand into +1/-1 strand for FeatureLocation'''
    return int(strand + str(1))




def getPegRNA(chromosome, seqStart, seqEnd, seq, guideStart, guideEnd, guideStrand, variantStart, variantEnd, replacementSeq, pbsLength=13, rtTemplateLength=15):
    '''
    All coordinates are based on chromosome reference, 0-based, as in BED files

    Notes on design from Workshop 11/15:
        Nicking in the "+" direction from nick site works better in general, but "-" ok too
        PE3b is ideal - low indels and also high efficiency.  It nicks closer
            Does PE3b still work for insertion and deletion?  It might not form the heteroduplex and so might not work?
        If not PE3b, nick the bottom strand 50-100 bp downstream 

        What will happen when you have a bunch of guides in a single cell?  Can they replace each others' edits?

        PE2 to PE3 increase in efficiency is higher for deletions

        RT Template Length: 
            - GC content could play a role, but the pattern is not clear
            - Don't have the flap stop at a G (C in the pegRNA 3' extension)
            - Ideally RT should be at least 7bp past the end of the edit

            - LoxP insertion was 44bp insertion + 34bp downstream (78bp total)

        End point of the template is more important than the absolute length

        GoldenGate cloning system makes it easy to do 96 at a time

        Neurons survived for 2 weeks after lenti induction, suggesting it's not toxic

        "HDR mode" in CRISPRESSO2
    '''

    featureStrand = getFeatureStrand(guideStrand)

    sequence = Seq(seq)

    guide = GenomicRange(chromosome, guideStart, guideEnd, guideStrand)
    
    if (guideStart < seqStart) or (guideEnd > seqEnd):
        raise ValueError("Guide bounds are outside sequence bounds")

    guideSeq = (FeatureLocation(guide.start, guide.end, strand=featureStrand)+(-seqStart)).extract(sequence)
    if guideSeq == Seq(''):
        raise ValueError("Guide bounds are outside sequence bounds")

    variant = GenomicRange(chromosome, variantStart, variantEnd, guideStrand)

    if PrimeEditor.GetPositionRelativeToNick(variant.get5pEnd(), guide) < 0:
        raise ValueError("Variant start is behind nick position and cannot be edited with this guide.")

    pbs = guide.slop(pbsLength-NICK_POSITION - guide.length(), NICK_POSITION, considerStrand=True)
    pbsSeq = (FeatureLocation(pbs.start, pbs.end, strand=featureStrand)+(-seqStart)).extract(sequence)

    #    rtRange1 = guide.slop(-guide.length()-NICK_POSITION, variant.get5pEnd()-guide.get3pEnd(), considerStrand=True)
    #    rtRange2 = guide.slop(guide.get5pEnd() - variant.get3pEnd(), variant.get3pEnd()+rtTemplateLength-guide.get3pEnd(), considerStrand=True)
    rtRange1 = guide.slop(-guide.length()-NICK_POSITION, variant.five_prime().distanceFrom(guide.three_prime()), considerStrand=True)
    rtRange2 = guide.slop(guide.five_prime().distanceFrom(variant.three_prime()), variant.three_prime().distanceFrom(guide.three_prime()), considerStrand=True)
    rtRange2 = rtRange2.slop(0, rtTemplateLength, considerStrand=True)

    rtSeq1 = (FeatureLocation(rtRange1.start, rtRange1.end, strand=featureStrand)+(-seqStart)).extract(sequence)
    rtSeq2 = (FeatureLocation(rtRange2.start, rtRange2.end, strand=featureStrand)+(-seqStart)).extract(sequence)

    if (guideStrand == "-"):
        rtTemplate = GenomicRange(chromosome, rtRange2.start, rtRange1.end, guideStrand)
    else:
        rtTemplate = GenomicRange(chromosome, rtRange1.start, rtRange2.end, guideStrand)

    ## Check that last base of rtSeq2 is not G

    if (guideStrand == "-"):
        rtSeq = rtSeq1 + Seq(replacementSeq).reverse_complement() + rtSeq2
    else:
        rtSeq = rtSeq1 + Seq(replacementSeq) + rtSeq2

    results = { "guideSeq": guideSeq,
             "pbsSeq" : pbsSeq.reverse_complement(),
             "rtSeq" : rtSeq.reverse_complement(),
             "pbsSeqInGenome": pbsSeq,
             "rtSeqInGenome": rtSeq }

    return PrimeEditor(results['guideSeq'], results['pbsSeq'], results['rtSeq'], guide, variant, pbs, rtTemplate)


def getPegRNAFromPandas(row):
    return PrimeEditor(
        spacer=row['spacer'],
        primerBindingSite=row['primerBindingSite'], 
        rtTemplate=row['rtTemplate'], 
        guideRegion=GenomicRange(row['chr'],int(row['guideStart']),int(row['guideEnd']),row['guideStrand']), 
        variantRegion=GenomicRange(row['chr'],int(row['variantStart']),int(row['variantEnd']),row['guideStrand']), 
        pbsRegion=GenomicRange(row['chr'],int(row['pbsStart']),int(row['pbsEnd']),row['guideStrand']),
        rtTemplateRegion=GenomicRange(row['chr'],int(row['rtTemplateStart']),int(row['rtTemplateEnd']),row['guideStrand'])
        )


def getAllPegRNAs(
    chromosome, 
    seqStart, 
    seqEnd, 
    seq, 
    guideStart, 
    guideEnd, 
    guideStrand, 
    variantStart, 
    variantEnd, 
    replacementSeq, 
    minPbsLength,
    maxPbsLength,
    minRTPastEdit,
    maxRTPastEdit,
    maxRTTemplateLength):
    
    results = []
    for pbsLength in range(minPbsLength, maxPbsLength+1):
        for rtLengthPastEdit in range(minRTPastEdit, maxRTPastEdit+1):
            try:
                pegRNA = getPegRNA(chromosome, seqStart, seqEnd, seq, guideStart, guideEnd, guideStrand, variantStart, variantEnd, replacementSeq, pbsLength, rtLengthPastEdit)
                if len(pegRNA.rtTemplate) <= maxRTTemplateLength:
                    curr = pegRNA.toPandas()
                    results.append(curr)
            except ValueError:
                continue

    if (len(results) > 0):
        return pd.concat(results, axis=1).transpose()
    else:
        return pd.DataFrame()


def reverseStrand(strand):
    if (strand == "+"):
        return "-"
    elif (strand == "-"):
        return "+"
    else:
        raise ValueError("reverseStrand accepts only + and -")


def getNickingGuides(guides, chromosome, seqStart, seqEnd, seq, guideStart, guideEnd, guideStrand, minPE3NickDistance, maxPE3NickDistance):
    guide = GenomicRange(chromosome, guideStart, guideEnd, guideStrand)
    editnick = guide.three_prime().shift(NICK_POSITION, considerStrand=True)

    ## To do: Speed up code by replacing pd.DataFrame().append() with pd.concat()
    results = pd.DataFrame()
    
    if maxPE3NickDistance > 0:
        for index, row in guides.iterrows():
            nick2 = GenomicRange(row['chr'], row['start'], row['end'], row['strand']).three_prime().shift(NICK_POSITION, considerStrand=True)
            nickDistance = nick2.distanceFrom(editnick)

            if (nick2.chromosome == guide.chromosome) & (nick2.strand == reverseStrand(guideStrand)) & (abs(nickDistance) <= maxPE3NickDistance) & (abs(nickDistance) >= minPE3NickDistance):
                row = row.rename(lambda x: "PE3nick_" + str(x))
                row['PE3nick_distance'] = nickDistance
                results = results.append(row)

    return results

