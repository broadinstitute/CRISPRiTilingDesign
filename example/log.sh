## Jesse Engreitz
## 10/30/19
## Design for CRISPRi pool using dead gRNAs with Cas9 mouse

## Gene list and other overall design info: https://docs.google.com/document/d/1aJNGAz5ExeedPdJVLvgp6V3jktqP1Z1qACwnv5uWPFA/edit
## Eye-balled some enhancers in the regions (but am certain that ABC would find them significant)
##   Ly96 does not have good enhancer cnadidates
##   Tlr4 has some good candidates too but not chosen due to space constraints
##   > ChosenEnhancers.bed

## Set up GeneList.txt


###########################################################################
## Parameters to set by command line in future
PROJECT=/seq/lincRNA/Jesse/CRISPR_Screen/191030_BMDC_Cas9_CRISPRi/
GENELIST=$PROJECT/01_ChooseRegions/GeneList.txt

CODEDIR=$PROJECT/CRISPRiTilingDesign/

###########################################################################
## Parameters for genome to factor out in future

## mm9
SIZES=/seq/lincRNA/data/mm9/sizes
OFF_TARGET_BITS=/seq/lincRNA/CRISPR/OffTargets/mm9.CRISPR.sorted.bit
GENOME_FASTA=/seq/lincRNA/data/mm9/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa
TSS500BP=/seq/lincRNA/data/mm9/RefSeqCurated.180427.bed.CollapsedGeneBounds.TSS.500bp.bed

## hg19 
#SIZES=/seq/lincRNA/data/hg19/sizes
#OFF_TARGET_BITS=/seq/lincRNA/CRISPR/OffTargets/hg19.CRISPR.bit
#GENOME_FASTA=/seq/lincRNA/data/hg19/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
#TSS500BP=/seq/lincRNA/data/hg19/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS.500bp.bed


## Key scripts
use BEDTools


#######################################################################################
## 01_ChooseRegions
## Get promoter regions
join -1 4 -2 1 -o 1.1,1.2,1.3,1.4  \
 <(sort -k 4 $TSS500BP | bedtools slop -s -l -150 -r 50 -g $SIZES) \
 <(sort $GENELIST) | tr ' ' '\t' \
 > $PROJECT/01_ChooseRegions/GenePromoters.bed


TARGETS=$PROJECT/01_ChooseRegions/AllRegions.bed
cat $PROJECT/01_ChooseRegions/GenePromoters.bed $PROJECT/01_ChooseRegions/ChosenEnhancers.bed > $TARGETS


#######################################################################################
## SIDE STEP:  PICK CODING GUIDES
use .python-3.5.1; source /seq/lincRNA/Ben/VENV_MIP/bin/activate
python $CODEDIR/src/GetCodingGuides.py \
        --genes ../01_ChooseRegions/GeneList.txt \
        --outfile ../01_ChooseRegions/GeneList.CodingGuides.txt \
        --guideLibrary $CODEDIR/data/broadgpp-brie-library-contents.txt.gz


#######################################################################################
## 02_RunCRISPRDesigner

DIR=$PROJECT/02_RunCRISPRDesigner/
mkdir -p $DIR $DIR/log/ $DIR/design/

CMD="reuse Java-1.7; java -Xmx12g -jar $CODEDIR/src/CRISPRDesigner.jar TARGETS=$TARGETS OUTPUT_DIR=$DIR/design/ GENOME_FASTA=$GENOME_FASTA LENIENT=false OFF_TARGETS=$OFF_TARGET_BITS SKIP_PAIRING=true DIVIDE_AND_CONQUER=true MAX_DIVIDE=600 GRID_ENGINE=true QUEUE=gsa"
quick-qsub -o $DIR/log/o.runCD.qq -s $DIR/log/s.runCD.qq $CMD

## Charlie's script
## Parameters to consider modifying:  MinStartDistance (i.e., minimum spacing between gRNA start sites)
## Add that parameter at top to main script
reuse Python-2.7  # this is needed for Charlie's scripts
$CODEDIR/charlie/GetTileGuides.py \
  --infile $DIR/design/filteredGuides.bed \
  --outprefix $DIR/design/filteredGuides.bed \
  --chrom All --start 0 --end 0 -T 4 --OTMin 50 --MinStartDistance 5

## Make version openable in IGV for looking at locations
sed 1d $DIR/design/filteredGuides.bed.preDesign.bed > $DIR/design/filteredGuides.bed.getGuides.thin.bed


#######################################################################################
## 03_Subpools

mkdir $PROJECT/03_Subpools/

python $CODEDIR/src/MakeGuidePool.py \
    --PoolID CRISPRiFullGuide \
    --input $PROJECT/02_RunCRISPRDesigner/design/filteredGuides.bed.preDesign.bed \
    --outdir $PROJECT/03_Subpools/ \
    --nCtrls 50 \
    --negCtrlList $CODEDIR/data/Weissman1000.negative_control.20bp.design.txt \
    --nGuidesPerElement 15 \
    --vector sgMS2 \
    --vectorDesigns $CODEDIR/data/CloningDesigns.txt \

python $CODEDIR/src/MakeGuidePool.py \
    --PoolID CodingFullGuide \
    --input $PROJECT/01_ChooseRegions/GeneList.CodingGuides.txt \
    --outdir $PROJECT/03_Subpools/ \
    --nGuidesPerElement 4 \
    --vector sgMS2 \
    --vectorDesigns $CODEDIR/data/CloningDesigns.txt 


## Make dead guides for sgMS2-KRAB test
    python $CODEDIR/src/GetDeadGuides.py \
        --design $PROJECT/03_Subpools/CRISPRiFullGuide.design.txt \
        --outfile $PROJECT/03_Subpools/CRISPRiDeadGuide.design.txt \
        --PoolID CRISPRiDeadGuide

python $CODEDIR/src/GetDeadGuides.py \
    --design $PROJECT/03_Subpools/CodingFullGuide.design.txt \
    --outfile $PROJECT/03_Subpools/CodingDeadGuide.design.txt \
    --PoolID CodingDeadGuide


## 11/11/19 TODO:

## TO do (much of which appears to be in Charlie's script)
##   clean up code and move to github
##   Determine whether we are going to subpool this with other guides 
##   Write script to combine multiple subpools and add handles
##   CHECK and QC!

