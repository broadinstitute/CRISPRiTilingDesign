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
GCP_PROJECT=landerlab-atacseq-200218
PROJECT=/seq/lincRNA/thouis/GuideScanTest
GENELIST=GeneList.txt
CODEDIR=$PROJECT/CRISPRiTilingDesign/

LOGGING_BUCKET=gs://landerlab-guidescan/logging
GUIDESCAN_BUCKET=gs://landerlab-guidescan-index
TMP_BUCKET=gs://landerlab-guidescan/test

###########################################################################
## Parameters for genome to factor out in future

## hg38
SIZES=/seq/lincRNA/data/hg38/sizes
TSS500BP=/seq/lincRNA/data/hg38/RefSeqCurated.200730.UniqueTSS.CollapsedGeneBounds.500bp.bed

## Key scripts
use BEDTools

#######################################################################################
## 01_ChooseRegions
## Get promoter regions
mkdir -p $PROJECT/01_ChooseRegions

join -1 4 -2 1 -o 1.1,1.2,1.3,1.4  \
 <(sort -k 4 $TSS500BP | bedtools slop -s -l -150 -r 50 -g $SIZES) \
 <(sort $GENELIST) | tr ' ' '\t' \
 > $PROJECT/01_ChooseRegions/GenePromotersLocations.bed


# add unique designator to multi-TSS genes
source /seq/lincRNA/thouis/VENV36/bin/activate
python $CODEDIR/src/bed_add_unique_suffix.py TSS $PROJECT/01_ChooseRegions/GenePromotersLocations.bed $PROJECT/01_ChooseRegions/GenePromoters.bed



TARGETS=$PROJECT/01_ChooseRegions/AllRegions.bed
cat $PROJECT/01_ChooseRegions/GenePromoters.bed $PROJECT/01_ChooseRegions/ChosenEnhancers.bed > $TARGETS


#######################################################################################
## SIDE STEP:  PICK CODING GUIDES
source /seq/lincRNA/thouis/VENV36/bin/activate
python $CODEDIR/src/GetCodingGuides.py \
        --genes $GENELIST \
        --outfile $PROJECT/01_ChooseRegions/GeneList.CodingGuides.txt \
        --guideLibrary $CODEDIR/data/broadgpp-brie-library-contents.txt.gz


#######################################################################################
## 02_RunGuideScan

reuse -q Google-Cloud-SDK
gsutil cp $TARGETS $TMP_BUCKET/targets.bed

/seq/lincRNA/thouis/VENV36/bin/dsub --project $GCP_PROJECT \
     --zones "us-*" \
     --logging $LOGGING_BUCKET \
     --input GS_BAM=$GUIDESCAN_BUCKET/cas9_hg38_all_guides.bam \
     --input GS_BAI=$GUIDESCAN_BUCKET/cas9_hg38_all_guides.bam.bai \
     --input TARGETS=$TMP_BUCKET/targets.bed \
     --output GUIDES_CSV=$TMP_BUCKET/guidescan_guides.csv \
     --image us.gcr.io/$GCP_PROJECT/guidescan \
     --wait \
     --command "mkdir guidescan_outdir; \
                guidescan_guidequery \
                -b \$GS_BAM \
                --target within \
                -o guidescan_outdir \
                --batch \$TARGETS -n 8 \
                --annot \$TARGETS \
                --select score; \
                cp guidescan_outdir/*.csv \$GUIDES_CSV"

gsutil cp $TMP_BUCKET/guidescan_guides.csv $PROJECT/guidescan_guides.csv
     
## Charlie's script - removed by Ray August 4, 2020
##
## Charlie's sript removed any guides with 4 or more poly T,
## Off-Target min 50 score (not sure how that compares betwen
## CRISPRDesigner and Guidescan), and made sure they were separated by
## at least 5bp in their start location.
DIR=$PROJECT/02_RunCRISPRDesigner/
mkdir -p $DIR $DIR/design/

python $CODEDIR/src/FilterGuides.py \
  --infile $PROJECT/guidescan_guides.csv \
  --outfile $DIR/design/filteredGuides.bed.preDesign.bed \
  --polyTmax 4 --minStartDistance 5


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

