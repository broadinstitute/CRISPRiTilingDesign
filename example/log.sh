#!/bin/bash

## Run for all Genes in hg38 list


###########################################################################
## Parameters to set by command line in future
PROJECT=/mnt/DATA/GuideScanTest
GENELIST=GeneList.txt
CODEDIR=/mnt/DATA/CRISPR_design_scripts/CRISPRiTilingDesign/
GUIDESCAN_INDEX=/mnt/DATA/CRISPR_design_scripts/guidescan_index
NGUIDESPERELEMENT=15


###########################################################################
## Parameters for genome to factor out in future

## hg38
## TODO - remove hardcoded paths
SIZES=/mnt/DATA/CRISPR_design_scripts/hg38/sizes
TSS500BP=/mnt/DATA/CRISPR_design_scripts/hg38/RefSeqCurated.200730.UniqueTSS.CollapsedGeneBounds.500bp.bed

# get Python 3.7 or higher
eval "$(conda shell.bash hook)"
conda activate guidescan_env


#######################################################################################
## 01_ChooseRegions
## Get promoter regions
mkdir -p $PROJECT/01_ChooseRegions

join -1 4 -2 1 -o 1.1,1.2,1.3,1.4  \
 <(sort -k 4 $TSS500BP | bedtools slop -s -l -150 -r 50 -g $SIZES) \
 <(sort $GENELIST) | tr ' ' '\t' \
 > $PROJECT/01_ChooseRegions/GenePromotersLocations.bed


# add unique designator to multi-TSS genes
python $CODEDIR/src/bed_add_unique_suffix.py TSS $PROJECT/01_ChooseRegions/GenePromotersLocations.bed $PROJECT/01_ChooseRegions/GenePromoters.bed

TARGETS=$PROJECT/01_ChooseRegions/AllRegions.bed
cat $PROJECT/01_ChooseRegions/GenePromoters.bed $PROJECT/01_ChooseRegions/ChosenEnhancers.bed > $TARGETS


#######################################################################################
## SIDE STEP:  PICK CODING GUIDES
python $CODEDIR/src/GetCodingGuides.py \
        --genes $GENELIST \
        --outfile $PROJECT/01_ChooseRegions/GeneList.CodingGuides.txt \
        --guideLibrary $CODEDIR/data/broadgpp-brie-library-contents.txt.gz


#######################################################################################
## 02_RunGuideScan via docker


GS_COMMAND="/usr/local/bin/guidescan_guidequery \
    -b /mnt/guidescan/cas9_hg38_all_guides.bam \
    --target within \
    -o /mnt/project \
    -n 20 \
    --batch /mnt/project/01_ChooseRegions/AllRegions.bed \
    --annot /mnt/project/01_ChooseRegions/AllRegions.bed \
    --select score"

docker run \
       -v ${GUIDESCAN_INDEX}:/mnt/guidescan \
       -v $PROJECT:/mnt/project \
       guidescan \
       $GS_COMMAND

## Charlie's script - removed by Ray August 4, 2020
##
## Charlie's sript removed any guides with 4 or more poly T,
## Off-Target min 50 score (not sure how that compares betwen
## CRISPRDesigner and Guidescan), and made sure they were separated by
## at least 5bp in their start location.
DIR=$PROJECT/02_RunCRISPRDesigner/
mkdir -p $DIR $DIR/design/

python $CODEDIR/src/filter_guides.py \
  --infile $PROJECT/GuideScan_batch_output.csv \
  --outfile $DIR/design/filteredGuides.bed.preDesign.bed \
  --polyTmax 4 --minStartDistance 5


## Make version openable in IGV for looking at locations
sed 1d $DIR/design/filteredGuides.bed.preDesign.bed > $DIR/design/filteredGuides.bed.getGuides.thin.bed


#######################################################################################
## 03_Subpools

mkdir -p $PROJECT/03_Subpools/

echo Processing Target Guides...
python $CODEDIR/src/MakeGuidePool.py \
    --PoolID CRISPRiFullGuide \
    --input $PROJECT/02_RunCRISPRDesigner/design/filteredGuides.bed.preDesign.bed \
    --outdir $PROJECT/03_Subpools/ \
    --nCtrls 50 \
    --negCtrlList $CODEDIR/data/Weissman1000.negative_control.20bp.design.txt \
    --nGuidesPerElement $NGUIDESPERELEMENT \
    --vector sgOpti \
    --seqCol gRNA \
    --noPAM \
    --vectorDesigns $CODEDIR/data/CloningDesigns.txt

echo ""
echo Processing Coding Guides...
python $CODEDIR/src/MakeGuidePool.py \
    --PoolID CodingFullGuide \
    --input $PROJECT/01_ChooseRegions/GeneList.CodingGuides.txt \
    --outdir $PROJECT/03_Subpools/ \
    --nGuidesPerElement 4 \
    --vector sgOpti \
    --vectorDesigns $CODEDIR/data/CloningDesigns.txt 
echo ""

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

