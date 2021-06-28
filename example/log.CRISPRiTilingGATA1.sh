###########################################################################
## Jesse Engreitz
## 6/14/21
## CRISPRi library targeting all DNase peaks in GATA1 locus in K562, GM12878, and HCT116 cell lines
## This log.sh file lists commands that will need to be modified for future projects,
##  and is meant to be run by copy/pasting each line sequentially into your shell

###########################################################################
## Set key variables for the project
PROJECT=$OAK/Projects/CRISPRDesign/210614_GATA1LocusCRISPRi/; cd $PROJECT
CODEDIR=$PROJECT/CRISPRiTilingDesign/
SIZES=/oak/stanford/projects/genomics-refs/refs/hg19/hg19.chrom.sizes
GENOME_FASTA=/oak/stanford/projects/genomics-refs/refs/hg19/male.hg19.fa
OFF_TARGET_BITS=$OAK/Projects/CRISPR/OffTargets/hg19.CRISPR.bit
PROMOTERS=$GROUP_HOME/Software/ABC-Enhancer-Gene-Prediction/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS500bp.bed


## Download CRISPR designer code
git clone git@github.com:broadinstitute/CRISPRiTilingDesign.git

## If necessary, install the conda environment (already installed on Stanford Sherlock in $GROUP_HOME/Software/anaconda3/envs/):
## conda env create -f CRISPRiTilingDesign/CRISPRDesigner.yml
conda activate EngreitzLab


###########################################################################
## Set up directories
mkdir -p $PROJECT/{01_ChooseRegions,02_RunCRISPRDesigner,03_Subpools,04_CombinePools}


#######################################################################################
## 01_ChooseRegions
## Output of this section:  
##      - A BED file containing chr,start,end,name (no header) listing all of the regions 
##        in which to design gRNAs

ABCLIST=/oak/stanford/groups/engreitz/Projects/ABC/191216_ABC/config/210221.ABCCellTypeTable.ABCPaperV3.LocusPlots.txt
CELLTYPES="K562,HCT116,GM12878"  ## Matching MergedCellType column of $ABCLIST
RANGE=1000000  ## Distance on either side of gene to include DNase peaks
GENES="GATA1"  ## comma-delimited

## This lightweight snakemake pipeline implements the bash code below for making CRISPRi tiling screens:
snakemake \
  -s CRISPRiTilingDesign/SelectRegions/Snakefile \
  --directory 01_ChooseRegions/ \
  --config abcParams=$ABCLIST cellTypes=$CELLTYPES distance=$RANGE genes=$GENES promoters=$PROMOTERS sizes=$SIZES \
  -n 


{
  ## This is unnecessary if running the snakemake pipeline above
  ## Included as an example in case you want to design guides using slightly different rules

  ## Set genes to tile:
  echo $GENES | tr ',' '\n' > $PROJECT/01_ChooseRegions/ChosenGenes.txt

  ## Grab the promoter regions for each of the selected genes
  grep -w -f $PROJECT/01_ChooseRegions/ChosenGenes.txt $PROMOTERS | awk -v OFS=$'\t' '{ $4=$4 "_TSS"; print $0 }'> $PROJECT/01_ChooseRegions/ChosenGenes.Promoters.bed

  ## Create a ChosenGenes.Loci.bed by selecting 1 Mb on either side of the chosen genes
  cat $PROJECT/01_ChooseRegions/ChosenGenes.Promoters.bed | bedtools slop -i stdin -b $RANGE -g $SIZES > $PROJECT/01_ChooseRegions/ChosenGenes.Loci.bed

  ## Grab 300bp regions centered on DNase peak summits from the ABC candidate regions in each of the cell types
  for celltype in $CELLTYPES; do
    ABCDIR=$(cat $ABCLIST |  csvtk grep -t -f MergedCellType -p $celltype | csvtk cut -t -f Neighborhoods | tail -1)
    cat $ABCDIR/../Peaks/macs2_summits.bed | bedtools intersect -u -a stdin -b $PROJECT/01_ChooseRegions/ChosenGenes.Loci.bed | \
      bedtools intersect -a stdin -b $ABCDIR/EnhancerList.bed | \
      bedtools slop -b 150 -i stdin -g $SIZES | \
      bedtools sort -i stdin -faidx $SIZES | \
      cut -f 1-3 | \
      awk -v celltype=$celltype '{ print $1 "\t" $2 "\t" $3 "\t" celltype "|" $1 ":" $2 "-" $3 }' > $PROJECT/01_ChooseRegions/ChosenGenes.${celltype}.summits.bed
  done

  ## Add in full 500-bp promoter regions, not just 300-bp regions around summits, and combine multiple cell types into one file
  cat $PROJECT/01_ChooseRegions/ChosenGenes.Promoters.bed $PROJECT/01_ChooseRegions/ChosenGenes.*.summits.bed | \
    cut -f 1-4 | bedtools sort -i stdin -faidx $SIZES > \
    $PROJECT/01_ChooseRegions/ChosenGenes.AllRegions.bed
}


## To do: Add some code reporting total number of regions and region sizes -- want to avoid any regions that are too big

## Final result: a BED file with chr start end regionName, for input into the CRISPRDesigner
## e.g.:
#chr1    11623   12123   DDX11L1
#chr1    29120   29620   WASH7P
#chr1    35831   36331   FAM138F
#chr1    35831   36331   FAM138A
#chr1    68840   69340   OR4F5

#######################################################################################
## 02_RunCRISPRDesigner

## Checkout the snakemake pipeline above
cd $PROJECT/02_RunCRISPRDesigner/
git clone git@github.com:EngreitzLab/CRISPRDesigner.git
cd $PROJECT/02_RunCRISPRDesigner/CRISPRDesigner/
mkdir logs

## Set up snakemake pipeline by editing CRISPRDesigner/config/config.yaml
## Then, run the CRISPRDesigner snakemake pipeline locally:
snakemake &

## Or, to run on cluster for a long job with lots of gRNA scoring:
snakemake \
  --config java_memory=16g \
  --cores 1 \
  --jobs 50 \
  --cluster "sbatch -n 1 -c 1 --mem 16G -t 12:00:00 -p owners -J CRISPRDesigner_{rule} -o logs/{rule}_{wildcards} -e logs/{rule}_{wildcards}" &

cp $PROJECT/02_RunCRISPRDesigner/CRISPRDesigner/results/GuideDesign/designGuides.txt $PROJECT/02_RunCRISPRDesigner/designGuides.txt
cd $PROJECT

#######################################################################################
## 03_Subpools

## First, edit and copy in the list of guides to force inclusion (e.g., previously validated gRNAs for GATA1)
## Here, for GATA1 locus:
cp $PROJECT/../CRISPR_Screen/200115_GATA1HyPR/01_WorkingGuides/GuidesToInclude.txt $PROJECT/03_Subpools

python $CODEDIR/src/MakeGuidePool.py \
      --PoolID 210615_GATA1CRISPRiTiling \
      --input $PROJECT/02_RunCRISPRDesigner/designGuides.txt \
      --outdir $PROJECT/03_Subpools/ \
      --nCtrls 463 \
      --negCtrlList $CODEDIR/data/Weissman1000.negative_control.20bp.design.txt \
      --nSafeCtrls 200 \
      --safeCtrlList $CODEDIR/data/Tycko2019SafeTargeting.txt \
      --nGuidesPerElement 15 \
      --selectMethod even \
      --trimElements 5 \
      --vector Crop-Opti \
      --vectorDesigns $CODEDIR/data/CloningDesigns.txt \
      --forceGuides $PROJECT/03_Subpools/GuidesToInclude.txt


######################################################################################
## 04_CombinePools
## First:  Edit CombinePoolsConfig.txt, and pull primer pairs from $CODEDIR/SubpoolPrimers.txt

mkdir -p $PROJECT/04_CombinePools/

dos2unix $PROJECT/04_CombinePools/CombinePoolsConfig.txt
python $CODEDIR/src/CombineGuidePools.py \
    --config $PROJECT/04_CombinePools/CombinePoolsConfig.txt \
    --outbase $PROJECT/04_CombinePools/210615_GATA1CRISPRiTiling  \
    > $PROJECT/04_CombinePools/210615_GATA1CRISPRiTiling.out 

