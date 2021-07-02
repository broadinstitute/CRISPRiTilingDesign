###########################################################################
## Jesse Engreitz
## 7/1/21
## Prime editor gRNAs targeting a few positive control sites near GATA1 and PPIF, selected by Gabbie
## This log file lists commands that will need to be modified for future projects,
##  and is meant to be run by copy/pasting each line sequentially into your shell

###########################################################################
## Set key variables for the project
PROJECT=$OAK/Projects/CRISPRDesign/210701_GATA1PPIFpegRNAs/; cd $PROJECT
CODEDIR=$PROJECT/CRISPRiTilingDesign/
REF=/oak/stanford/projects/genomics-refs/refs/
SIZES=$REF/hg19/hg19.chrom.sizes
GENOME_FASTA=$REF/hg19/male.hg19.fa
OFF_TARGET_BITS=$OAK/Projects/CRISPR/OffTargets/hg19.CRISPR.bit

mkdir -p $PROJECT/{01_ChooseRegions,log,02_RunCRISPRDesigner,03_DesignPrimeEditor}

## First time only:
git clone git@github.com:broadinstitute/CRISPRiTilingDesign.git

## First time running pipeline: 
## conda env create -f CRISPRiTilingDesign/CRISPRDesigner.yml
## This should now be installed in $GROUP_HOME/Software/anaconda3/envs/
conda activate CRISPRDesigner

## Set up Variants.txt file

#######################################################################################
## 01_ChooseRegions

REGIONS=$PROJECT/01_ChooseRegions/Variants.50bp.bed

## Take the variants and expand by 50bp on either side for finding potential spacers/pegRNAs
cd $PROJECT/01_ChooseRegions/
sed 1d GabbieVariants.txt | bedtools slop -b 50 -g $SIZES -i stdin | cut -f 1-3,7 | sort | uniq > $REGIONS


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
  --config java_memory=15g \
  --cores 1 \
  --jobs 50 \
  --cluster "sbatch -n 1 -c 1 --mem 16G -t 12:00:00 -p owners -J CRISPRDesigner_{rule} -o logs/{rule}_{wildcards} -e logs/{rule}_{wildcards}" &

cp $PROJECT/02_RunCRISPRDesigner/CRISPRDesigner/results/GuideDesign/designGuides.txt $PROJECT/02_RunCRISPRDesigner/designGuides.txt
cd $PROJECT



#######################################################################################
## 03_DesignPrimeEditor

DIR_PE=$PROJECT/03_DesignPrimeEditor/

python $CODEDIR/src/DesignPrimeEditor.py \
  --edits $PROJECT/01_ChooseRegions/GabbieVariants.txt \
  --outfile $DIR_PE/RT30/RT30 \
  --guides $PROJECT/02_RunCRISPRDesigner/designGuides.txt \
  --fasta $GENOME_FASTA \
  --minPbsLength=11 --maxPbsLength=11 \
  --minRTPastEdit=10 --maxRTPastEdit=10 \
  --maxRTTemplateLength 30 \
  --minPbsGcContent 30 --minRTGcContent 30 \
  --maxUUcount 4 \
  --vectorDesigns $PROJECT/CRISPRiTilingDesign/data/CloningDesigns.txt \
  --vector pegRNAsgOptiGibsonOrigScaffold \
  --ignoreRegions



#######################################################################################
## 04_SelectPegRNAs
## To do:  Stats about # variants covered, length of PBS + RT, etc.

conda activate EngreitzLab

## Select pegRNAs across entire (large, expanded by 500bp) regions
## Ran this originally -- not used anymore
python $CODEDIR/src/SelectPrimeEditorGuides.py \
  -i $DIR_PE/RT30/RT30.full.tsv \
  -i $DIR_PE/RT50/RT50.full.tsv \
  -i $DIR_PE/RT50-LenientGuides/RT50Lenient.full.tsv \
  --edits $PROJECT/01_ChooseRegions/Amplicons200115.plus500.tile5.txt \
  -o $DIR_PE/merged/merged \
  --guidesPerVariant 1 \
  --spacerChoice mixed \
  --splitOutputByRegion


## Select pegRNAs for the selected PCR amplicons 
python $CODEDIR/src/SelectPrimeEditorGuides.py \
  -i $DIR_PE/RT30/RT30.full.tsv \
  -i $DIR_PE/RT50/RT50.full.tsv \
  -i $DIR_PE/RT50-LenientGuides/RT50Lenient.full.tsv \
  --edits $PROJECT/01_ChooseRegions/Amplicons210205.tile5.txt \
  --strictPegRegions $PROJECT/01_ChooseRegions/Amplicons210124.GuideRegions.bed \
  -o $DIR_PE/Amplicons210205/Amplicons210205 \
  --guidesPerVariant 1 \
  --spacerChoice mixed \
  --splitOutputByRegion


## Repeat for sgOpti scaffold
python $CODEDIR/src/SelectPrimeEditorGuides.py \
  -i $DIR_PE/RT30-sgOpti/RT30.full.tsv \
  -i $DIR_PE/RT50-sgOpti/RT50.full.tsv \
  -i $DIR_PE/RT50-LenientGuides-sgOpti/RT50Lenient.full.tsv \
  --edits $PROJECT/01_ChooseRegions/Amplicons210205.tile5.txt \
  --strictPegRegions $PROJECT/01_ChooseRegions/Amplicons210124.GuideRegions.bed \
  -o $DIR_PE/Amplicons210205-sgOpti/Amplicons210205 \
  --guidesPerVariant 1 \
  --spacerChoice mixed \
  --splitOutputByRegion


## Select pegRNAs for the single-nucleotide variants
python $CODEDIR/src/SelectPrimeEditorGuides.py \
  -i $DIR_PE/VariantsRT30/VariantsRT30.full.tsv \
  -i $DIR_PE/VariantsRT30-Lenient/VariantsRT30-Lenient.full.tsv \
  --edits $PROJECT/01_ChooseRegions/Variants.plus30.mutagenesis.txt \
  --strictPegRegions $PROJECT/01_ChooseRegions/Amplicons210124.GuideRegions.bed \
  -o $DIR_PE/Variants/Variants \
  --guidesPerVariant 1 \
  --spacerChoice closest \
  --splitOutputByRegion


## Repeat for sgOpti scaffold
python $CODEDIR/src/SelectPrimeEditorGuides.py \
  -i $DIR_PE/VariantsRT30-sgOpti/VariantsRT30.full.tsv \
  -i $DIR_PE/VariantsRT30-Lenient-sgOpti/VariantsRT30-Lenient.full.tsv \
  --edits $PROJECT/01_ChooseRegions/Variants.plus30.mutagenesis.txt \
  --strictPegRegions $PROJECT/01_ChooseRegions/Amplicons210124.GuideRegions.bed \
  -o $DIR_PE/Variants-sgOpti/Variants \
  --guidesPerVariant 1 \
  --spacerChoice closest \
  --splitOutputByRegion




######################################################################################
## 04_CombinePools
## First:  Edit CombinePoolsConfig.txt, and pull primer pairs from $CODEDIR/SubpoolPrimers.txt

mkdir -p $PROJECT/04_CombinePools/

dos2unix -c mac $PROJECT/04_CombinePools/CombinePoolsConfig.txt
python $CODEDIR/src/CombineGuidePools.py \
    --config $PROJECT/04_CombinePools/CombinePoolsConfig.txt \
    --outbase $PROJECT/04_CombinePools/210205_IL2RA_PPIF \
    --fillToOligoPoolSize 7500 > $PROJECT/04_CombinePools/210205_IL2RA_PPIF.out ## Smallest size from Agilent

## Agilent array sizes: 7.5K, 15K, 60K, 100K, 244K

