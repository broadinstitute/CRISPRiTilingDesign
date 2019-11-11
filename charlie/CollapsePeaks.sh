
# USAGE: CollapsePeaks.sh PEAKLIST PREFIX BUFFER MINSIZE
# PEAKLIST is txt file of peak bed files or a bedfile itself

CODEDIR=/seq/lincRNA/cfulco/bin/scripts/
SIZES=/seq/lincRNA/cfulco/Data/hg19/hg19.sizes
PROJECTDIR=.

PEAKLIST=$1
PREFIX=$2
BUFFER=$3
MINSIZE=$4
EXTRA=$5

# ----------------------------------------------------------#
# 1. Make list of all potential peaks

# count columns in PEAKLIST
COLS=$(awk 'NR==1{print NF;exit}' $PEAKLIST)

# 1a. PEAKLIST is list of files. assemble into single bed
if [ $COLS == 1 ]; then 
$CODEDIR/catFilesInFile.sh $PEAKLIST $PREFIX.allPeaks.bed.tmptmp
fi

# 1b. PEAKLIST is already a single bed
if [ $COLS -ge 3 ]; then
cp $PEAKLIST $PREFIX.allPeaks.bed.tmptmp
fi

# ----------------------------------------------------------#
# 2. Merge the peaks (could come from close together peaks or the same peak in different cell types)
cat $PREFIX.allPeaks.bed.tmptmp | cut -f 1-3 |\
bedtools slop -i stdin -b $EXTRA -g $SIZES |\
bedtools slop -i stdin -b $BUFFER -g $SIZES |\
bedtools sort -i stdin |  bedtools merge -i stdin |\

# ----------------------------------------------------------#
# 3. Remove buffer but enforce the minimum size
awk -v OFS="\t" -v BUFFER=$BUFFER -v MINSIZE=$MINSIZE \
'{ $3 = $3 - BUFFER; $2 = $2 + BUFFER; if ($3 - $2 < MINSIZE) { adjust=(MINSIZE-($3-$2))/2; $2 = int($2 - adjust); $3 = int($3 + adjust) } print $0 }' > $PREFIX.merged.min$MINSIZE.bed

# ----------------------------------------------------------#

rm $PREFIX.allPeaks.bed.tmptmp

