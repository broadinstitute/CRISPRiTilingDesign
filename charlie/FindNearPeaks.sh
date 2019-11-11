

# Finds all the peaks within radius of gene (mostly just runs bedtools window)
# USAGE: FindNearPeaks.sh GENES PEAKS RADIUS OUTFILE

RADIUS=$3
GENES=$1
PEAKS=$2
OUTFILE=$4

cut -f1,2,3,4 $GENES |\
bedtools window -w $RADIUS -a stdin -b $PEAKS |\
awk -v OFS="\t" '{print $5,$6,$7,$4}' > $OUTFILE


