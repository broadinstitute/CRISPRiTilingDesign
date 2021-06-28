
# For regions with the guides already desined, grab the guides and report any regions that lack guides

# USAGE: GetGuides TARGETS GUIDES OUTFILE

TARGETS=$1
GUIDES=$2
OUTFILE=$3

# intersect to get the guides in the target regions
bedtools intersect -wa -a $GUIDES -b $TARGETS > $OUTFILE

# report the number of guides per target region
bedtools intersect -c -wa -a $TARGETS -b $GUIDES > $OUTFILE.count 

