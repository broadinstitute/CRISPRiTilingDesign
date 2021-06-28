
# PickGuides
# USAGE: PickGuides.sh TARGETS PROJECTDIR

TARGETS=$1
PROJECT=$2


CODEDIR=/seq/lincRNA/cfulco/bin/scripts/
CDDIR=/seq/lincRNA/Jesse/bin/scripts/

## Run CRISPR designer
OFF_TARGET_BITS=/seq/lincRNA/CRISPR/OffTargets/hg19.CRISPR.bit
GENOME_FASTA=/seq/lincRNA/data/hg19/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa

mkdir -p $PROJECT/log/

CMD="java -Xmx12g -jar $CDDIR/CRISPRDesigner.jar TARGETS=$TARGETS OUTPUT_DIR=$PROJECT/design/ GENOME_FASTA=$GENOME_FASTA LENIENT=true OFF_TARGETS=$OFF_TARGET_BITS SKIP_PAIRING=true DIVIDE_AND_CONQUER=true MAX_DIVIDE=1000 GRID_ENGINE=true QUEUE=gsa"

NAME=$(basename $TARGETS .bed)
$CODEDIR/quick-qsub-gsa -o $PROJECT/log/o.$NAME.PICKGUIDES.qq -s $PROJECT/log/s.$NAME.PickGuides.qq $CMD



