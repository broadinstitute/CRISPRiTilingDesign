# CRISPR Tiling Design - Engreitz Lab
Scripts for designing CRISPR tiling screens, including CRISPRi and prime editing.  

## Authors

* Jesse Engreitz (@engreitz)
* Please cite: [Fulco\*, Nasser\* et al. Nature Genetics 2019](https://www.nature.com/articles/s41588-019-0538-0)

## Description

This repository includes scripts used to design CRISPRi and prime editing tiling screens for the Engreitz Lab.

The basic steps in the CRISPRi design process are:  
1. Choosing regions (possibly overlapping) to tile with CRISPRi gRNAs  
2. Designing all possible gRNAs in those regions.  
3. Scoring and filtering gRNAs to eliminate off-target effects and poor efficacy gRNAs.  
4. Choosing gRNAs within each region, e.g. by even selection across the region.  
5. Collating different sets of gRNAs ("subpools") into a single oligo "pool", including oligo sequences for ordering

The repository also includes helper scripts for designing dead gRNAs, prime editing gRNAs,
selecting protein-coding gRNAs from Brie/Brunello, and annotating gRNAs from previous designers. 

Depends on CRISPR design code repository using the MIT specificity score, located [here](https://github.com/EngreitzLab/CRISPRDesigner/tree/master).

## Setup

### Step 1: Clone this github repository

[Clone](https://help.github.com/en/articles/cloning-a-repository) this to your local system, into the place where you want to perform the data analysis.

Example:

	git clone git@github.com:broadinstitute/CRISPRiTilingDesign.git

### Step 2: Install conda environment if needed

Install Snakemake and conda environment using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    conda env create --file envs/EngreitzLab.yml  
    conda activate EngreitzLab

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

Engreitz Lab members working on Sherlock can skip this step.


## Usage (CRISPRi)

A full example workflow for designing CRISPRi guides tiling across all DNase peaks around GATA1 is provided in `example/log.CRISPRiTilingGATA1.sh`

### Step 1: Select regions in which to design gRNAs

The first step is to generate a BED file containing the regions in which we want to design gRNAs.  

Do this according to the needs of your experiment. The end result should be BED4 file with `chr	start	end	name`.

For CRISPRi-FlowFISH experiments, you can use the provided Snakemake workflow to select DNase peaks from specific cell types around specified genes (this excerpt points to files on Engreitz Lab cluster; for others, the relevant files are included in the `data/` and `example/` directories of this repository):

	PROMOTERS=$GROUP_HOME/Software/ABC-Enhancer-Gene-Prediction/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS500bp.bed
	ABCLIST=/oak/stanford/groups/engreitz/Projects/ABC/191216_ABC/config/210221.ABCCellTypeTable.ABCPaperV3.LocusPlots.txt
	CELLTYPES="K562,HCT116,GM12878"  ## Matching MergedCellType column of $ABCLIST
	RANGE=1000000  ## Distance on either side of gene to include DNase peaks
	GENES="GATA1"  ## comma-delimited
	snakemake \
	  -s CRISPRiTilingDesign/SelectRegions/Snakefile \
	  --directory 01_ChooseRegions/ \
	  --config abcParams=$ABCLIST cellTypes=$CELLTYPES distance=$RANGE genes=$GENES promoters=$PROMOTERS sizes=$SIZES \
	  -n 

For CRISPRi screens targeting promoters, you can do something simpler, like grabbing the regions from our ABC promoter region files:

	cat $PROMOTERS | csvtk grep -H -t -f 4 -P 01_ChooseRegions/ChosenGenes.txt > 01_ChooseRegions/ChosenGenes.Promoters.txt
	#Output:
	#chrX    48644731        48645231        GATA1   0       +

### Step 2: Design and score gRNAs

The second step involves sending the chosen regions to a CRISPR design tool, which will enumerate all NGG PAM sites and apply a specificity score filter.

To do so, first clone the [CRISPRDesigner repository](https://github.com/EngreitzLab/CRISPRDesigner/tree/master):

	git clone git@github.com:EngreitzLab/CRISPRDesigner.git
	cd CRISPRDesigner/
	mkdir logs

Set up snakemake pipeline by editing CRISPRDesigner/config/config.yaml to point to your region file.  
Optionally, provide a list of pre-scored gRNAs to speed the design process (see [here](https://github.com/EngreitzLab/CRISPRDesigner/tree/master/resources)).

Then, run the CRISPRDesigner snakemake pipeline:

	snakemake 


### Step 3: Select gRNAs and format for oligo pool synthesis

Example:

	python $CODEDIR/src/MakeGuidePool.py \
	      --PoolID 210615_GATA1CRISPRiTiling \
	      --input 02_RunCRISPRDesigner/designGuides.txt \
	      --outdir 03_Subpools/ \
	      --nCtrls 463 \
	      --negCtrlList data/Weissman1000.negative_control.20bp.design.txt \
	      --nSafeCtrls 200 \
	      --safeCtrlList data/Tycko2019SafeTargeting.txt \
	      --nGuidesPerElement 15 \
	      --selectMethod even \
	      --trimElements 5 \
	      --vector Crop-Opti \
	      --vectorDesigns data/CloningDesigns.txt \
	      --forceGuides 03_Subpools/GuidesToInclude.txt


Documentation:

	usage: MakeGuidePool.py [-h] --input INPUT --outdir OUTDIR --nGuidesPerElement
                        NGUIDESPERELEMENT [--negCtrlList NEGCTRLLIST]
                        [--nCtrls NCTRLS] [--safeCtrlList SAFECTRLLIST]
                        [--nSafeCtrls NSAFECTRLS]
                        [--nCtrls NCTRLS] [--safeCtrlList SAFECTRLLIST]
                        [--nSafeCtrls NSAFECTRLS]
                        [--vectorDesigns VECTORDESIGNS] [--vector VECTOR]
                        [--PoolID POOLID] [--seqCol SEQCOL] [--noPAM]
                        [--PAM PAM] [--selectMethod {score,even}]
                        [--forceGuides FORCEGUIDES]
                        [--trimElements TRIMELEMENTS]
                        [--excludeRestrictionSites EXCLUDERESTRICTIONSITES]
                        [--barcodes BARCODES]

	Select gRNAs in regions, add PCR handles for cloning, and create sequences for oligo pool synthesis.

	Optional arguments:
	  -h, --help            show this help message and exit
	  --input INPUT         preDesign.bed file from GetTileGuides.py (default: None)
	  --outdir OUTDIR
	  --nGuidesPerElement NGUIDESPERELEMENT
	                        Number of guides to choose per guide set ('guideSet' column). Set to 0 to select all gRNAs. (default: None)
	  --negCtrlList NEGCTRLLIST
	                        Formatted guide file of negative controls. (default: data/Weissman1000.negative_control.20bp.design.txt)
	  --nCtrls NCTRLS       Number of negative control gRNAs to add to the pool (default: 0)
	  --safeCtrlList SAFECTRLLIST
	                        Formatted guide file for safe-targeting controls. (default: data/Tycko2019SafeTargeting.txt)
	  --nSafeCtrls NSAFECTRLS
	                        Number of safe-targeting gRNAs to add to the pool (default: 0)
	  --vectorDesigns VECTORDESIGNS
	                        Reference table with gibson arm sequences for various plasmid designs (default: data/CloningDesigns.txt)
	  --vector VECTOR       Name of vector to index into vectorDesigns (default: sgOpti)
	  --PoolID POOLID       Unique name of pool or subpool for naming oligos - e.g. 191110_GATA1 (default: MyPool)
	  --seqCol SEQCOL       Name of column with guide sequence (default: GuideSequenceWithPAM)
	  --noPAM               Whether to trim PAM from guide sequence in seqCol (default: False)
	  --PAM PAM             Only NGG PAM is currently supported (default: NGG)
	  --selectMethod {score,even}
	                        Method to select guides within elements. 'even' selects every Nth gRNA evenly across the region [recommended]. 'score' selects by best off-target score [not recommended]. (default: score)
	  --forceGuides FORCEGUIDES
	                        Pass a file containing one guide name per line to force selection of these gRNAs (useful for including gRNAs that we know work) (default: None)
	  --trimElements TRIMELEMENTS
	                        Remove guides corresponding to elements that have fewer than this many guides (e.g., because not enough guides existed to choose from, repeats, etc.). (default: 0)
	  --excludeRestrictionSites EXCLUDERESTRICTIONSITES
	                        Comma-delimited list of subsequences (e.g. restriction enzyme sites) to exclude gRNAs. (default: )
	  --barcodes BARCODES   File with one unique barcode sequence per line, e.g. for HyPR screens. Will substitute these sequences in place of '[NNNNN]' sequence from vector design file.  Errors out if there are more unique guides than barcodes. (default: None)


### Step 4: Combine multiple guide pools, add PCR handles, and output files to order oligo pools

First, set up a config file (see `example/CombinePoolsConfig.txt`) that lists paths to all of the guide subpools to process (and if desired, combine).

Then, run:

	python $CODEDIR/src/CombineGuidePools.py \
	    --config $PROJECT/04_CombinePools/CombinePoolsConfig.txt \
	    --outbase $PROJECT/04_CombinePools/210615_GATA1CRISPRiTiling  \
	    > $PROJECT/04_CombinePools/210615_GATA1CRISPRiTiling.out 

If you do not have multiple subpools that need specific subpool primers, leave those columns in the CombinePoolsConfig.txt file blank.

Documentation:

	usage: CombineGuidePools.py [-h] --config CONFIG --outbase OUTBASE
	                            [--fillToOligoPoolSize FILLTOOLIGOPOOLSIZE]
	                            [--includeReverseComplements]

	Combine different oligo subpools output by MakeGuidePool.py into a final oligo pool for ordering. 
	Terminology: "subpool" refers to a set of oligos with the same PCR handles on the outside, i.e. a set of oligos that will be amplified together
	Before running, set up a config file (e.g. called CombinePoolsConfig.txt) that lists the paths to each subpool, and the unique primers to use. 

	optional arguments:
	  -h, --help            show this help message and exit
	  --config CONFIG       Config file with columns: subpool pool DesignFile Multiply FwdPrimer RevPrimer (default: None)
	  --outbase OUTBASE
	  --fillToOligoPoolSize FILLTOOLIGOPOOLSIZE
	                        Total size of oligo pool; if total number of guides is less than this number, will output oligo order file with this total number of gRNAs.  If total guides is more than this number, will truncate. (default: -1)
	  --includeReverseComplements
	                        Include this flag to adding reverse complement oligos (e.g., as a hedge against strand-specific errors in synthesis). Suggest skipping this if this array includes tiling sequences (or HyPR barcodes) (default: False)

	Input config file:
	subpool                 Name of the subpool of gRNAs with unique primer handles to amplify 
	pool                    Arbitrary name of a set of gRNAs to include
	DesignFile              Path to design file for this pool (e.g., the *.design.txt file output by MakeGuidePool.py).  
	Multiply                Integer, default to 1. Use this column to balance the subpools by printing 
	                          the same sequences multiple times (so that all of the subpools are at 
	                          approximately the same abundance and require the same number of PCR cycles to amplify.)
	                          e.g. Multiply=2 will print all sequences for this subpool twice relative to subpools where Multiply=1
	FwdPrimer               Primer 1 for amplifying the subpool (order this sequence).  Leave blank if no subpool handles are required.
	RevPrimer               Primer 2 for amplifying the subpool (order this sequence).  Leave blank if no subpool handles are required.

	Key columns in output file:
	GuideSequenceMinusG     Sequence of the gRNA minus leading G - this should be completed by all guides, including negative controls
	CoreOligo               Sequence of gRNA + flanking promoter + scaffold sequences
	OligoSequence           Full oligo sequence to order, including subpool primers
	subpool                 Name of the subpool of gRNAs with unique primer handles to amplify 
	FwdPrimer               Primer 1 for amplifying the subpool (order this sequence)
	RevPrimer               Primer 2 for amplifying the subpool (order this sequence)
	guideSet                Gene / enhancer / target of this gRNA

	Current behavior is, if design file has identical oligo sequences, to collapse these before re-duplicating to fill the space.


## Usage (Prime editing)

A full example workflow for designing pegRNAs guides for a few control sites: `example/log.pegRNAsGATA1.sh`

### Step 1. Set up a table listing desired variants

Example:

	chr     start   end     name    ref     alt     region
	chr10   81107192        81107197        chr10:81107193: CCAAT>AGCCA     CCAAT   AGCCA   PPIF_CCAAT_box

Notes:
- This file must follow BED file format conventions: coordinates are zero-based for the coordinate start and one-based for the coordinate end.  So, the first 3 bases of chromosome 1 is `chr1	0	3`.
- The reference sequence listed in the 'ref' column must match the reference genome sequence provided (to check against coordinate position errors, e.g. see previous point).  


