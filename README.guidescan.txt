To get the needed dependencies:

docker pull xerez/guidescan
docker tag xerez/guidescan guidescan
conda env create -f guidescan.yml

Put a copy of the guidecscan .bam and .abi in:
/mnt/DATA/CRISPR_design_scripts/guidescan_index

For the example, put 
- RefSeqCurated.200730.UniqueTSS.CollapsedGeneBounds.500bp.bed
- sizes (hg38)
in /mnt/DATA/CRISPR_design_scripts/hg38/.
