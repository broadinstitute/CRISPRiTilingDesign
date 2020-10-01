To get the needed dependencies:

docker pull quay.io/biocontainers/guidescan:1.2--0
docker tag  quay.io/biocontainers/guidescan:1.2--0 guidescan

conda env create -f guidescan.yml

Put a copy of the guidecscan .bam and .abi in:
/mnt/DATA/CRISPR_design_scripts/guidescan_index

For the example (all TSSes) these files are needed in /mnt/DATA/CRISPR_design_scripts/hg38/:
  - RefSeqCurated.200730.UniqueTSS.CollapsedGeneBounds.500bp.bed
  - sizes (hg38)

