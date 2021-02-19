First, pull in the dependencies with the following commands:

- install docker
- pull image
  docker pull quay.io/biocontainers/guidescan:1.2--0
  docker tag  quay.io/biocontainers/guidescan:1.2--0 guidescan

- install conda
- set up conda env
  conda env create -f guidescan.yml

- fetch guidescan indices from http://www.guidescan.com/help#Databases
  copy cas9_hg38_all_guides.bam and cas9_hg38_all_guides.bai to
  a directory, such as /mnt/DATA/CRISPR_design_scripts/guidescan_index

After these are in place, you may need to edit paths in
example/design_guides.sh

Run with:
design_guides.sh regions.bed projectdir

For example:
mkdir -p ./all_TSSes
design_guides.sh RefSeqCurated.200730.UniqueTSS.CollapsedGeneBounds.500bp.bed ./all_TSSes/
