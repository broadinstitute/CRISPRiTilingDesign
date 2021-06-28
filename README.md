# CRISPRi Tiling Design - Engreitz Lab
Scripts for designing CRISPRi tiling screens.  

## Authors

* Jesse Engreitz (@engreitz)

## Description

This repository includes scripts used to design CRISPRi tiling screens for the Engreitz Lab.

The major steps in the design process are:  
1. Choosing regions (possibly overlapping) to tile with CRISPRi gRNAs  
2. Designing all possible gRNAs in those regions.  
3. Scoring and filtering gRNAs to eliminate off-target effects and poor efficacy gRNAs.  
4. Choosing gRNAs within each region, e.g. by even selection across the region.  
5. Collating different sets of gRNAs ("subpools") into a single oligo "pool" for ordering

The repository also includes helper scripts for designing dead gRNAs, selecting protein-coding 
gRNAs from Brie/Brunello, and annotating gRNAs from previous designers

Depends on CRISPR design code repository using the MIT specificity score, located [here](https://github.com/EngreitzLab/CRISPRDesigner/tree/master).

## Setup

### Step 1: Clone this github repository

[Clone](https://help.github.com/en/articles/cloning-a-repository) this to your local system, into the place where you want to perform the data analysis.

### Step 2: Install conda environment if needed

Install Snakemake and conda environment using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    conda env create --file envs/EngreitzLab.yml  

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

Engreitz Lab members working on Sherlock can skip this step.


## Usage

A full example workflow for designing CRISPRi guides tiling across all DNase peaks around GATA1 is provided in `example/log.sh`

### Step 1: Select regions in which to design gRNAs
