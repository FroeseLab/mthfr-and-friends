# mthfr-and-friends
[![Snakemake](https://img.shields.io/badge/snakemake-≥8.11.3-brightgreen.svg)](https://snakemake.github.io)

A repository for the ap ms analysis of MTHFR interactions for the 
paper: "Evidence for Interaction of 5,10-Methylenetetrahydrofolate Reductase (MTHFR) with Methylenetetrahydrofolate Dehydrogenase (MTHFD1) and General Control Nonderepressible 1 (GCN1)"
available as preprint on: 
https://www.biorxiv.org/content/10.1101/2024.08.22.609157v1

## How to run the analysis
### Development
This analysis requires `conda` to run.

Install it eg by following the `miniforge` tutorial: https://github.com/conda-forge/miniforge

#### Preparation
Install the `mthfr-and-friends` conda environment using:

```
mamba env create -n mthfr-and-friends -f environment.yml
```
Activate the environment using

```
conda activate mthfr-and-friends
```

### Run the analysis

To run the complete analysis activate the conda environment and run

```
snakemake -c4 --sdm conda
```

To create the anlysis report run:

```
snakemake -c4 --sdm conda --report report.html
```
