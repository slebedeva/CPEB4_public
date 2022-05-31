# CPEB4_public

Code for CLIP data analysis for the CPEB4 manuscript.


## If you would like to reproduce figures from the processed data

0. Look into config.yml. It should work with default values. However if you have, for example, the hg19 genome you can change the path to it to shorten the runtime and to avoid downloading it again. You can also specify custom paths to the plots and results folders. Default paths are relative to the directory where manuscript_figures.R is.

1. Get the R container and install relevant packages (if you do not have R installation).

Presuming that you have no root access and are on linux:

A. [Install miniconda3](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html): follow the instructions for your system.

You will need to [download](https://docs.conda.io/en/latest/miniconda.html#linux-installers) and run a bash installation script.

`bash Miniconda3-latest-Linux-x86_64.sh`

B. Install singularity in a conda environment

```
eval "$($HOME/miniconda3/bin/conda shell.bash hook)"
conda create --name singularity singularity=3.7.1
conda activate singularity
```

C. Pull the latest R bioconductor container

`singularity pull bioc.sif docker://bioconductor/bioconductor_docker:devel`

D. Install the necessary packages 

singularity exec bioc.sif Rscript install_packages.R

E. Run the script to produce figures 

(It requires memory and time to run. If it crashes, just run it again until it finishes).

`singularity exec bioc.sif Rscript manuscript_figures.R`

## Explanations for directories:

### Directory data:

Contains the data encessary to run main R scripts, including bam and bed files of CPEB4 CLIP that result from the mapping pipeline and omniCLIP.

### Directory scripts:

Auxiliary R scripts (should be called from within manuscript_figures.R).

### Main directory: 

R code to to create figures: manuscript_figures.R. 

The script should produce a "plots" folder with pdf images of raw figure panels as well as "results" folder with the R objects and "source_data" subfolder with the csv files that were used for plots.

---

## If you want to look into analysis of the raw CLIP sequencing data

### Directory mapping_pipeline

Contains the pipelines used for mapping 1) p32 and 2) ir CPEB4 CLIP data as well as conda environment file with software version specifications.

### Directory omniCLIP_scripts

Contains the scripts used to run omniCLIP on CPEB4 data as well as multiple published HeLa CLIP data. The omniCLIP version used for this, as well as scripts to generate the database are in this [omniCLIP fork](https://github.com/slebedeva/omniCLIP). The general pipeline to map any CLIP data is in [this repository](https://github.com/slebedeva/CLIP_mapping).
