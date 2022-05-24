# CPEB4_public

Code for CLIP data analysis for the CPEB4 manuscript.

- If you would like to reproduce the figures from the paper: start with manuscript_figures.R.

- If you would like to start with raw fastq files: look into mapping_pipeline and omniCLIP_scripts.

## Reproducing figures

1. Look into config.yml. It should work as it is already. If you have, for example, the hg19 genome you can change the path to it to shorten the runtime and to avoid downloading it again. You can also specify paths to the results folders. Default paths are relative to the directory where manuscript_figures.R is.

2. Run manuscript_figures.R. Give it enough memory (10G). Especially the genome browser screenshot figures require memory and time to run. 

### Directory data:

Contains the data encessary to run main R scripts, including bam and bed files of CPEB4 CLIP that result from the mapping pipeline and omniCLIP.

### Directory scripts:

Auxiliary R scripts (should be called from within manuscript_figures.R).

### Main directory: R code to to create figures.

Start with manuscript_figures.R. Make sure that you have Bioconductor installed. 

The script will try to install R packages that are missing (but may fail due to dependencies).

The script should produce a "plots" folder with pdf images of raw figure panels as well as "results" folder with the R objects and "source_data" subfolder with the csv files that were used for plots.

---

## Restarting from raw data

### Directory mapping_pipeline

Contains the pipelines used for mapping 1) p32 and 2) ir CPEB4 CLIP data as well as conda environment file with software version specifications.

### Directory omniCLIP_scripts

Contains the scripts used to run omniCLIP on CPEB4 data as well as multiple published HeLa CLIP data. The omniCLIP version used for this, as well as scripts to generate the database are in this [omniCLIP fork](https://github.com/slebedeva/omniCLIP). The general pipeline to map any CLIP data is in [this repository](https://github.com/slebedeva/CLIP_mapping).
