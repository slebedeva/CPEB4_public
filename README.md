# CPEB4_public

Code for CLIP data analysis for the CPEB4 manuscript.

### Directory mapping_pipeline

Contains the pipelines used for mapping 1) p32 and 2) ir CPEB4 CLIP data as well as conda environment file with software version specificaitons.

### Directory omniCLIP_scripts

Contains the scripts used to run omniCLIP on CPEB4 data as well as multiple published HeLa CLIP data. The omniCLIP version used for this, as well as scripts to generate the database are in this [omniCLIP fork](https://github.com/slebedeva/omniCLIP). The general pipeline to map any CLIP data is in [this repository](https://github.com/slebedeva/CLIP_mapping).

### Directories data and annotaiton:

Contain data encessary to run main R scripts, including bam and bed files of CPEB4 CLIP that result from the mapping pipeline and omniCLIP.

### Main directory: R code to to create figures.

Start with manuscript_figures.R. Make sure that you have Bioconductor installed. 

The script will try to install R packages that are missing (but may fail due to dependencies).

The script should produce a "plots" folder with pdf images of raw figure panels as well as "source_data" folder with the csv files that were used for plots.
