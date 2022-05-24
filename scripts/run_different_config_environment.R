## to run different environment

Sys.setenv(R_CONFIG_ACTIVE = "old_DE")

if(!config::is_active("old_DE")){message:"make sure your environment is correct!"}

config=config::get()

## make sure those directories exist
plotdir=config$plotdir ##Plots
datadir=config$datadir ##data
ann_dir=config$ann_dir ##annotation
scriptdir=config$scriptdir ##scripts
resultdir=config$resultdir ## results
for(mydir in c(ann_dir,plotdir,resultdir)){
  if(!dir.exists(mydir)) {dir.create(mydir)}
}

## folder for R scripts
## warning if cannot find scripts
if(! dir.exists(scriptdir)){message("cannot find scripts, please check your scripts directory path in the config file!")}


############### generate data ###################

# important: they have to be executed in this order
if(!file.exists(file.path(ann_dir,"myregions.RData"))){source(file.path(scriptdir,"generate_annotation.R"))}
#if(!file.exists(file.path(resultdir,"bed_reduced.RData"))){source(file.path(scriptdir,"process_CLIP_peaks.R"))}
#if(!file.exists(file.path(resultdir,"DEseq.RData"))){source(file.path(scriptdir,"Differential_expression_analysis.R"))}

## change what you want to change here

