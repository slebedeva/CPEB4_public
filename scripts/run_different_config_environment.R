## test other DEs (thick=sumDE)

Sys.setenv(R_CONFIG_ACTIVE = "new_DE")

if(!config::is_active("new_DE")){message:"make sure your environment is correct!"}

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
if(!file.exists(file.path(resultdir,"bed_reduced.RData"))){source(file.path(scriptdir,"process_CLIP_peaks.R"))}
if(!file.exists(file.path(resultdir,"DEseq.RData"))){source(file.path(scriptdir,"Differential_expression_analysis.R"))}


## load old data and do a different "thick"

ctrl2$thick <- IRanges(start=ifelse(!is.na(ctrl2$thickDESum),ctrl2$thickDESum,
                                    (start(ranges(ctrl2))+round(width(ranges(ctrl2))/2))),
                       width=1)

rmd2$thick <- IRanges(start=ifelse(!is.na(rmd2$thickDESum),rmd2$thickDESum,
                                    (start(ranges(rmd2))+round(width(ranges(rmd2))/2))),
                       width=1)


export.bed(ctrl2, file.path(resultdir, "CPEB4_Ctrl_hg19.bed"))
export.bed(rmd2, file.path(resultdir, "CPEB4_RMD_hg19.bed")) 

save(rmd, ctrl,rmd2,ctrl2,file= file.path(resultdir, "bed_reduced.RData")) 
