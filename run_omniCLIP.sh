#!/bin/bash

#$-cwd #start from current directory
#$ -l m_mem_free=20G
##$-l h_rt=36:0:0 #runtime ## cluster downtime otherwise comment out
#$-V #export all the environmental variables into the context of the job
#$-j yes #merge the stderr with the stdout
#$-o logs #stdout, job log
#$-m eas # send email beginning, end, and suspension
#$-M svetlana.lebedeva@mdc-berlin.de
#$-pe smp 4 ## crushes and abuses the cpu (not sure why)
#$-l data
#$-N 'omni3'
##$ -l longrun #if needs more than 4 days something is wrong


########## !!! change these values every time: #############
## 1. data dir
## 20200907 new irCLIP (qc fail sample)
#mydir="/fast/AG_Ohler/svetlana/Stoecklin/irCLIP_20200904"
## 20201014 resequence with input 
mydir="/fast/AG_Ohler/svetlana/Stoecklin/irCLIP_20201008/pipeline_clip_mapping_hg19" ##hg19
#mydir="/fast/AG_Ohler/svetlana/Stoecklin/irCLIP_20201008" ##hg38 ##currently broken - OK change chr1 1  to chr1 in fasta files
## need bam files form irCLIP and p32 clip together
mydir2="/fast/AG_Ohler/svetlana/Stoecklin/p32_CLIP"

## 2. sample name
myname="Ctrl"
#myname="RMD"

## 3. where to output 
#outdir="${mydir}/omniCLIP_results/input/${myname}" ##input
#outdir="${mydir}/omniCLIP_results/CLIP/${myname}" ##CLIP
outdir="/fast/AG_Ohler/svetlana/Stoecklin/omniCLIP_results/p32_ir_together/${myname}"
tmpdir="${outdir}/tmp"
mkdir -p $tmpdir

## 4. database to use (hg19 or hg38)
### how to create db see make_omni3db script, use gff3 from gencode #python data_parsing/CreateGeneAnnotDB.py INPUT.gff OUTPUT.gff.db
## hg19
myGenome="/data/ohler/svetlana/genomes/hg19"
myAnno="gencode.v19.annotation.gff3.new.db"
## hg38
#myGenome="/fast/AG_Ohler/Genomes/human/hg38/genomes/GRCh38.p13.genome.fa.split"
#myAnno="/fast/AG_Ohler/Genomes/human/hg38/annotation/gencode.v35.annotation.gff3.db"

## 5. Background bam files ## change accordingly: RNAseq as background samples (use input as another fake CLIP and substract these peaks as false positives)
#inputdir="/fast/AG_Ohler/svetlana/Stoecklin/RNA_seq_Fabian/pipeline_star_single_hg38/star/" ##hg38
inputdir="/fast/AG_Ohler/svetlana/Stoecklin/RNA_seq_Fabian/star/" ##hg19
input_rmd=""
input_ctrl=""
for file in $( ls ${inputdir}/*lane1R*plus ); do
    input_rmd=$(echo "$input_rmd --bg-files $file")
done
for file in $( ls ${inputdir}/*lane1C*plus ); do
    input_ctrl=$(echo "$input_ctrl --bg-files $file")
done

## 5a. which background to use
backgr=$input_ctrl
#backgr=$input_rmd

## 6. Actual CLIP bam files - change CLIP or input!!! input should be ".plus" as well!
### list all bam files and assign them to vars https://unix.stackexchange.com/a/197183
clipinput=""
#for file in $( ls ${mydir}/merged/${myname}*_input.merged_uniq.bam.plus ); do clipinput=$(echo "$clipinput --clip-files $file"); done
for file in $( ls ${mydir}/merged/${myname}*_CLIP.merged_uniq.bam ); do clipinput=$(echo "$clipinput --clip-files $file"); done
clipinput=$(echo "$clipinput --clip-files ${mydir2}/merged/CPEB4_CLIP_r1_CLIP_4SU_DMSO.merged_uniq.bam") ##ctrl

echo "$(date) processing ${myname} from ${mydir}..nd ${mydir2}. CLIP is ${clipinput}, background is ${backgr}, output is in ${outdir}"

eval "$(/home/lebedeva/miniconda3/bin/conda shell.bash hook)" ## both workstation and cluster
conda activate omniCLIP3
export PYTHONNOUSERSITE=True


eval "python omniCLIP.py $clipinput --annot  $myAnno --genome-dir $myGenome  --save-tmp $backgr --out-dir $outdir --tmp-dir $tmpdir --filter-snps --max-it 10 --bck-var --nb-cores 4 --nr_mix_comp 5 --diag_event_mod DirchMultK --diag-bg --norm_class --bg-type Coverage_bck --mask-miRNA --rev_strand 0 --seed 42 --verbosity 1"
## if restarting failed run, add #--use-precomp-CLIP-data --use-precomp-bg-data --restart-from-iter
## params from Philipp were  --nr_mix_comp 10 and no --max-mismatch (default is 2) and no seed and pv

exit 0
