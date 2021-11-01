#!/bin/bash


## give it a full dir name with CLIP bam files as argument

mydir=$1
myname=$(basename ${mydir})

echo "processing ${myname} from ${mydir}"

outdir="/scratch/AG_Ohler/svetlana/HeLa_CLIPs/omniCLIP_results_bulk/${myname}"
#outdir="/data/ohler/svetlana/pipelines/omniCLIP/example_data_result" ##testdata
#outdir="${mydir}/omniCLIPresults"

tmpdir="${outdir}/tmp"
mkdir -p $tmpdir

echo "output directory: ${outdir}"

#$-cwd #start from current directory
###$-l h_vmem=20G ## old ## now default is 4G ##h_vmem=m_mem_free*1,25!
#$ -l m_mem_free=30G
##$-l h_rt=48:0:0 #runtime ## cluster downtime otherwise comment out
#$-V #export all the environmental variables into the context of the job
#$-j yes #merge the stderr with the stdout
##$-o omni_bulk.log #stdout, job log
#$-m eas # send email beginning, end, and suspension
#$-M svetlana.lebedeva@mdc-berlin.de
##$-pe smp 1
#$-l data
#$-N 'omni_bulk'
##$ -l longrun


eval "$(/home/lebedeva/miniconda3/bin/conda shell.bash hook)" ## both workstation and cluster
conda activate omniCLIP
export PYTHONNOUSERSITE=True

### create db ###python data_parsing/CreateGeneAnnotDB.py INPUT.gff OUTPUT.gff.db ## use new gffutils in a separate environment
myGenome="/data/ohler/svetlana/genomes/hg19"
myAnno="gencode.v19.annotation.gff3.new.db"

#input_rmd="/data/ohler/svetlana/Stoecklin/RNA_seq_Fabian/C8Y28ACXX_FP_May16_16s002853-1-2_Ibberson_lane1R1_sequence.aligned_plus.bam"
input_dmso="/data/ohler/svetlana/Stoecklin/RNA_seq_Fabian/C8Y28ACXX_FP_May16_16s002853-1-2_Ibberson_lane1C1_sequence.aligned_plus.bam"

### list all bam files and assign them to vars https://unix.stackexchange.com/a/197183

clipinput=""
#for file in $( ls ${mydir}/PUM2*_chr1.bam ); do ##test data
for file in $( ls ${mydir}/*_uniq.bam ); do
    clipinput=$(echo "$clipinput --clip-files $file")
done

eval "python omniCLIP.py ${clipinput} --annot  ${myAnno} --genome-dir ${myGenome}  --save-tmp --bg-files ${input_dmso} --out-dir ${outdir} --tmp-dir ${tmpdir} --filter-snps --max-it 10 --bck-var --nb-cores 1 --nr_mix_comp 10 --diag_event_mod DirchMultK --diag-bg --norm_class --trans-model binary_bck --bg-type Coverage_bck --mask-miRNA --rev_strand 0 --use-precomp-CLIP-data --use-precomp-bg-data > ${outdir}/${myname}.log 2>&1"

### add this if restarting failed run
###--use-precomp-CLIP-data --use-precomp-bg-data --restart-from-iter


###test data
#eval "python omniCLIP.py $clipinput --annot example_data/gencode.v19.annotation.chr1.gtf.db --genome-dir example_data/hg37/ --bg-files example_data/RZ_rep1_chr1.bam --bg-files example_data/RZ_rep2_chr1.bam --out-dir ${outdir} --collapsed-CLIP --bck-var"

exit 0
