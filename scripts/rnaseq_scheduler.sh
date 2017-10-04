#!/bin/bash

echo "
You're running the RNAseq scheduler script.
$(date)

...

great
"

### Preliminaries and niceties

# Your email string
emailString="dhm267@nyu.edu"

# What is the id of your experiment? Give me something short and
# without weird characters, non-descriptive is great (ie dme211).
experimentName="testing"

# What kind of sequencing are you trying to analyze?
isPE=1  #0 if SE
isUMI=0 #1 if UMIs used

# How many basepairs?
# readLength = 75
# ( could be useful for dynamic trimming, filtering, match scores )

### Module definitions
trim_galore="trim_galore/0.4.4"
cutadapt="cutadapt/intel/1.12"
fastqc="fastqc/0.11.5"
hisat2="hisat2/intel/2.0.5"
samtools="samtools/intel/1.3.1"
picard="picard/2.8.2"
htseq="htseq/intel/0.6.1p1"
pysam="pysam/intel/0.11.2.2"
r="r/intel/3.3.2"
qualimap="qualimap/2.2.1"
rseqc="rseqc/intel/2.6.4"

### Paths
# Input data, where the fastqs you're trying to align
INPUT_DATA_DIR='/scratch/cgsb/gencore/out/Gresham/2017-08-09_HN2WTAFXX/merged'
# HISAT2 index file
HISAT_INDEX='/genomics/genomes/Public/Fungi/Saccharomyces_cerevisiae/NCBI/R64/GCF_000146045.2_R64_genomic_hisat_index/genome'
# Reference GFF file, for counting and things
REF_GFF='/genomics/genomes/Public/Fungi/Saccharomyces_cerevisiae/NCBI/R64/GCF_000146045.2_R64_genomic.gff'
#gtf=
#ref_flat=/scratch/cgsb/gresham/LABSHARE/Resources/Picard_files/ncbi_yeast_flatfile.txt
#ribosome_intervals=/scratch/cgsb/gresham/LABSHARE/Resources/Picard_files/ribosome_intervals_tabs
# Log directory
LOGS="./logs"
if [ ! -d ${LOGS} ]; then mkdir ${LOGS}; fi

# First, identify the input files and copy them locally for 
# reproducibility.
# ... we'll do that later
# For now, we just list and do this direct
# List of input files
INPUT_FILES=($( ls $INPUT_DATA_DIR/*fastq.gz ))
ID_LIST=( $( echo ${INPUT_FILES[@]} | \
  perl -pe 's/[^\s]+n0[0-2]_(.+?)\.fastq\.gz/$1/g' | sort | uniq 
  ) )

echo "I see IDs:
  ${ID_LIST[@]}"
echo

for THIS_SAMPLE_ID in ${ID_LIST[@]} ; do

  echo "Scheduling for ID ${THIS_SAMPLE_ID}"

#  INPUT_FILE_LIST=$( ls $INPUT_DATA_DIR/*n0[0-2]_${THIS_SAMPLE_ID}*.fastq.gz )

#  echo "  I see ${INPUT_FILE_LIST[@]}"
#  echo

#  if [ isPE == 1]; then
    ## full paths to paired mate files
#    INPUT_1=$(ls $DIR/*n01_$ID.fastq.gz)
#    INPUT_2=$(ls $DIR/*n02_$ID.fastq.gz)
#
#    ##perform adaptor and quality trimming

done





#
#TRIM_RUNSTRING="sbatch 
#  --mail-type=BEGIN,END,FAIL --mail-user=${emailString} 
#  --job-name=${experimentName}_Trimming
#  --output=${LOGS}/%A_${experimentName}_Trimming.out 
#  --error=${LOGS}/%A_${experimentName}_Trimming.err 
#  --nodes=1 --ntasks-per-node=1 --mem=30GB --time=00:05:00
#  --wrap=\"module purge; module load trim_galore;
#    trim_galore --paired --fastqc ${INPUT_1} ${INPUT_2}\"
#  ;"
#echo "    Running this:
#${RUNSTRING}"
#RESPONSE=$(eval ${RUNSTRING} )                                      
#RESPONSE=($RESPONSE)                                                
#sbatch --dependency=afterok:${RESPONSE[-1]} --mail-type=BEGIN,END,FAIL --mail-user=dhm267@       nyu.edu --output=${LOGS}/%A_${STEPNAME}_chaser.out --error=${LOGS}/%A_${STEPNAME}_chaser.err --    job-name=${STEPNAME}_chaser --time=00:01:00 --wrap="mv ${LOGS}/running_${STEPNAME} ${LOGS}/        done_${STEPNAME}"
#
#
## This is going to have to take custom adapter for using TrumiSeq
## Darach can write that part
#WORKDIR=$(pwd) 
#mkdir ${WORKDIR}"/tmp/dme211"
## data directory, obviously, pointing at directory with the fastqs demultiplexed
##   by the core
#DATADIR="/data/cgsb/gencore/out/Gresham/2017-01-06_HGGNWBGX2/new/"
#
#module purge
#module load cutadapt/intel/1.12
#
## the below reads in an index, parses it into bash arrays, humor me okay?
#
#unset indicies
#declare -A indicies
#unset adapterName
#declare -A adapterName
#IFS=$'\n';
#for i in $(tail -n +2 data/dme211/dme211barcodeIndex.csv );do 
#  thisSampleName=$(echo -n $i | perl -ne '/^(.+?),(.+?),(.+?),(.+?)$/;print $1;' ); 
#  thisAdapterName=$(echo -n $i | perl -ne '/^(.+?),(.+?),(.+?),(.+?)$/;print $3;' ); 
#  adapterIndex=$(echo -n $i | perl -ne '/^(.+?),(.+?),(.+?),(.+?)$/;print $4;' ); 
#  indicies["${thisAdapterName}"]="${adapterIndex}";
#  adapterName["${thisSampleName}"]="${thisAdapterName}";
#done;
#echo "Read in:"
#echo ${!adapterName[@]}
#echo "as mapping to"
#echo ${adapterName[@]}
#echo "and"
#echo ${!indicies[@]}
#echo "as mapping to"
#echo "${indicies[@]}"
#echo 
#
#unset adaptSeq
#declare -A adaptSeq
#IFS=$'\n';
#for i in $(tail -n +2 data/dme211/trumiseqAdapters.csv );do 
#  thisSampleName=$(echo -n $i | perl -ne '/^(DG.+)_P7,(.+?)$/;print $1;' ); 
#  if [ "$thisSampleName" != "" ]; then
#    thisAdapterSeq=$(echo -n $i | perl -ne '/^(DGseq_.+?)_P7,(.+?)$/;print $2;' ); 
#    adaptSeq["$thisSampleName"]="${thisAdapterSeq}";
#  fi
#done;
#echo "Read in:"
#echo ${!adaptSeq[@]}
#echo "as mapping to"
#echo ${adaptSeq[@]}
#echo 
#
#for i in $(/bin/ls $DATADIR | grep "_[wq]1\?[0-9]_"); do
#  echo `date`
#  thisSampleName=$(echo -n $i | gawk -F _ '{print $3}');
#  thisAdapterName=${adapterName["$thisSampleName"]}
#  echo "doing file $i, which is $thisSampleName, which is $thisAdapterName"
#  thisAdaptSeq=${adaptSeq["$thisAdapterName"]}
#  thisIndex=${indicies["${adapterName["$thisSampleName"]}"]}
## below we grab the non-empty lines from the fastq, as there was a
## bioinformatics hiccup previously
#  runstring="
#cat ${DATADIR}$i | grep --regexp ^$ --invert-match > tmp/dme211/fixed$i;
#cutadapt -a A${thisAdaptSeq} --cut=0 -o tmp/dme211/dme211.${thisSampleName}.${adapterName["$thisSampleName"]}.$thisIndex.adapterTrimmed.fastq tmp/dme211/fixed$i"
#
#
##########################
## Alignment with HISAT2 #
##########################
#
#    ##alignment with hisat2
#    hisat2 -q -x ${hisat_ref} \
#      -1 HN2WTAFXX_n01_${ID}_val_1.fq.gz \
#      -2 HN2WTAFXX_n02_${ID}_val_2.fq.gz \
#      -S $ID.sam
#
#    ##file renaming, sorting and indexing
#    samtools sort -o ${ID}.bam ${ID}.sam
#    samtools index ${ID}.bam
#
#    ## PS - clean up intermediary .sam files
#    rm ${ID}.sam

#elif [ isPE == 0 ]
#
#then
#
#    ## full paths to single end files
#    INPUT_1=$(ls $DIR/*n01_$ID.fastq.gz)
#
#    ##perform adaptor and quality trimming
#    trim_galore --fastqc $INPUT_1
#
#    ##alignment with hisat2
#    hisat2 -q -x ${hisat_ref} -1 HN2WTAFXX_n01_${ID}_val_1.fq.gz -S $ID.sam
#
#    ##file renaming, sorting and indexing
#    samtools sort -o ${ID}.bam ${ID}.sam
#    samtools index ${ID}.bam
#
#    #PS - clean up intermediairy .sam files
#    rm ${ID}.sam
#
#    if [ isUMI=1 ]
#    
#    then
#
#	umi_tools dedup --umi-separator=: --method=directional --output-stats=${ID}.stats -I ${ID}.bam -S ${ID}.dedup.bam -L ${ID}.dedup.log
#
#    fi
#
#else
# 
#echo "Define analysis as either a paired end or single end sequencing protocol"
#
#fi


##############
# QC metrics #
##############

#Run Qualimap on individual BAM file
#http://qualimap.bioinfo.cipf.es/doc_html/command_line.html#command-line
#NB. for qualimap rnaseq set -pe if paired end sequencing used 

#if [isUMI == 0]
#
#then
#    qualimap bamqc -bam ${ID}.bam -gff ${GFF} -c
#    qualimap rnaseq -bam ${ID}.bam -gtf ${GTF} -outdir rnaseq_qc_results
#
#elif [isUMI == 1]
#
#then
#
#    qualimap bamqc -bam ${ID}.dedup.bam -gff ${GFF} -c
#    qualimap rnaseq -bam ${ID}.bam -gtf ${GTF} -outdir rnaseq_qc_results
#
#fi
#
#####Steps below can only be run once all alignment and QC  steps are complete.
#
#############################
## Generate table of counts #
#############################
#
## make this a local copy
#Rscript /scratch/cgsb/gresham/LABSHARE/Scripts/Pipelines/count_features.R






# This archived for 
#SBATCH --array=0-23
#SBATCH --output=logs/%A_%a.o
#SBATCH --error=logs/%A_%a.e
#SBATCH --job-name=%a
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --mail-user=${USER}@nyu.edu
#SBATCH --mail-type=FAIL
