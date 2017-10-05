#!/bin/bash

## Log directory
LOCAL_DIR_LOGS="./tmp/logs"
if [ ! -d ${LOCAL_DIR_LOGS} ]; then mkdir -p ${LOCAL_DIR_LOGS}; fi
LOG="${LOCAL_DIR_LOGS}/captainslog_$( date +%Y%m%d%H%M%S ).txt"

# I like the weird header things, it makes sense. Does anyone know
#   how to program a vim macro to make these easily? Shouldn't be hard

##############################
# Preliminaries and niceties #
##############################

# What is the id of your experiment? Give me something short and
# without weird characters, non-descriptive is great (ie dme211).
# I'm going to use it a bunch in bash scripting.
EXPERIMENT_NAME="testing"

# What kind of sequencing are you trying to analyze?
#isPE=1  #0 if SE
#isUMI=0 #1 if UMIs used
# How many basepairs?
# readLength = 75
# ( could be useful for dynamic trimming, filtering, match scores )

# Module definitions - which ones are you using? 
# You probably don't want to touch these.
GFFREAD="gffread/intel/0.9.8c"
SAMTOOLS="samtools/intel/1.3.1"
HISAT2="hisat2/intel/2.0.5"
#trim_galore="trim_galore/0.4.4"
#cutadapt="cutadapt/intel/1.12"
#fastqc="fastqc/0.11.5"
#picard="picard/2.8.2"
#htseq="htseq/intel/0.6.1p1"
#pysam="pysam/intel/0.11.2.2"
#r="r/intel/3.3.2"
#qualimap="qualimap/2.2.1"
#rseqc="rseqc/intel/2.6.4"

###################
# Let's get going #
###################

echo "
You're running the RNAseq scheduler script on experiment 
\"${EXPERIMENT_NAME}\". I've made a log file at ${LOG} for you.

...

great
" | tee -a ${LOG}
# See that `tee` usage? That pipes the echo to ${LOG} and to stdout,
# so you'll see it on the command line and have it recorded as well.

#######################
# Moving files around #
#######################
# Then, we make local copies of everything. Why? 
# Reproducibility, via ASF1 ( Atomic Semantic Fusion 1 folder ).
#
# Local data directory, so in this folder you find the local copy
# of the data, configuration files, indicies, other crap.
LOCAL_DIR_DATA="./data"
if [ ! -d ${LOCAL_DIR_DATA} ] ; then mkdir -p ${LOCAL_DIR_DATA} ; fi
#
# For scratch files
#
LOCAL_DIR_TMP="./tmp"
if [ ! -d ${LOCAL_DIR_TMP} ] ; then mkdir -p ${LOCAL_DIR_TMP} ; fi
#
# Input data, so this should aim at the place on the gencore scratch
#
INPUT_DIR_FASTQS='/scratch/cgsb/gresham/LABSHARE/Scripts/Pipelines/PE_RNAseq'
LOCAL_DIR_FASTQS="${LOCAL_DIR_DATA}/$( basename ${INPUT_DIR_FASTQS})"
#cp -r ${INPUT_DIR_FASTQS} ${LOCAL_DIR_FASTQS}
#
# This should aim at the folder of HISAT2 indicies of your reference
#
REF_DIR_HISAT2INDEX="/genomics/genomes/Public/Fungi/Saccharomyces_cerevisiae/NCBI/R64/GCF_000146045.2_R64_genomic_hisat_index"
LOCAL_DIR_HISAT2INDEX="${LOCAL_DIR_DATA}/$( basename $REF_DIR_HISAT2INDEX)"
LOCAL_HISAT2INDEX="${LOCAL_DIR_HISAT2INDEX}/genome"
#cp -r ${REF_DIR_HISAT2INDEX} ${LOCAL_DIR_HISAT2INDEX}
#
# This should aim at an appropriate GFF for counting
#
REF_GFF="/genomics/genomes/Public/Fungi/Saccharomyces_cerevisiae/NCBI/R64/GCF_000146045.2_R64_genomic.gff"
LOCAL_GFF="${LOCAL_DIR_DATA}/$( basename ${REF_GFF})"
#cp -r ${REF_GFF} ${LOCAL_GFF}
#
# This is the name of the local GTF file. We make it in step 1
#
LOCAL_GTF=$( echo -n ${LOCAL_GFF} | sed 's/\.gff/.gtf/' )

echo "
Okay, I copied
  ${INPUT_DIR_FASTQS} to ${LOCAL_DIR_FASTQS}
,
  ${REF_DIR_HISAT2INDEX} to ${LOCAL_DIR_HISAT2INDEX}
,
  ${REF_GFF} to ${LOCAL_GFF}
, and I expect us to make a ${LOCAL_GTF}.
" | tee -a ${LOG}

##############################
# Generate gtf file from gff #
##############################
# You're going to see these below alot. It's how you submit a string
# as a sbatch job. This means you can schedule jobs independently
# for different steps, and can tweak the resource needs appropriately
# for each step. 
# The --wrap bit is the actual command. I've stuck ;'s in there 
# because I think the shell wraps the lines. \" are necessary to keep
# it quoted.
RUNSTRING_GFFREAD="sbatch 
  --mail-type=BEGIN,END,FAIL --mail-user=${USER}@nyu.edu
  --job-name=${EXPERIMENT_NAME}_GFFread
  --output=${LOCAL_DIR_LOGS}/%A_${EXPERIMENT_NAME}_GFFread.out 
  --error=${LOCAL_DIR_LOGS}/%A_${EXPERIMENT_NAME}_GFFread.err 
  --nodes=1 --ntasks-per-node=1 --mem=30GB --time=00:05:00
  --wrap=\"module purge;
    module load ${GFFREAD};
    gffread -E ${LOCAL_GFF} -T -o- > ${LOCAL_GTF};
    \"
  ;"
echo "Attempting to run this: 
  ${RUNSTRING_GFFREAD}
" | tee -a ${LOG}
#RESPONSE_GFFREAD=( $(eval ${RUNSTRING_GFFREAD} ) )
#JOBID_GFFREAD=${RESPONSE_GFFREAD[-1]} 
JOBID_GFFREAD="12"
echo "... and quoth the scheduler:
\" ${RESPONSE_GFFREAD[@]} \"
" | tee -a ${LOG}

#####################
# Define the inputs #
#####################

# This is interpreting the ls command as an array
LOCAL_FASTQS=( $( ls $LOCAL_DIR_FASTQS/*fastq.gz ) )
# When we spit the array out and use David's perl one-liner to
# cut out just the sample ids, then look for unique ones.
LOCAL_SAMPLEIDS=( $( echo ${LOCAL_FASTQS[@]} | \
  perl -pe 's/[^\s]+n0[0-9]_(.+?)\.fastq\.gz/$1/g' | sort | uniq 
  ) )

echo "Looking in ${LOCAL_DIR_FASTQS}, I see IDs:
  ${LOCAL_SAMPLEIDS[@]}

I'm going to start launching jobs for each of these.
" | tee -a ${LOG}

# Next, we iterate over each sample ID
for THIS_SAMPLE_ID in ${LOCAL_SAMPLEIDS[@]} ; do

  echo "Scheduling for ID ${THIS_SAMPLE_ID}
    " | tee -a ${LOG}

# Here, I'm only supporting reads that are 01 and 02.
# We find all the FASTQ files of that ID.
  THIS_SAMPLE_FILES=( $( \
    ls ${LOCAL_DIR_FASTQS}/*n0[12]_${THIS_SAMPLE_ID}.fastq.gz ) )
# And that let's us calculate if it's paired end or not.
  isPE=$(( ${#THIS_SAMPLE_FILES[@]} - 1 ))
  
  echo "I see files ${THIS_SAMPLE_FILES[@]} with that sample ID.
    " | tee -a ${LOG}
  echo -n "Since there are ${#THIS_SAMPLE_FILES[@]}, I conclude that
    the files are " | tee -a ${LOG}
  if [ "${isPE}" = 0 ] ; then
    echo "single ended reads.
      " | tee -a ${LOG}
  elif [ "${isPE}" = 1 ] ; then
    echo "paired end reads.
      " | tee -a ${LOG}
  else
    echo "... err I dunno. Giving up on this sample." | tee -a ${LOG}
    (>&2 echo "
    ERROR: Sample ${THIS_SAMPLE_ID} saw ${#THIS_SAMPLE_FILES[@]} 
    files, so gave up. Suggested action?
    " | tee -a ${LOG} )
    continue
  fi
# That error message above is in case we get something really weird.

  if [ "${isPE}" = 1 ] ; then
    THIS_SAMPLE_FASTQ_PAIR1=$( \
      ls ${LOCAL_DIR_FASTQS}/*01_${THIS_SAMPLE_ID}.fastq.gz )
    echo "Paired end file 1 is : ${THIS_SAMPLE_FASTQ_PAIR1}
      " | tee -a ${LOG}
    THIS_SAMPLE_FASTQ_PAIR2=$( \
      ls ${LOCAL_DIR_FASTQS}/*02_${THIS_SAMPLE_ID}.fastq.gz )
    echo "Paired end file 2 is : ${THIS_SAMPLE_FASTQ_PAIR2}
      " | tee -a ${LOG}
  elif [ "${isPE}" = 0 ] ; then
    THIS_SAMPLE_FASTQ=$( \
      ls ${LOCAL_DIR_FASTQS}/*01_${THIS_SAMPLE_ID}.fastq.gz )
    echo "Single end file is : ${THIS_SAMPLE_FASTQ}
      " | tee -a ${LOG}
  fi

#########################
# Alignment with HISAT2 #
#########################

# I'm aligning before trimming because this is a demo.

  RUNSTRING_HISAT2="sbatch 
    --mail-type=BEGIN,END,FAIL --mail-user=${USER}@nyu.edu
    --dependency=afterok:${JOBID_GFFREAD}
    --job-name=${EXPERIMENT_NAME}_HISAT2
    --output=${LOCAL_DIR_LOGS}/%A_${EXPERIMENT_NAME}_HISAT2.out 
    --error=${LOCAL_DIR_LOGS}/%A_${EXPERIMENT_NAME}_HISAT2.err 
    --nodes=1 --ntasks-per-node=8 --mem=30GB --time=00:30:00
    --wrap=\"module purge
      ;
      module load ${HISAT2}
      ;
      hisat2 -q -x ${LOCAL_HISAT2INDEX} 
        -1 ${THIS_SAMPLE_FASTQ_PAIR1}
        -2 ${THIS_SAMPLE_FASTQ_PAIR2}
        -S ${LOCAL_DIR_TMP}/hisat2aligned_${THIS_SAMPLE_ID}.sam 
      ;
      \"
    ;"
  echo "Attempting to run this: 
    ${RUNSTRING_HISAT2}
  " | tee -a ${LOG}
  RESPONSE_HISAT2=( $(eval ${RUNSTRING_HISAT2} ) )
  JOBID_HISAT2=${RESPONSE_HISAT2[-1]} 
  echo "... and quoth the scheduler:
  \" ${RESPONSE_HISAT2[@]} \"
  " | tee -a ${LOG}

  RUNSTRING_CLEANHISAT2="sbatch 
    --mail-type=BEGIN,END,FAIL --mail-user=${USER}@nyu.edu
    --dependency=afterok:${JOBID_HISAT2}
    --job-name=${EXPERIMENT_NAME}_CLEANHISAT2
    --output=${LOCAL_DIR_LOGS}/%A_${EXPERIMENT_NAME}_CLEANHISAT2.out 
    --error=${LOCAL_DIR_LOGS}/%A_${EXPERIMENT_NAME}_CLEANHISAT2.err 
    --nodes=1 --ntasks-per-node=1 --mem=30GB --time=00:05:00
    --wrap=\"module purge
      ;
      module load ${SAMTOOLS}
      ;
      samtools sort 
        ${LOCAL_DIR_TMP}/hisat2aligned_${THIS_SAMPLE_ID}.sam >
        ${LOCAL_DIR_TMP}/hisat2aligned_${THIS_SAMPLE_ID}.bam 
      ;
      samtools index 
        ${LOCAL_DIR_TMP}/hisat2aligned_${THIS_SAMPLE_ID}.bam 
      ;
      \"
    ;"
  echo "Attempting to run this: 
    ${RUNSTRING_CLEANHISAT2}
  " | tee -a ${LOG}
  RESPONSE_CLEANHISAT2=( $(eval ${RUNSTRING_CLEANHISAT2} ) )
  JOBID_CLEANHISAT2=${RESPONSE_CLEANHISAT2[-1]} 
  echo "... and quoth the scheduler:
  \" ${RESPONSE_CLEANHISAT2[@]} \"
  " | tee -a ${LOG}

  exit

done



exit

#sbatch --dependency=afterok:${RESPONSE[-1]} --mail-type=BEGIN,END,FAIL --mail-user=dhm267@       nyu.edu --output=${LOGS}/%A_${STEPNAME}_chaser.out --error=${LOGS}/%A_${STEPNAME}_chaser.err --    job-name=${STEPNAME}_chaser --time=00:01:00 --wrap="mv ${LOGS}/running_${STEPNAME} ${LOGS}/        done_${STEPNAME}"

#  --wrap=\"module purge; module load trim_galore;
#    trim_galore --paired --fastqc ${INPUT_1} ${INPUT_2}\"


## This is going to have to take custom adapter for using TrumiSeq
## Darach can write that part
#WORKDIR=$(pwd) 
#mkdir ${WORKDIR}"/tmp/dme211"
## data directory, obviously, pointing at directory with the fastqs demultiplexed
##   by the core

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


