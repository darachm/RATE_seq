#!/bin/bash

echo "
You're running the RNAseq scheduler script.
$(date)

...

great
"

# Module definitions
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


#Sequence format
isPE=1 #0 if SE
isUMI=0 #1 if UMIs used

#File locations
DIR='/scratch/cgsb/gencore/out/Gresham/2017-08-09_HN2WTAFXX/merged/'
hisat_ref=/genomics/genomes/Public/Fungi/Saccharomyces_cerevisiae/NCBI/R64/GCF_000146045.2_R64_genomic_hisat_index/genome
gff=/genomics/genomes/Public/Fungi/Saccharomyces_cerevisiae/NCBI/R64/GCF_000146045.2_R64_genomic.gff
#gtf=
#ref_flat=/scratch/cgsb/gresham/LABSHARE/Resources/Picard_files/ncbi_yeast_flatfile.txt
#ribosome_intervals=/scratch/cgsb/gresham/LABSHARE/Resources/Picard_files/ribosome_intervals_tabs
count_script=/scratch/cgsb/gresham/LABSHARE/Scripts/Pipelines/count_features.R

## list of input files
FILES=($(ls $DIR/*fastq.gz | perl -pe 's/^.+n0\d_(.+)\.fastq\.gz/$1/g' | sort | uniq))

## $ID will hold the sample_id
ID=${FILES[$SLURM_ARRAY_TASK_ID]}

## create directory for sample and change into it to organize output by sample                                                   
mkdir $ID
cd $ID

#############
# Alignment #
#############

if [ isPE == 1]

then
    
    ## full paths to paired mate files
    INPUT_1=$(ls $DIR/*n01_$ID.fastq.gz)
    INPUT_2=$(ls $DIR/*n02_$ID.fastq.gz)

    ##perform adaptor and quality trimming
    trim_galore --paired --fastqc $INPUT_1 $INPUT_2

    ##alignment with hisat2
    hisat2 -q -x ${hisat_ref} -1 HN2WTAFXX_n01_${ID}_val_1.fq.gz -2 HN2WTAFXX_n02_${ID}_val_2.fq.gz -S $ID.sam

    ##file renaming, sorting and indexing
    samtools sort -o ${ID}.bam ${ID}.sam
    samtools index ${ID}.bam

    ## PS - clean up intermediary .sam files
    rm ${ID}.sam

elif [ isPE == 0 ]

then

    ## full paths to single end files
    INPUT_1=$(ls $DIR/*n01_$ID.fastq.gz)

    ##perform adaptor and quality trimming
    trim_galore --fastqc $INPUT_1

    ##alignment with hisat2
    hisat2 -q -x ${hisat_ref} -1 HN2WTAFXX_n01_${ID}_val_1.fq.gz -S $ID.sam

    ##file renaming, sorting and indexing
    samtools sort -o ${ID}.bam ${ID}.sam
    samtools index ${ID}.bam

    #PS - clean up intermediairy .sam files
    rm ${ID}.sam

    if [ isUMI=1 ]
    
    then

	umi_tools dedup --umi-separator=: --method=directional --output-stats=${ID}.stats -I ${ID}.bam -S ${ID}.dedup.bam -L ${ID}.dedup.log

    fi

else
 
echo "Define analysis as either a paired end or single end sequencing protocol"

fi


##############
# QC metrics #
##############

#Run Qualimap on individual BAM file
#http://qualimap.bioinfo.cipf.es/doc_html/command_line.html#command-line
#NB. for qualimap rnaseq set -pe if paired end sequencing used 

if [isUMI == 0]

then
    qualimap bamqc -bam ${ID}.bam -gff ${GFF} -c
    qualimap rnaseq -bam ${ID}.bam -gtf ${GTF} -outdir rnaseq_qc_results

elif [isUMI == 1]

then

    qualimap bamqc -bam ${ID}.dedup.bam -gff ${GFF} -c
    qualimap rnaseq -bam ${ID}.bam -gtf ${GTF} -outdir rnaseq_qc_results
fi

####Steps below can only be run once all alignment and QC  steps are complete.

############################
# Generate table of counts #
############################

#Run from an R script
R --vanilla ${count_script}






