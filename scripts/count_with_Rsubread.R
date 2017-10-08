#generate counts using Rsubread in R

args <- commandArgs(trailingOnly=TRUE)
# trailingOnly means after "--args"

message(paste(args))

if ( length(args) != 8 ) stop("ZERP. Gimme 7 args. Read the script.")
# arg 1 is the script or call, I think

#bam file to be counted
bam_file <- args[2]
#output filename
output_file <- args[3]
#annotation file
annotation.file <- args[4]
# Are we not strand specific?
strandSpecific <- as.logical(args[5])
# Are we not paired end?
pairedEnd <- as.logical(args[6])
# What feature are we counting?
gtfFeatureType <- args[7]
# What attribute type are we calling it?
gtfAttributeType <- args[8]

message("Reading in ",bam_file)

#source("http://bioconductor.org/biocLite.R")
#biocLite("Rsubread")

library("Rsubread")
#use featureCounts to count reads to count transcripts per gene
counts <- featureCounts(files=bam_file
  ,annot.ext=annotation.file, useMetaFeatures=T
  ,isGTFAnnotationFile=T
  ,GTF.featureType=gtfFeatureType, GTF.attrType=gtfAttributeType
  ,strandSpecific=strandSpecific, isPairedEnd=pairedEnd
  )

message("Spitting out counts in ",output_file)

#write counts out in a table
write.table( data.frame( counts$counts, stringsAsFactors=FALSE)
  ,file=output_file ,quote=FALSE, sep="\t"
  )
