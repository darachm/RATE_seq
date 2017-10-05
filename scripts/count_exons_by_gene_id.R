#generate counts using Rsubread in R

args <- commandArgs(trailingOnly=TRUE)

if ( length(args) != 3 ) stop("ZERP. Gimme 3 args. Read the script.")

#bam file to be counted
bam_file <- args[1]

#output filename
output_file <- args[2]

#annotation file
annotation.file <- args[3]

strandSpecific <- as.logical(args[4])

pairedEnd <- as.logical(args[5])


#source("http://bioconductor.org/biocLite.R")
#biocLite("Rsubread")

#use featureCounts to count reads to count transcripts per gene
counts <- Rsubread::featureCounts( files=bam_file
  ,annot.ext=annotation.file, useMetaFeatures=T
  ,isGTFAnnotationFile=T
  ,GTF.featureType="exon", GTF.attrType="gene_id"
  ,strandSpecific=strandSpecific, isPairedEnd=pairedEnd
  )

#write counts out in a table
write.table( data.frame( counts$counts, stringsAsFactors=FALSE)
  ,file=output_file ,quote=FALSE, sep="\t"
  )
