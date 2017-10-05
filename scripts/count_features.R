#generate counts using Rsubread in R

#annotation file
annotation.file <- /genomics/genomes/Public/Fungi/Saccharomyces_cerevisiae/NCBI/R64/GCF_000146045.2_R64_genomic.gff
#bam files to be analyzed
bam_files <- list.files(pattern = ".bam")
#outputfile name
output_file <- "counts.txt"

#source("http://bioconductor.org/biocLite.R")
#biocLite("Rsubread")
library(Rsubread)

#use featureCounts to count reads to count transcripts per gene
counts<-featureCounts(files=bam_files,annot.ext=annotation.file,useMetaFeatures=T,isGTFAnnotationFile=T,GTF.featureType="exon",GTF.attrType="gene_id",strandSpecific=0,isPairedEnd=T)

#write counts out in a table
write.table(data.frame(counts$counts, stringsAsFactors=FALSE), file=output_file, quote=FALSE, sep="\t")
