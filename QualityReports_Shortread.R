
library(ShortRead) # Bioconductor Package
library(Cairo)


#########################
#   Read FASTQ Files    #
#########################

# set the directory and read the list of fastq.gz files
setwd("") # file path to fastq.gz files
fastqDir = getwd()
fastqPath = list.files(fastqDir, pattern = ".fastq.gz$", full = TRUE)


#########################
#   Quality Report     #
#########################

#depending on the number of files this could take a while
qaSummary = qa(fastqPath, type ="fastq")   # qa() runs the quality assesment
browseURL(report(qaSummary))               #opens the report in a URL



#########################
#  Read Summary Export  #
#########################


# over represented reads
freqreads = ShortRead:::.freqSequences(qaSummary, "read")
write.csv(as.data.frame(freqreads), file = "shortread_frequentreads.csv")

# reads across samples as data table
reads = ShortRead:::.ppnCount(qaSummary[["readCounts"]])
write.csv(as.data.frame(reads), file = "shortread_readcounts.csv")

#look for possible adaptor contamination
contam = ShortRead:::.ppnCount(qaSummary[["adapterContamination"]]) 
write.csv(as.data.frame(contam), file = "shortread_adapterContamination.csv")

#########################
#   Individual Plots    #
#########################


df_quality <- qaSummary[["readQualityScore"]]
df_dist <- qaSummary[["sequenceDistribution"]]
perCycle <- qaSummary[["perCycle"]]

#save the individual plots as high resolution images
CairoPNG(filename = "Shortread%03d.png", 
         units = "in",
         height = 8,
         width = 8, 
         pointsize = 12, 
         dpi = 300,
         bg = "white")
ShortRead:::.plotReadCount(qaSummary)  # read counts across samples
ShortRead:::.plotNucleotideCount(qaSummary)  # nucleotide distribution across samples
ShortRead:::.plotReadQuality(df_quality[df_quality$type=="read",])  # average quality scores summary
ShortRead:::.plotReadOccurrences(df_dist[df_dist$type=="read",], cex=.5)  # look for over represented sequences
ShortRead:::.plotCycleBaseCall(perCycle$baseCall)   # per cycle base calls of nucleotides
ShortRead:::.plotCycleQuality(perCycle$quality)    # plot quality across cycles

dev.off() 








