
library(tidyverse)
library(readxl) # part of the tidyverse but needs to be explictly loaded
library(DESeq2) # For differential expression - from bioconductor
library(AnnotationDbi)  # for gene annotations - from bioconductor
library(org.Mm.eg.db) # for mouse annotation specifically - from bioconductor
library(Cairo) # to save the plots

#####################################
#### Read in the data for Deseq #####
#####################################

# read in your count table
count.table = read_csv("DESeq_Data/count_table.csv")
row.names(count.table) = count.table$ENSEMBL
count.table = dplyr::select(count.table, -ENSEMBL)

# read in the experimental design
experimental.design <- read_excel("DESeq_Data/experimental_design.xlsx")
row.names(experimental.design) = experimental.design$sample
experimental.design = dplyr::select(experimental.design, -sample)

# it is critical that the names are identical
# they should be in the same order  before setting names
rownames(experimental.design) = colnames(count.table) 
experimental.design$treatment = as.factor(experimental.design$treatment)
experimental.design$time = as.factor(experimental.design$time)


##########################################
#### Analysis for Multi-Factor Design ####
##########################################


#to make all the comparisons across groups
experimental.design$group = factor(paste0(experimental.design$treatment, experimental.design$time))

# Generate the DESeq Data Set 
dds = DESeqDataSetFromMatrix(countData = count.table,
                             colData = experimental.design,
                             design = ~ group)

# Run the differential expression
dds = DESeq(dds)
#resultsNames(dds)
res = results(dds, alpha=0.05)
#res = results(dds, alpha=0.05, lfcThreshold = 2)
summary(res)


###########################
#### Get the Contrasts ####
###########################

# extract contrasts from the groups
# format: c("factor", "experimental", "control")

cmatrix = tidyr::expand(experimental.design, nesting(treatment, time)) # get all possible combinations
cmatrix = dplyr::filter(cmatrix, treatment != "BMDM") # BMDM is our control - so we remove it from the table
cmatrix$group = paste0(cmatrix$treatment, cmatrix$time)

# we will split the time points
h24 = dplyr::filter(cmatrix, time == "h24")
h6 = dplyr::filter(cmatrix, time == "h6")

# Create a Matrix of Contrasts
l1 = c("group")
l2 = c("BMDMh24")
l3 = h24$group
m1 = cbind(l1, l2, l3)

# Matrix for 6hrs
l4 = c("group")
l5 = c("BMDMh6")
l6 = h6$group
m2 = cbind(l4, l5, l6)

# Put them both together
m = rbind(m2,m1)

# Write a little function to get our contrasts
contrast_matrix<-function(x)
{
  results(dds, contrast = x) 
}

# Get a list of DESeq Results from Our Contrasts of Interest
res_list= apply(m, 1, FUN = contrast_matrix)
names(res_list) <- paste("Control_V_", m[,3], sep = "")


# Write a little function to get Result Summaries
apply(m, 1, function(x) summary(results(dds, contrast = x, alpha = 0.05)))

# Save the Summarized Results with capture.output
cat("DESeq Summarized Results ", file = "summarized_results.txt")
cat("\n\nComparisons Made in Order\n\n ", file = "summarized_results.txt", append = TRUE)
capture.output(names(res_list), file = "summarized_results.txt", append = TRUE)
# export anova test output
capture.output(apply(m, 1, function(x) summary(results(dds, contrast = x, alpha = 0.05))), file = "summarized_results.txt", append = TRUE)


#############################################
### Work with the Contrasts a little bit ####
#############################################

# Make a Sub-directory to save the outputs
main_dir <- getwd()
contrast_dir <- "DESeq_Contrasts_Results"
output_dir <- file.path(main_dir, contrast_dir)

# Check if the sub-directory exists 
if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Directory already exists!")
}

# Save all the files out
sapply(names(res_list), 
       function (x) write.csv(res_list[[x]], file=paste0(output_dir, "/", x, ".csv") )   )

setwd(output_dir)
list.filenames = list.files(path = output_dir, pattern = ".csv")
list.data = list()  #make an empty list to receive data
#create a loop to read in the data
for (i in 1:length(list.filenames))
{
  list.data[[i]] = read_csv(list.filenames[i] )
}
names(list.data)<-list.filenames  # add the names of your data to the list


# Do some filtering
sig = lapply(list.data, function(x) dplyr::filter(x, padj < 0.05))
up_sig = lapply(sig, function(x) dplyr::filter(x, log2FoldChange > 0))
down_sig = lapply(sig, function(x) dplyr::filter(x, log2FoldChange < 0))

# Save all the files out
sapply(names(sig), 
       function (x) write.csv(sig[[x]], file=paste0(output_dir, "/", "significant_", x) ) )
sapply(names(up_sig), 
       function (x) write.csv(up_sig[[x]], file=paste0(output_dir, "/", "significant_UP_", x) ) )
sapply(names(down_sig), 
       function (x) write.csv(down_sig[[x]], file=paste0(output_dir, "/", "significant_DOWN_", x) ) )



##################################################
#### Prep the Data for Export with Annotation ####
##################################################

# Get all the data as a DF
full_contrast_DE_results= as.data.frame(res_list)

# Get the full list of all possible arguments
#columns(org.Mm.eg.db)

full_contrast_DE_results$symbol <- mapIds(org.Mm.eg.db,
                                          keys = row.names(full_contrast_DE_results),
                                          column = "SYMBOL",
                                          keytype = "ENSEMBL",
                                          multiVals = "first")
full_contrast_DE_results$entrez <- mapIds(org.Mm.eg.db, 
                                          keys = row.names(full_contrast_DE_results),
                                          column = "ENTREZID",
                                          keytype = "ENSEMBL",
                                          multiVals = "first")
full_contrast_DE_results$name <- mapIds(org.Mm.eg.db, 
                                          keys = row.names(full_contrast_DE_results),
                                          column = "GENENAME",
                                          keytype = "ENSEMBL",
                                          multiVals = "first")


#change the first column to the gene ID
names(full_contrast_DE_results)[2] = "TotalBaseMean"
#remove all the columns with basemean - these are redundant
full_contrast_DE_results = full_contrast_DE_results[, -grep("baseMean_*", colnames(full_contrast_DE_results))]
# write this out for later
write.csv(full_contrast_DE_results, file= "Annotated_Full_Contrast_DE_Results.csv")


###########################
#### Descriptive Plots ####
###########################

plots_dir <- "DESeq_Descriptive_Plots"
plots_output_dir <- file.path(main_dir, plots_dir)

if (!dir.exists(plots_output_dir)){
  dir.create(plots_output_dir)
} else {
  print("Directory already exists!")
}


# Plot Counts of Specific Genes
# Most significant gene for Example
plotCounts(dds, gene=which.min(res$padj), intgroup="treatment")
d = plotCounts(dds, gene=which.min(res$padj), intgroup="treatment",
               returnData = TRUE)


#PCA Plots
rld = rlog(dds) # vst transformation might also be approriate 
pcaData= plotPCA(rld, intgroup= c("treatment", "time"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
colnames(pcaData) = c("PC1", "PC2", "Group", "Treatment", "Time",  "Name") # Make them look nicer


#save the individual plots as high resolution images
setwd(plots_output_dir)
CairoPNG(filename = "DESeq_Descriptive_Plots%03d.png", 
         units = "in",
         height = 8,
         width = 8, 
         pointsize = 12, 
         dpi = 300,
         bg = "white")
# Plot Counts of Specific Genes as ggplot object
ggplot(d, aes(x=treatment, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))
# PCA Plot as ggplot object
ggplot(pcaData, aes(x = PC1, y = PC2, color = Treatment, shape = Time)) +
  geom_point(size = 8) +
  scale_fill_brewer(palette="Spectral") +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_minimal() +
  theme(title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(face = "plain", size = 12),
        legend.title = element_text(face = "bold", hjust = 1),
        legend.text = element_text(face = "plain", size = 12),
        axis.line.x = element_line(size = .8, lineend = "square"),
        axis.line.y = element_line(size = .8, lineend = "butt"),
        panel.border = element_rect(fill = NA))
# Dispersion Estimates Plot
plotDispEsts(dds)
# MA Plot
ylim <- c(-5,5)
drawLines <- function() abline(h=c(-2,2),col="dodgerblue",lwd=2)
plotMA(res, ylim=ylim); drawLines()
dev.off() 

