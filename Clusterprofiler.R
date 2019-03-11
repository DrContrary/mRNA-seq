
library(tidyverse)
library(readxl) # part of the tidyverse but needs to be explictly loaded
library(DESeq2) # For differential expression - from bioconductor
library(AnnotationDbi)  # for gene annotations - from bioconductor
library(org.Mm.eg.db) # for mouse annotation specifically - from bioconductor
library(Cairo) # to save the plots
library(clusterProfiler)


###########################
#### Read In the Files ####
###########################

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

# Move into the Directory and Start Loading
setwd(output_dir)
list.filenames = list.files(path = output_dir, pattern = '^Control.*\\.csv')
list.data_files = list()  #make an empty list to receive data
#create a loop to read in the data
for (i in 1:length(list.filenames))
{
  list.data_files[[i]] = read_csv(list.filenames[i], 
                                  col_names = c("ENSEMBL", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
}
names(list.data_files)<-list.filenames  # add the names of your data to the list
setwd(main_dir)

###################################################
#### Use AnnotationDbi to Get Gene Annotations ####
###################################################

## get the keys
keys <- keys(org.Mm.eg.db)
## lookup gene symbol and unigene ID for the 1st 6 keys
annotate = select(org.Mm.eg.db, keys=keys, columns = c("SYMBOL","GENENAME", "ENSEMBL"))
# DESeq Uses ENSEMBL IDs so any that are NA can be omitted
# Not neccesary but reduces extraneous data
annotate = annotate %>%  na.omit()  

# Join our Annotation df to each of the contrast result tables
dat = lapply(list.data_files, function(x) dplyr::left_join(annotate, x, by = "ENSEMBL"))
# Filter these so up and down regualted genes can be analyzed individually
dat_sig = lapply(dat, function(x) dplyr::filter(x, padj < 0.05))
dat_up = lapply(dat_sig, function(x) dplyr::filter(x, log2FoldChange > 1))
dat_down = lapply(dat_sig, function(x) dplyr::filter(x, log2FoldChange < 1))

########################################
#### Background List for Enrichment ####
########################################

#background list
setwd("./DESeq_Data")
count.table = read_csv("count_table.csv")
row.names(count.table) = count.table$ENSEMBL
background.list = dplyr::select(count.table, -ENSEMBL)

# Map the ENSEMBL IDs to ENTREZIDS
background.list$ENTREZID <- mapIds(org.Mm.eg.db,
                                 keys = row.names(background.list),
                                 column = "ENTREZID",
                                 keytype = "ENSEMBL",
                                 multiVals = "first")

# Get the Row Sums
background.list= dplyr::mutate(background.list, rowsum = rowSums(background.list[,1:36]))
# Remove rows with less than one count
background.list = dplyr::filter(background.list, rowsum > 0)

#######################################
#### ClusterProfiler Gene Ontology ####
#######################################

# Can also use BP and MF
# "all" can also be used but memory intensive for large data
ontology = "CC"

# Write a little function to perform the enrichment
go_enrich<-function(x)
{
  enrichGO(gene          = x$ENTREZID,
           universe      = background.list$ENTREZID,
           OrgDb         = org.Mm.eg.db,
           ont           = ontology,
           pAdjustMethod = "BH",
           pvalueCutoff  = 0.01,
           qvalueCutoff  = 0.05,
           readable      = TRUE) 
}

# Perform the Enrichment for GO terms
# For GO - UP and Down are Often Split
ego_down = lapply(dat_down, FUN = go_enrich)
ego_up = lapply(dat_up, FUN = go_enrich)


##############################
#### Save the CSV Results ####
##############################

# Make a Sub-directory to save the outputs
setwd(main_dir)

ClusterProfiler_dir <- "ClusterProfiler_Analysis"
cluster_output_dir <- file.path(main_dir, ClusterProfiler_dir)

# Check if the sub-directory exists 
if (!dir.exists(cluster_output_dir)){
  dir.create(cluster_output_dir)
} else {
  print("Directory already exists!")
}
setwd(cluster_output_dir)


# Write out all the Results of the S4 object to file
# Save the down Regulated Enrichments
sapply(names(ego_down), 
       function (x) write.csv(ego_down[[x]]@result, 
                              file=paste0(cluster_output_dir, "/", "significant_down_", ontology, "_", x) ) ) 
# Save the Up Regulated Enrichments
sapply(names(ego_up), 
       function (x) write.csv(ego_down[[x]]@result, 
                              file=paste0(cluster_output_dir, "/", "significant_up_", ontology, "_", x) ) ) 

############################
#### Generate the Plots ####
############################

# Get the individual Dotplots
plots_down = lapply(ego_down, function(x) dotplot(x))
plots_up = lapply(ego_up, function(x) dotplot(x))

# Other Plotting Options
#cnetplot(ego, categorySize="pvalue", foldChange=geneList)
## categorySize can be scaled by 'pvalue' or 'geneNum'
#emapplot(ego)


# Change the Plot Names to the List Names
# Makes it easier to read later
for (i in 1:length(plots_down))
{
  plots_down[[i]][["labels"]]$title  = names(plots_down)[i]
}
#Run for the Up Regulated Enrichments
for (i in 1:length(plots_up))
{
  plots_up[[i]][["labels"]]$title  = names(plots_up)[i]
}


# Save the plots as a multipage PDF  
p1 <- marrangeGrob(plots_down, nrow=1, ncol=1)
p2 <- marrangeGrob(plots_up, nrow=1, ncol=1)
ggsave(filename = paste0("significant_down_", ontology, "_", "ClusterProfiler_Plots", ".pdf"), p1)
dev.off()
ggsave(filename = paste0("significant_up_", ontology, "_", "ClusterProfiler_Plots", ".pdf"), p2)
dev.off()


#########################################
#### ClusterProfiler KEGG Enrichment ####
#########################################

# Write a little function to perform the enrichment
KEGG_enrich<-function(x)
{
  enrichKEGG(gene         = x$ENTREZID,
             universe     = background.list$ENTREZID,
             organism     = 'mmu',
             pvalueCutoff = 0.05) 
}

# Perform the Enrichment for KEGG Pathways
# KEGG up and down are performed together - only significant genes
kk = lapply(dat_sig, FUN = KEGG_enrich)

# Write out all the Results of the S4 object to file
# Save the down Regulated Enrichments
sapply(names(kk), 
       function (x) write.csv(kk[[x]]@result, 
                              file=paste0(cluster_output_dir, "/", "significant_", "KEGG", "_", x) ) ) 


plots_kk = lapply(kk, function(x) dotplot(x))

# Change the Plot Names to the List Names
# Makes it easier to read later
for (i in 1:length(plots_kk))
{
  plots_kk[[i]][["labels"]]$title  = names(plots_kk)[i]
}

p3 <- marrangeGrob(plots_kk, nrow=1, ncol=1)
ggsave(filename = paste0("significant_KEGG_", "ClusterProfiler_Plots", ".pdf"), p3)
dev.off()





