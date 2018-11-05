
library(readr) #to import the tab deliminated data from HTSeq

#########################
#   Read Tabular Files  #
#########################

list.filenames = list.files(pattern = ".tabular")
list.data = list()  #make an empty list to receive data
#create a loop to read in the data
for (i in 1:length(list.filenames))
{
  list.data[[i]] = read_delim(list.filenames[i], "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE )
}

names(list.data)<-list.filenames  # add the names of your data to the list
count_table = as.data.frame(list.data)


#########################
#     Format Table    #
#########################


#change the first column name to  gene_id
names(count_table)[1] = "gene_id"
count_table = count_table[, -grep(".X1", colnames(count_table))] #remove all the columns with .X2 - these are extra gene IDs


columns = colnames(count_table)
columns = gsub(pattern = ".*X6h", "", columns, perl = TRUE)
column = gsub(pattern = ".tabular.X2*", "", columns, perl = TRUE)
columns = paste0(column, '_6hr')

#now take them and add them back to the column names
names(count_table) = columns
names(count_table)[1] = "gene_id"
head(count_table)

write.csv(count_table, "count_table.csv")

