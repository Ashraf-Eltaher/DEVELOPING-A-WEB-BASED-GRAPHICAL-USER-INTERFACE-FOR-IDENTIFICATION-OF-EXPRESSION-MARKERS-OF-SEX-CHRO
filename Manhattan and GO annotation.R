library("biomaRt")
# import ensemble annotation data set
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

annot.df = getBM(attributes = c('start_position','entrezgene_id'),filters=c('entrezgene_id'),values=exprs_data[,'EntrezID'],mart=ensembl,uniqueRows = T)
  
#write.table(bp,"D:/MXP project/Datasets MXP3/annotation_columns.txt")

#merge the two data sets' columns and order genes according to Entrez Id
annotated <- unique(merge(exprs_data,annot.df,by.x = 'EntrezID', by.y = 'entrezgene_id' ,all.x = T))
annotated <- annotated[,c(2:9,1,10,11,98,12:ncol(annotated)-1)]

# update column names to use in Manhattan plot
colnames(annotated)[colnames(annotated) == 'Probe_Id'] <- 'SNP'
colnames(annotated)[colnames(annotated) == 'start_position'] <- 'BP'

# save data set in a table
write.table(annotated,"D:/MXP project/Datasets MXP3/exp. data annotated.txt",sep = '\t', row.names = F)

# extract extra genes in a new table
mart.export <- readline(prompt="Enter biomart exported file path: ")
rest_genes <- unique(read.csv(mart.export,header=TRUE))
rest_genes <- rest_genes[unique(rest_genes$SNP),]
write.table(rest_genes,"mart_export.txt", row.names = F)

# merge two tables by BP (gene start position)
qw <- unique(merge(expr_annotated,rest_genes,by = 'BP' ,all.x = T))
write.table(qw,"qw.txt",sep = '\t', row.names = F) # save #

# extract our original data set genes only to make the gene ontology gene list file
g <- rest_genes[expr_annotated$SNP %in% rest_genes$SNP,]


gene_list <- exprs_data[,1:11]
write.csv(gene_list,"gene list.gmt") #save gmt format#