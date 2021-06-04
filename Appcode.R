library(limma)
library(ggplot2)
library( taRifx )
library(ComplexHeatmap)
library(DT)
library(M3C)
library(qqman)
Sys.setenv(RGL_USE_NULL = TRUE)
library(plot3Drgl)
rgl.open()
library(circlize)
library(tibble)
library(fgsea)
library(dplyr)
library(topGO)
library(Rgraphviz)

#user input files
expression.file <- readline(prompt="Enter expression values file path: ")
phenodata.file <- readline(prompt="Enter pheno data file path: ")
BP_annotated.file <- readline(prompt="Enter expression values annotated with BP file path: ")

# Reading expression data and pheno data and filtering irrelevant rows or columns. 
exprs_data <- unique(read.csv(expression.file,sep="\t",header=TRUE,dec=","))  ## Reading expression data ##
pdata = read.csv(phenodata.file,sep="\t",header=TRUE) ## Reading phenoData ##
pdata = subset(pdata,Mutation != '2x' & Mutation != 'Unkown')
expr_annotated <- read.csv(BP_annotated.file,sep=",",header=TRUE,dec=",")  ## Reading annotated expression data for Manhattan plot##
exprs_data <- cbind(exprs_data,expr_annotated$BP)
colnames(exprs_data)[colnames(exprs_data) == 'expr_annotated$BP'] <- 'BP'
exprs_data <- exprs_data[,c(1:11,ncol(exprs_data),(12:ncol(exprs_data)-1))]


# A function to create expression values dataframe of the given samples
createValuesDF <- function(sample1,sample2,chromosome = 'all'){# chromosome -> to choose chrom. X,autosomal or all(default)
  sampleList <<- c(sample1,sample2)
  samNames <<- paste(sample1, sample2,sep = ' with ')
  
  #Create design matrix from pdata
  filter_samples = subset(pdata,Sex %in% sampleList) ## Filtering only for given samples
  groups <<- filter_samples$Sex
  sex <<- factor(groups, levels = sampleList) # Creating factor of sample types for plots' legends
  d.matrix <- model.matrix( ~ 0 + sex) # Creating design matrix #
  colnames(d.matrix) <- sampleList
  design <<- d.matrix
  
  # Create the dataframe of expression values of the desired samples to be fitted
  myData <<- subset(exprs_data, !duplicated(exprs_data$Probe_Id))
  if(chromosome == 'autosomal'){
    myData <<- subset(myData,!Chromosome %in% 'X' & !Chromosome %in% 'Y')
  }else if(chromosome == 'X'){
    myData <<- subset(myData,Chromosome %in% 'X')
  }
  
  chr <<- chromosome
  df = myData[,12:ncol(exprs_data)] # Filteration of expression values only for desired comparisons #
  rownames(df) = myData$Probe_Id # Specify matrix rownames as Probe ID #
  df = df[,colnames(df) %in% filter_samples$SampleIDs] # Specify matrix colnames as sample IDs #
}
                      ####################################

## A function to combine exp. values together with pheno data information ##
showTable <- function(dataframe){
  # extract sex and age data
  c <- subset(pdata,SampleIDs %in% colnames(dataframe))
  ts <- t(data.frame(c$Sex))
  ta <- t(data.frame(c$Age))
  # join both in a dataframe
  v <- rbind(ts,ta)
  colnames(v) <- colnames(dataframe)
  v = data.frame(v)
  # join sex-age dataframe(v) to values dataframe (dataframe)
  v2 <- plyr::rbind.fill(v,dataframe)
  rownames(v2) <- c('Sex','Age',rownames(dataframe))
  rows <<- myData[myData$Probe_Id %in% rownames(dataframe),]
  v2$Chromosome <- c('','',as.character(rows$Chromosome))
  v2 <- v2[,c(ncol(v2),1:ncol(v2)-1)]
}  
                        ##################################

## A function fit expression data of given samples into linear and create a statistical dataframe ##
fitData <- function(dataframe){ 
  # Fitting data to linear model and performing empirical Bayes statistics
  fit <- lmFit(dataframe, design) # Initial Fit model without contrast #
  contrast <- gsub(' with ','-',samNames)
  cont.matrix <- makeContrasts(contrasts = contrast, levels=design) # Creating contrast to identify DEGs # 
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2, trend=FALSE)
  fit.df <- topTable(fit2,number = nrow(fit2)) # Creating a data frame of Bayes statistics' values ##
  fit.df <- signif(fit.df[,c('logFC','P.Value','adj.P.Val')],4)
}

                              ######################
## A function to plot a volcano plot of samples log fold change against p value to distinguish ##
## differently expressed genes ##
volcanoPlotData <- function(dataframe,p.value,logFC){
  #Define colours
  colour <- ifelse(dataframe[,"P.Value"] < p.value & abs(dataframe[,"logFC"]) > logFC,"red",
                   (ifelse(dataframe[,"P.Value"] < p.value ,"green", (ifelse(abs(dataframe[,"logFC"]) > logFC,'blue',"grey")))))

  #Plot
  par(mar=c(4,4,4,8),xpd = T)
  with(dataframe, plot(dataframe[,"logFC"], -log10(dataframe[,"P.Value"]), pch=20, xlim = c(-max(abs(dataframe[,"logFC"]))-0.1, max(abs(dataframe[,"logFC"]))+0.1),
                       main= paste("Diff.expr. genes of ",samNames,'in ', chr, " Chromosome(s)"),
                       xlab = "log FC",ylab = "-log10 p value",col =colour))
  abline(h= -log10(p.value), v=c(-logFC,logFC),col= 'blue',lty=3) # Specify margins
  legend('topright',inset = c(-0.22,0),legend = c('differently expressed','above p value threshold',
                               'above logFC threshold','normally expressed'),pch = 16,
                              col = c('red','green','blue','grey'),cex = 0.7)
}

                            ######################
## A function to convert the data to pca ##
createPCA <- function(dataframe){
  # Change values dataframe to numeric
  dataframe <- japply(dataframe, which(sapply(dataframe, class)=="factor"), as.character )
  dataframe <- japply(dataframe, which(sapply(dataframe, class)=="character"), as.numeric )
  
  # Create PCA data
  fdata <- filter.fitDF(dataframe,pvalue.value = 0.05) #filter DEGs with p value <= 0.05 only
  dataframe <- subset(dataframe,rownames(dataframe) %in% rownames(fdata))
  
  pca <- prcomp(t(dataframe),center = T,scale = T)
  pca.var <- pca$sdev^2  # variance##
  pca.var.per <<- round(pca.var/sum(pca.var)*100,1) # variance percentage##
  pca.data <<- data.frame(X = pca$x[,1],Y = pca$x[,2],Z = pca$x[,3])
}
                            ######################
## A function to 2D plot PCA ##
pcaPlotData2 <- function(pcaData){
  myColors <- c('Female' = 'red', 'Klinefelter' = '#36c5e4', 'Male' = 'blue', 'Turner' = '#fa05c0')
  #myShapes <<- c('Female' = 21, 'Klinefelter' = 7, 'Male' = 17, 'Turner' = 18)
  ggplot(data = pcaData,aes(x=X,y=Y,colour = sex))+ geom_point(alpha = 0.8,size = 3,pch = 20)+
    xlab(paste("PC",1, " => ",pca.var.per[1],"%",sep = ""))+
    ylab(paste("PC",2, " => ",pca.var.per[2],"%",sep = ""))+
    scale_color_manual(values = myColors)+
    ggtitle(paste("PCA of Expression data of ", chr, " chromosome(s) genes"))
}
                            ######################
## A function to 3D plot PCA ##
pcaPlotData3 <- function(pcaData){
 myColors <- c('Female' = 'red', 'Klinefelter' = '#36c5e4', 'Male' = 'blue', 'Turner' = '#fa05c0')
 rgl.close()
 pch3d(x=pcaData$X,y=pcaData$Y,z=pcaData$Z, color = subset(myColors,names(myColors) %in% sex),
    xlab= paste("PC",1, " => ",pca.var.per[1],"%",sep = ""),
    ylab= paste("PC",2, " => ",pca.var.per[2],"%",sep = ""),
    zlab= paste("PC",3, " => ",pca.var.per[3],"%",sep = ""),pch=20,colkey = F,cex = 0.6)
 legend3d("right", legend = sampleList, pch = 20,col = c(myColors[sampleList[1]],myColors[sampleList[2]]), cex=0.55)
 
}
                            #######################
## A function to generate a tSNE plot ##
tsnePlot <- function(dataframe){
  fdata <- filter.fitDF(dataframe,pvalue.value= 0.05) #filter DEGs with p value <= 0.05 only
  dataframe <- subset(dataframe,rownames(dataframe) %in% rownames(fdata))
  tsne(dataframe,perplex = 1,labels = sex, seed = 1,legendtextsize = 10)
}
                            #######################
## A function to generate a heatmap values of the given sample types ##
heatmapData <- function(dataframe){
  # scale data values
  dataframe <- scale(dataframe)
  dataframe <- t(scale(t(dataframe)))
  defitDF <- filter.fitDF(dataframe,pvalue.value = 0.05) #filter DEGs with p value <= 0.05 only
  dataframe <- subset(dataframe,rownames(dataframe) %in% rownames(defitDF)) 
  
  #convert dataframe to matrix
  dataframe <- as.matrix(dataframe)
  
  #colour-code row names
  myColors <- c('Female' = 'red', 'Klinefelter' = '#36c5e4', 'Male' = 'blue', 'Turner' = '#fa05c0')
  colour <- character()
  for(sam in groups){colour <- append(colour,myColors[[sam]])}
  
  Heatmap(t(dataframe),# transpose matrix to make samples more readeable at rows
          column_title = paste(chr,' chromosome(s) diff.expressed Genes'), row_title = 'Samples',
          row_names_side = 'left',row_names_gp = gpar(col = colour, cex = 0.5),
          row_dend_width = unit(3,'cm'),column_names_gp = gpar(cex = 0.4),
          show_column_names = FALSE,row_labels =groups, column_labels = as.numeric(row.names(dataframe)),
          show_column_dend = F, column_dend_reorder = F, row_dend_reorder = F,
          heatmap_legend_param = list(
                                title = "Transcription \n",
                                legend_height = unit(3, "cm")))

 }

                    ###################################
## A function to generate a Manhattan plot for the given samples (containing location information)
manPlotData <- function(dataframe){
  manTable <- data.frame(rownames(dataframe))
  colnames(manTable) <- 'SNP'     
  rows <- subset(rows,Probe_Id %in% manTable$SNP)  # defining probe-ID as SNP column
  manTable <- cbind(manTable,rows[,c('Chromosome','BP')])  # adding chromosome and gene location columns
  
  #chromosome data as numerical and proper naming
  manTable$Chromosome <- as.character(manTable$Chromosome)
  manTable$Chromosome <- plyr::revalue(manTable$Chromosome,c('X'= 23,'Y'= 24))
  manTable$Chromosome[!manTable$Chromosome %in% c(1:25)] <- 0
  manTable$Chromosome <- as.numeric(manTable$Chromosome)

  manTable$P.Value = fitData(dataframe)$P.Value  #adding p values column
  manTable <- manTable[!is.na(manTable$BP),]
  manTable <<- manTable[manTable$Chromosome != 0,] #filtering false chromosomes and gene locations
  manhattan(manTable, chr="Chromosome", bp="BP", snp="SNP", p="P.Value",chrlabs =  c(1:22, "X", "Y", "MT"))
}

                    ###################################
## A function to filter tables' samples and volcano plot according to logFC & p value given threshold##
filter.fitDF <- function(dataframe,logFC.value = NULL,pvalue.value = NULL){
  f <- fitData(dataframe)
  if(! is.null(logFC.value)){
    logFC.value <- abs(logFC.value)
    f <- subset(f, logFC <= - logFC.value | logFC >= logFC.value)
  }
  if(! is.null(pvalue.value)){
    f <- subset(f, P.Value <= pvalue.value)
  }
  f
}

                    ####################################
## A function for generating gene ontology by TopGo method#
goData <- function(dataframe){
  fdat <- filter.fitDF(dataframe,pvalue.value =0.05) #filter DEGs with p value <= 0.05 only
  geneList = c(fdat$P.Value)
  names(geneList)= rownames(fdat)
  data(geneList)
  
  sampleGOdata <- new("topGOdata",
                      description = "Simple session", ontology = "BP",
                      allGenes = geneList, geneSel = topDiffGenes,
                      annot = annFUN.db, affyLib = "illuminaHumanv4.db")
  
  resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
  resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
  resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")
  
  allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                     classicKS = resultKS, elimKS = resultKS.elim,
                     orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)
  pValue.classic <<- score(resultKS)
  pValue.elim <<- score(resultKS.elim)[names(pValue.classic)]
  gstat <<- termStat(sampleGOdata, names(pValue.classic))
  gSize <<- gstat$Annotated / max(gstat$Annotated) * 4
  allRes
}

                    ###################################
# valuesDF <- createValuesDF('Klinefelter','Turner')
# s <- showTable(valuesDF)
# fitDF <- fitData(valuesDF)
# volcanoPlotData(fitDF,0.05,0.5)
# pcaPlotData2(createPCA(valuesDF))
# pcaPlotData3(createPCA(valuesDF))
# tsnePlot(valuesDF)
# heatmapData(valuesDF)
# manPlotData(valuesDF)
# f <- filter.fitDF(valuesDF,pvalue.value =0.05)
# g = goData(valuesDF)
# goPlotData(g)
