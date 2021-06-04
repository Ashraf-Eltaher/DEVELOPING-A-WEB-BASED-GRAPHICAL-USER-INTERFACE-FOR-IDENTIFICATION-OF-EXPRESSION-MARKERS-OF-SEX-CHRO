rm(list=ls())

library(lumi)
library(gdata)
library(lumiHumanIDMapping)
library(arrayQualityMetrics)

#user input files
DIR <- readline(prompt="Enter directory folder path: ")
probe <- readline(prompt="Enter raw probes file path: ")

#read expression file and annotation
dat = lumiR(file.path(DIR, probe), 
            annotationColumn=c('SEARCH_KEY', 'SYMBOL', 'CHROMOSOME', 'ILMN_GENE', 'DEFINITION', 'SYNONYMS', "ENTREZ_GENE_ID", "ONTOLOGY_COMPONENT", "ONTOLOGY_PROCESS", "ONTOLOGY_FUNCTION"), lib.mapping='lumiHumanIDMapping', dec=",", sep="\t") 
pdat = read.csv(file.path(DIR,"annotation.txt"),sep="\t")
pdat = subset(pdat,pdat$Ave.or.Ind == "Ind")
rownames(pdat) = pdat$SampleIDs


#Quality control
library(affy)
#plot logFC against log conc. 
MAplot(dat)
#Generate quality metrics
arrayQualityMetrics(dat, outdir=file.path(DIR, "QualityControlRaw/"), do.logtransform=TRUE) 


#lumiExpresso includes pre-processing packages from raw data till expression values 
#(BG correction,VST, Normalization and QC)
dat.norm = lumiExpresso(dat, varianceStabilize.param=list(method="vst"))
pData(dat.norm) = pdat
sampleNames(dat.norm) = colnames(dat.norm)
arrayQualityMetrics(dat.norm, outdir=file.path(DIR, "QualityControlNorm/"), force=TRUE)

dat.nAP = detectionCall(dat.norm, type="matrix")
take = c()
for(s in unique(pData(dat.norm)$Donor.Num.Tissue)){
  take = union(take, which(rowMeans(dat.nAP[, pData(dat.norm)$SampleIDs[which(pData(dat.norm)$Donor.Num.Tissue == s)], drop=F] == "P") == 1)) # take only probes, which are present in each patient for one group
}
dat.norm = dat.norm[take,] 

#Mapping IDs
mappingInfo = nuID2EntrezID(featureNames(dat.norm), lib.mapping='lumiHumanIDMapping', returnAllInfo=TRUE)
mappingInfo2 = nuID2IlluminaID(featureNames(dat.norm), species = "Human",lib.mapping='lumiHumanIDMapping')
mappingInfo3 = nuID2RefSeqID(featureNames(dat.norm), lib.mapping='lumiHumanIDMapping')
# only consider probes mapping to Entrez gene IDs
exprs_data = exprs(dat.norm)
exprs_data = merge(mappingInfo,exprs_data,by.x="row.names",by.y="row.names")
rownames(exprs_data) = exprs_data$Row.names
exprs_data = exprs_data[,-1]
exprs_data = merge(mappingInfo2,exprs_data,by.x="row.names",by.y="row.names")
rownames(exprs_data) = exprs_data$Row.names
exprs_data = exprs_data[,-1]
exprs_data = exprs_data[which(exprs_data$EntrezID!=""),]

# pre-processed expression values file
write.table(exprs_data,"Expression Corrected normalized vst all probes.txt",sep="\t",row.names=FALSE,quote=FALSE)