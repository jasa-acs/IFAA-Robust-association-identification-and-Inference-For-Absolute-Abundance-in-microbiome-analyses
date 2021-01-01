## 
## 
## R functions to run the second simulation
##
##

# setwd("C:/Dropbox (UFL)/papers/IFAA/Rpack/R")

library(IFAA)

#-------------------------------------------------------------------
## define function to read data
#-------------------------------------------------------------------

getData=function(dataFile="IFAA_Rev1.csv"){
   microdata=data.matrix(read.csv(file=dataFile,header=T,na.strings=c(NA,"")))
   return(microdata)
  }

#-------------------------------------------------------------------
## define function to get the sample variable name
#-------------------------------------------------------------------
getRepName=function(){return("rep")}

#-------------------------------------------------------------------
## define function to extract sample
#-------------------------------------------------------------------
 getSample=function(data=getData(),repName=getRepName(),sample){
   micro=data[data[,repName]==sample,]
   if(sum(is.na(micro))>0)stop("There are missing data.")
   rm(data)
   return(micro)
  }

#-------------------------------------------------------------------
## define function to extract covariates variable
#-------------------------------------------------------------------
 getCovsPrefix=function(){
   return("x")
  }

#-------------------------------------------------------------------
## define function to get the prefix of microbiome variable
#-------------------------------------------------------------------
 getMprefix=function(){return("rawCount")}

#-------------------------------------------------------------------
## define function to get the id variable name
#-------------------------------------------------------------------
 getIDname=function(){return("id")}

#-------------------------------------------------------------------
#
## run zero-inflated gaussian analysis
#
#-------------------------------------------------------------------

runIFAA_simu=function(
  sample,
  refTaxa=NULL,
  x1permut,
  sequentialRun=F,
  nRefMaxForEsti=10,
  Mprefix=getMprefix(),
  covsPrefix=getCovsPrefix(),
  fwer=0.2,
  nonZeroThresh=0,
  varThresh=0,
  seed=3
){
  results=list()
  data=getSample(sample=sample)

  # get taxa variable names
  
  MVarNamLength=nchar(Mprefix)
  micros = sapply(substr(colnames(data),1,MVarNamLength), function(x) {grep(Mprefix, x)})
  microPositions=which(micros == 1)
  rm(micros)
  
  taxaNames=colnames(data)[microPositions]
  rm(microPositions)

  w=data[,c("id",taxaNames),drop=F]

  # get predictor data
  xVarNamLength=nchar(covsPrefix)
  predics = sapply(substr(colnames(data),1,xVarNamLength), function(x) {grep(covsPrefix, x)})
  predPositions=which(predics == 1)
  rm(predics)
  predNames=colnames(data)[predPositions]
  rm(predPositions,xVarNamLength)
  
  xData=data[,c("id",predNames),drop=F]
  dataSpars=mean(data[,"dataSpars"])
  
  rm(data)
  
  cat("dataSpars:",dataSpars,"\n")
 
  IFAAresul=IFAA(sample=sample,MicrobData=w,
         CovData=xData,linkIDname="id",
         testCov=c("x1"),ctrlCov=c("x2","x3"),
         nRef=40,
         refTaxa=refTaxa,
         nPermu=40,x1permut=x1permut,
         reguMethod="mcp",paraJobs=8,
         nRefMaxForEsti=nRefMaxForEsti,
         standardize=T,sequentialRun=sequentialRun,
         refReadsThresh=nonZeroThresh,
         SDThresh=varThresh,
         SDquantilThresh=0,
         balanceCut=0,
         fwerRate=fwer,bootB=500,seed=seed)

  fileName=paste0("IFAA_simu_sample_",sample,".RData")
  
  IFAAresul$dataSpars=dataSpars

  save(IFAAresul,file=fileName)
  return(IFAAresul)
}

#
## run simulation on 50 data sets
#
for(i in 1:50){
 results=runIFAA_simu(sample=i,x1permut=T,fwer=0.2)
 }
