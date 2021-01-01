## 
## 
## R functions to obtain MLEs for the linear model 
## of the microbiome data on Arsenic based on zero-inflated 
## normal distribution
##
##

# setwd("C:/Dropbox (UFL)/papers/taxaSelect/R")

source("./generatData/TruParam_minorR.R")

#
## function to load packages
#

loadPacks=function(){
library(IFAA) ## for function IFAA() 
library("ancom.R") ## for ancom analysis
library(metagenomeSeq) ## for ZIG analysis
library(biomformat)    ## for data processing for ZIG analysis
library(DESeq2) ## for Deseq analysis
library(edgeR) ## for edgeR analysis
}

loadPacks()

#-------------------------------------------------------------------
## define function to read data
#-------------------------------------------------------------------
 getData=function(dataFile="../dataGeneratR/IFAAancom1.csv"){
     microdata=data.matrix(read.csv(file=dataFile,header=T))
   return(microdata)
  }

# class(getData())

#-------------------------------------------------------------------
## define function to get the rep variable name
#-------------------------------------------------------------------
getRepName=function(){return("rep")}

#getRepName()

#-------------------------------------------------------------------
## define function to extract sample
#-------------------------------------------------------------------
 getSample=function(data=getData(),repName=getRepName(),sample=sample){
   micro=data[data[,repName]==sample,]
   if(sum(is.na(micro))>0)stop("There are missing data.")
   rm(data)
   return(micro)
  }

# class(getSample(sample=1))

#-------------------------------------------------------------------
## define function to extract covariates variable
#-------------------------------------------------------------------
 getCovsPrefix=function(){
   return("x")
  }

# getCovsPrefix()

#-------------------------------------------------------------------
## define function to get the prefix of microbiome variable
#-------------------------------------------------------------------
 getMprefix=function(){return("rawCount")}

# getMprefix()

#-------------------------------------------------------------------
## define function to get the id variable name
#-------------------------------------------------------------------

 getIDname=function(){return("id")}

# getIDname()


#-------------------------------------------------------------------
#
## data sparsity check
#
#-------------------------------------------------------------------

dataSparsCheck=function(
data,
Mprefix=getMprefix()
){
  results=list()

  # get the original sample size
  nSub=nrow(data)
  MVarNamLength=nchar(Mprefix)

  # get taxa variable names
  micros = sapply(substr(colnames(data),1,MVarNamLength), function(x) {grep(Mprefix, x)})
  microPositions=which(micros == 1)
  rm(micros)

  taxaNames=colnames(data)[microPositions]
  rm(microPositions)

  w=data[,taxaNames]
  rm(data,taxaNames)

  # check zero taxa and subjects with zero taxa reads
  numTaxaNoReads=length(which(Matrix::colSums(w)==0))
  print("which(Matrix::colSums(w)==0):")
  print(which(Matrix::colSums(w)==0))
    if(numTaxaNoReads>0){
    warning(paste("There are",numTaxaNoReads,"taxa without any sequencing reads"))
   }
  rm(numTaxaNoReads)

  numSubNoReads=length(which(Matrix::rowSums(w)==0))
  if(numSubNoReads>0){
    warning(paste("There are",numSubNoReads,"subjects without any sequencing reads"))
   }
 rm(numSubNoReads,w)
}

# dataSparsCheck(data=getSample(sample=1))


#-------------------------------------------------------------------
#
## take basic data information
#
#-------------------------------------------------------------------

dataInfo=function(
data,
Mprefix=getMprefix(),
covsPrefix=getCovsPrefix(),
binPredInd=getBinPredInd(),
refReadsThresh=0.1,
SDThresh=0.001,
SDquantilThresh=0,
balanceCut=0.1
){
  results=list()

  # get the original sample size
  nSub=nrow(data)

  dataSparsLvl=data[1,(colnames(data)=="dataSpars")]

  MVarNamLength=nchar(Mprefix)
  # get taxa variable names
  micros = sapply(substr(colnames(data),1,MVarNamLength), function(x) {grep(Mprefix, x)})
  microPositions=which(micros == 1)
  rm(micros)

  nTaxa = length(microPositions)
  taxaNames=colnames(data)[microPositions]
  rm(microPositions)

  # get the biggest index number of the taxa name
  taxaNameNum=rep(NA,nTaxa)

  for (i in 1:nTaxa){
   taxa.i.NamLength=nchar(taxaNames[i])
   taxaNameNum[i]=as.numeric(substr(taxaNames[i],(MVarNamLength+1),taxa.i.NamLength)) 
   }
  maxTaxaNameNum=max(taxaNameNum)
  rm(taxaNameNum)

  w=data[,taxaNames]

  # check zero taxa and subjects with zero taxa reads
  # print("percentage of taxa present for at least 25% subjects:")
  # print(sum((Matrix::colSums(w>0)/nSub)>=0.25))

  results$non0perct=Matrix::colSums(w>0)/nSub

  taxaOverThresh=taxaNames[(Matrix::colSums(w>0)>=nSub*refReadsThresh)]
  if(length(taxaOverThresh)==0){
    stop(paste("There are no taxa with presence over the threshold:",refReadsThresh,
      ". Try lower the reference taxa reads threshold"))
    }

  # check the sd threshold
  sdTaxaOverThresh=rep(0,length(taxaOverThresh))
  for (i in 1:length(taxaOverThresh)){
    taxa.i=w[,taxaOverThresh[i]]
    if(sum(taxa.i>0)>1){
     sdTaxaOverThresh[i]=sd(taxa.i[(taxa.i>0)])
     }
    }

  results$sdTaxa=sdTaxaOverThresh
 
  TaxaOverSdThresh=taxaOverThresh[(sdTaxaOverThresh>=SDThresh)]
  if(length(TaxaOverSdThresh)==0){
    stop(paste("There are no taxa with SD over the SD threshold:",SDThresh,
      ". Try lower the SD threshold"))
    }
  rm(taxa.i,taxaOverThresh)

  # check the sd quantile threshold
  sdAllTaxa=rep(0,nTaxa)
  for (i in 1:nTaxa){
    taxaAll.i=w[,taxaNames[i]]
    posTaxaAll.i=taxaAll.i[(taxaAll.i>0)]
    if(length(posTaxaAll.i)>1){sdAllTaxa[i]=sd(posTaxaAll.i)}
    }
  goodRefTaxaCandi=TaxaOverSdThresh[(TaxaOverSdThresh>=quantile(sdAllTaxa,probs=SDquantilThresh))]
  rm(sdAllTaxa,posTaxaAll.i,TaxaOverSdThresh)

  if(length(goodRefTaxaCandi)==0){
    stop(paste("There are no taxa with SD over the SD quantile threshold:",SDquantilThresh,
      ". Try lower the SD quantile threshold"))
    }

  nSubNoReads=length(which(Matrix::rowSums(w>0)==0))
  rm(w)

  # get predictor data
  xVarNamLength=nchar(covsPrefix)

  predics = sapply(substr(colnames(data),1,xVarNamLength), function(x) {grep(covsPrefix, x)})
  predPositions=which(predics == 1)
  predNames=colnames(data)[predPositions]
  EName=predNames[1]
  nPredics=length(predNames)
  rm(predics,predPositions)

 # find the pairs of binary preds and taxa for which the assocaiton is not identifiable
 if(length(binPredInd)>0){
  firstBinPredNam=paste0(covsPrefix,binPredInd)
  binPredStart=which(predNames%in%firstBinPredNam)
  binPredStop=length(predNames)
  allBinPred=predNames[binPredStart:binPredStop]
  nBinPred=length(allBinPred)
  rm(predNames)

  taxaAndBinIndexNoInt=vector()
  taxaNoBin=c()
  taxaBalanceBin=c()

  for(i in 1:nTaxa){
    for(j in 1:nBinPred){
      twoColumns.ij=data[,c(taxaNames[i],allBinPred[j])]
      nNonZero=length(which(twoColumns.ij[,1]>0))
      sumOfBin=sum(twoColumns.ij[(twoColumns.ij[,1]>0),2])
      rm(twoColumns.ij)
      if(sumOfBin%in%c(0,nNonZero)){
        taxaNoBin=c(taxaNoBin,taxaNames[i])
        index.ij=(i-1)*nPredics+binPredInd+j-1
        taxaAndBinIndexNoInt=c(taxaAndBinIndexNoInt,index.ij)
        }
      if(min(sumOfBin,(nNonZero-sumOfBin))>=balanceCut*nNonZero){
        taxaBalanceBin=c(taxaBalanceBin,taxaNames[i])
        }
      }
    }
   results$taxaAndBinIndexNoInt=taxaAndBinIndexNoInt
   rm(taxaAndBinIndexNoInt,allBinPred,data)
   

   # remove unbalanced taxa across binary variables
   goodRefTaxaCandi=goodRefTaxaCandi[!(goodRefTaxaCandi%in%taxaNoBin)]

   # keep balanced taxa
   goodRefTaxaCandi=goodRefTaxaCandi[(goodRefTaxaCandi%in%taxaBalanceBin)]

   rm(taxaNoBin)
  }
  # return 
  results$taxaNames=taxaNames
  rm(taxaNames)
  results$goodRefTaxaCandi=goodRefTaxaCandi
  rm(goodRefTaxaCandi)
  results$nTaxa=nTaxa
  results$nSub=nSub
  results$nSubNoReads=nSubNoReads
  results$nPredics=nPredics
  results$maxTaxaNameNum=maxTaxaNameNum
  results$EName=EName
  results$dataSparsLvl=dataSparsLvl
  return(results)
}

# data.info=dataInfo(data=getSample(sample=1))


#-------------------------------------------------------------------
## function for cross validation using glmnet package
#-------------------------------------------------------------------

runAncom=function(
  data,
  taxaNames,
  groupVar,
  imputeValue=0.001, #ancom paper used 0.001
  sig=0.2, 
  multcorr=1, 
  tau=0.02, 
  theta=0.1
){
  
  results=list()
  
  origin.taxaNames=taxaNames
  rm(taxaNames)
  # check zero taxa and subjects with zero taxa reads
  TaxaNoReads=which(Matrix::colSums(data[,origin.taxaNames])==0)
  subsetTaxaNames=origin.taxaNames
  if(length(TaxaNoReads>0)){
    subsetTaxaNames=origin.taxaNames[-TaxaNoReads]
  }
  rm(TaxaNoReads)
  
  # remove rows with no sequencing reads
  subNoReads=which(rowSums(data[,subsetTaxaNames])==0)
 
   if(length(subNoReads)>0){
    data=data[-subNoReads,]
  }
  rm(subNoReads)
  
  ancomData=as.data.frame(data[,c(subsetTaxaNames,groupVar)])
  rm(data)
  
  ancomData[,subsetTaxaNames][ancomData[,subsetTaxaNames]==0]=imputeValue

  ancomResu<-ANCOM(OTUdat=ancomData,sig=sig,multcorr=multcorr,tau=tau,theta=theta)
  
  rm(ancomData)
  print("W:")
  print(ancomResu$W)
  
  TaxaNamByAncom=ancomResu$detected
  
  TaxaByAncom=(origin.taxaNames%in%TaxaNamByAncom)+0
  rm(ancomResu)
  
  results$TaxaNamByAncom=TaxaNamByAncom
  rm(TaxaNamByAncom)
  
  results$TaxaByAncom=TaxaByAncom
  rm(TaxaByAncom)
  
  return(results)
}

# runAncom(data=getSample(sample=1),sig=0.3,taxaNames=data.info$taxaNames,groupVar=data.info$EName)


#-------------------------------------------------------------------
#
## run zero-inflated gaussian analysis
#
#-------------------------------------------------------------------

runIFAA=function(
  data,
  Mprefix=getMprefix(),
  covsPrefix=getCovsPrefix(),
  fwer=0.2,
  seed=1
){
  results=list()
  
  # get the original sample size
  MVarNamLength=nchar(Mprefix)
  
  # get taxa variable names
  micros = sapply(substr(colnames(data),1,MVarNamLength), function(x) {grep(Mprefix, x)})
  microPositions=which(micros == 1)
  rm(micros)
  
  taxaNames=colnames(data)[microPositions]
  rm(microPositions)

  # make data in the format for ZIG
  w=data[,c("id",taxaNames),drop=F]

  # get predictor data
  xVarNamLength=nchar(covsPrefix)
  predics = sapply(substr(colnames(data),1,xVarNamLength), function(x) {grep(covsPrefix, x)})
  predPositions=which(predics == 1)
  rm(predics)
  predNames=colnames(data)[predPositions]
  rm(predPositions,xVarNamLength)
  
  # annotate the predictor data
  xData=data[,c("id",predNames),drop=F]

    IFAAresul=IFAA(MicrobData=as.data.frame(w),CovData=as.data.frame(xData),
               linkIDname="id",
               testCov=c("x1"),ctrlCov=c(),
               nRef=40,
               nPermu=40,reguMethod="mcp",paraJobs=8,
               refReadsThresh=0,
               SDThresh=0,
               SDquantilThresh=0,
               balanceCut=0,
               fwerRate=fwer,bootB=500,seed=seed)

  results$IFAASelect=IFAAresul$analysisResults$selecTaxaFWER
  results$allIFAAresults=IFAAresul
  rm(IFAAresul)
  
  # return 
  return(results)
}

# IFAAresu=runIFAA(data=getSample(sample=1))

#-------------------------------------------------------------------
#
## run zero-inflated gaussian analysis
#
#-------------------------------------------------------------------

runZig=function(
data,
Mprefix=getMprefix(),
covsPrefix=getCovsPrefix(),
adjAlphaLevel=0.2
){
  results=list()

  # get the original sample size
  MVarNamLength=nchar(Mprefix)

  # get taxa variable names
  micros = sapply(substr(colnames(data),1,MVarNamLength), function(x) {grep(Mprefix, x)})
  microPositions=which(micros == 1)
  rm(micros)

  origin.taxaNames=colnames(data)[microPositions]
  rm(microPositions)
  taxaNames=origin.taxaNames

  # check zero taxa and subjects with zero taxa reads
  TaxaNoReads=which(Matrix::colSums(data[,origin.taxaNames])==0)
  if(length(TaxaNoReads>0)){
   taxaNames=origin.taxaNames[-TaxaNoReads]
   }
  rm(TaxaNoReads)
  nTaxa = length(taxaNames)

  # get the index numbers of taxa names
  taxaNameNum=rep(NA,nTaxa)
  for (i in 1:nTaxa){
   taxa.i.NamLength=nchar(taxaNames[i])
   taxaNameNum[i]=as.numeric(substr(taxaNames[i],(MVarNamLength+1),taxa.i.NamLength)) 
   }
  rm(MVarNamLength)

  # remove rows with no sequencing reads
  subNoReads=which(rowSums(data[,taxaNames])==0)
  if(length(subNoReads)>0){
    data=data[-subNoReads,]
    }
  rm(subNoReads)
  nSub=nrow(data)

  # make data in the format for ZIG
  w=data[,taxaNames]
  countsData=as.data.frame(t(w))
  rm(w)
  colnames(countsData)=paste0("sub",seq(1:nSub))
  rm(nSub)

  # annotate taxa names
  taxaNum=as.data.frame(taxaNameNum)
  rm(taxaNameNum)
  rownames(taxaNum)=taxaNames
  rm(taxaNames)
  taxaAnnot=AnnotatedDataFrame(taxaNum)
  rm(taxaNum)

  # get predictor data
  xVarNamLength=nchar(covsPrefix)
  predics = sapply(substr(colnames(data),1,xVarNamLength), function(x) {grep(covsPrefix, x)})
  predPositions=which(predics == 1)
  rm(predics)
  predNames=colnames(data)[predPositions]
  rm(predPositions,xVarNamLength)

  # annotate the predictor data
  xData=as.data.frame(data[,predNames])
  rm(data)
  rownames(xData)=colnames(countsData)
  colnames(xData)=predNames
  xDataAnnot=AnnotatedDataFrame(xData)

  print("dim(countsData):")
  print(dim(countsData))
  
  # create the final data for ZIG
  newData=newMRexperiment(counts=countsData,phenoData=xDataAnnot,featureData=taxaAnnot)
  rm(countsData,xDataAnnot,taxaAnnot)

  # normalization
  dataP=cumNormStat(newData)
  newData=cumNorm(newData, p=dataP)
  rm(dataP) 
  normFactor=normFactors(newData)
  normFactor=log2(normFactor/(median(normFactor))+1)

  # build the model
  mod=model.matrix(~data.matrix(xData))
  colnames(mod)=c("(intercept)",predNames)
  rm(predNames)
  settings=zigControl(maxit = 10, verbose = TRUE)

  newData=cumNorm(newData)
  fit = fitZig(newData, mod)

  ZIGresul = MRcoefs(fit, number = nTaxa)
  rm(fit,nTaxa)
  selectTaxa=rownames(ZIGresul)[(ZIGresul$adjPvalues<=adjAlphaLevel)]
  rm(ZIGresul)
  zigSelect=(origin.taxaNames%in%selectTaxa)+0
  rm(origin.taxaNames,selectTaxa)

  # return 
  results$zigSelect=zigSelect
  rm(zigSelect)
  return(results)
 }
 
 # zigResu=runZig(data=getSample(sample=1))


 # zigResu$zigSelect

 # length(zigResu$zigSelect)

#-------------------------------------------------------------------
#
## function to run edgeR analysis
#
#-------------------------------------------------------------------
runEdgeR=function(
data,
Mprefix=getMprefix(),
covsPrefix=getCovsPrefix(),
adjAlphaLevel=0.2
){
  results=list()

  # get the original sample size
  MVarNamLength=nchar(Mprefix)

  # get taxa variable names
  micros = sapply(substr(colnames(data),1,MVarNamLength), function(x) {grep(Mprefix, x)})
  microPositions=which(micros == 1)
  rm(micros)

  origin.taxaNames=colnames(data)[microPositions]
  rm(microPositions)
  taxaNames=origin.taxaNames

  # check zero taxa and subjects with zero taxa reads
  TaxaNoReads=which(Matrix::colSums(data[,origin.taxaNames])==0)
  if(length(TaxaNoReads>0)){
   taxaNames=origin.taxaNames[-TaxaNoReads]
   }
  rm(TaxaNoReads)

  # remove rows with no sequencing reads
  subNoReads=which(rowSums(data[,taxaNames])==0)
  if(length(subNoReads)>0){
    data=data[-subNoReads,]
    }
  rm(subNoReads)

  # make data in the format for DESeq2
  countsData=t(data[,taxaNames])
  rownames(countsData)=taxaNames
  rm(taxaNames)

  # get predictor data
  xVarNamLength=nchar(covsPrefix)
  predics = sapply(substr(colnames(data),1,xVarNamLength), function(x) {grep(covsPrefix, x)})
  predPositions=which(predics == 1)
  rm(predics)
  predNames=colnames(data)[predPositions]
  rm(predPositions)
  EName=predNames[1]
  rm(xVarNamLength)

  # annotate the predictor data
  xData=data[,predNames]
  rm(data)
  
print("dim(countsData):")
print(dim(countsData))

  # build the model
  model<-DGEList(counts=countsData,group=xData)  
  rm(countsData)
  
  z = edgeR:::calcNormFactors(model, method = "RLE")
  
  results$edgeRfail=0
  if (!all(is.finite(z$samples$norm.factors))) {
      results$edgeRfail=1
      results$edgerSelect=rep(0,length(origin.taxaNames))
      return(results)  
     }
  # Estimate dispersions
  z1 = estimateCommonDisp(z)
  z2 = estimateTagwiseDisp(z1)
  rm(z,z1)
  EdgerTestt = exactTest(z2)
  rm(z2)
  print("results$edgeRfail:")
  print(results$edgeRfail)

  # extract the analysis results
  fdrAdj=decideTestsDGE(EdgerTestt,adjust.method="fdr",p.value=adjAlphaLevel)
  rm(EdgerTestt)
  selectTaxa=rownames(fdrAdj)[(fdrAdj[,1]!=0)]
  rm(fdrAdj)
  edgerSelect=(origin.taxaNames%in%selectTaxa)+0
  rm(origin.taxaNames,selectTaxa)

  # return 
  results$edgerSelect=edgerSelect
  rm(edgerSelect)
  return(results)
 }
 
 # edgeRResu=runEdgeR(data=getSample(sample=1))

 # edgeRResu$edgerSelect

 # length(edgeRResu$edgerSelect)

#-------------------------------------------------------------------
#
## function to run DESeq2 analysis
#
#-------------------------------------------------------------------
runDeseq2=function(
data,
Mprefix=getMprefix(),
covsPrefix=getCovsPrefix(),
adjAlphaLevel=0.2
){
  results=list()

  # get the original sample size
  MVarNamLength=nchar(Mprefix)

  # get taxa variable names
  micros = sapply(substr(colnames(data),1,MVarNamLength), function(x) {grep(Mprefix, x)})
  microPositions=which(micros == 1)
  rm(micros)

  origin.taxaNames=colnames(data)[microPositions]
  rm(microPositions)
  taxaNames=origin.taxaNames

  # check zero taxa and subjects with zero taxa reads
  TaxaNoReads=which(Matrix::colSums(data[,origin.taxaNames])==0)
  if(length(TaxaNoReads>0)){
   taxaNames=origin.taxaNames[-TaxaNoReads]
   }
  rm(TaxaNoReads)

  # remove rows with no sequencing reads
  subNoReads=which(rowSums(data[,taxaNames])==0)
  if(length(subNoReads)>0){
    data=data[-subNoReads,]
    }
  rm(subNoReads)

  # make data in the format for DESeq2
  countsData=t(data[,taxaNames])+1
  rownames(countsData)=taxaNames
  rm(taxaNames)

  # get predictor data
  xVarNamLength=nchar(covsPrefix)
  predics = sapply(substr(colnames(data),1,xVarNamLength), function(x) {grep(covsPrefix, x)})
  predPositions=which(predics == 1)
  rm(predics)
  predNames=colnames(data)[predPositions]
  rm(predPositions)
  EName=predNames[1]
  rm(xVarNamLength)

  # annotate the predictor data
  xData=as.data.frame(data[,predNames])
  colnames(xData)=predNames
  rm(data)

  # build the model
  mod=model.matrix(~data.matrix(xData))
  colnames(mod)=c("(intercept)",predNames)
  rm(predNames)

  xData[,"x1"]=as.factor(xData[,"x1"])
  
  dds <- DESeqDataSetFromMatrix(countData=countsData,colData=xData,design=~x1)
  
  rm(countsData,xData,mod)

  # run analysis
  results$DESeqFail=0
  nbinom=0
  suppressWarnings(runDeseq<- try(DESeq(dds, quiet = TRUE), silent = TRUE))
  if (inherits(runDeseq, "try-error")) {
    # If the parametric fit failed, try the local.
    suppressWarnings(runDeseq<- try(DESeq(dds, fitType = "local", quiet = TRUE), 
                                silent = TRUE))
    if (inherits(runDeseq, "try-error")) {
      # If local fails, try the mean
      suppressWarnings(runDeseq<- try(DESeq(dds, fitType = "mean", quiet = TRUE), 
                                  silent = TRUE))
    }
    if (inherits(runDeseq, "try-error")) {
      dds <- estimateSizeFactors(dds)
      dds <- estimateDispersionsGeneEst(dds)
      dispersions(dds) <- mcols(dds)$dispGeneEst
      runDeseq <- suppressWarnings(try(nbinomWaldTest(dds,quiet = TRUE),
                                       silent = TRUE))
      nbinom=1
      }
    if (inherits(runDeseq, "try-error")) {
      # If still bad, quit with error.
      results$DESeqFail=1
      results$DESeqSelect=rep(0,length(origin.taxaNames))
      return(results)
    }
  }
  rm(dds)

  # extract the regression results
  if(!nbinom){
    res<-results(runDeseq,name=EName)
  }
  if(nbinom){
    res<-results(runDeseq)
  }
  rm(runDeseq,EName)
  selectTaxa=rownames(res)[(res$padj<=adjAlphaLevel)]
  rm(res)
  DESeqSelect=(origin.taxaNames%in%selectTaxa)+0
  rm(origin.taxaNames,selectTaxa)

  # return 
  results$DESeqSelect=DESeqSelect
  rm(DESeqSelect)
  return(results)
 }
 
 # DESeqResu=runDeseq2(data=getSample(sample=2))

 # DESeqResu$DESeqSelect

 # length(DESeqResu$DESeqSelect)



#-------------------------------------------------------------------
## Function to calculate recall, precision, F1 and type 1
#-------------------------------------------------------------------

PerformInd=function(coefEstVec,truVec){
  results=list()

  if(length(coefEstVec)!=length(truVec)){
    stop("Performance indices cannot be calculated due to different vector dimensions")
   }
  TruPos=sum(which(coefEstVec!=0)%in%which(truVec!=0))
  FalsePos=sum(which(coefEstVec!=0)%in%which(truVec==0))
  TruNeg=sum(which(coefEstVec==0)%in%which(truVec==0))
  FalseNeg=sum(which(coefEstVec==0)%in%which(truVec!=0))
  
  recall=TruPos/(TruPos+FalseNeg)
  precision=TruPos/(TruPos+FalsePos)
  F1=2*(recall*precision)/(recall+precision)
  typeI=FalsePos/(FalsePos+TruNeg)
  modelSize=length(which(coefEstVec!=0))
     
  if(modelSize==0){
    recall=0
    precision=1
    F1=0
    typeI=0
    }

  # return results
  results$recall=recall
  results$precision=precision
  results$F1=F1
  results$typeI=typeI
  results$modelSize=modelSize

  return(results)
 }

# PerformInd(c(1,2,0,0,3),c(3,4,5,3,5))



#-------------------------------------------------------------------
## Run regularization estimation
#-------------------------------------------------------------------

Regulariz=function(
  sample,
  compareMethods
 ){
  results=list()

  # load data
  data=getSample(sample=sample)
   
  # load data info
  dataSparsCheck(data=data)
  dataForEst=dataInfo(data=data,refReadsThresh=refReadsThresh)
  nSub=dataForEst$nSub
  results$taxaNames=dataForEst$taxaNames
  EName=dataForEst$EName
  nPredics=dataForEst$nPredics
  nTaxa=dataForEst$nTaxa
  rm(dataForEst)

  nPredicsWithInt=nPredics+1

  # load true regression coeficients
  
  truAlpha=getTruBetaEcho()
  truAlphaVec=c(t(truAlpha))
  rm(truAlpha)
  
  nAlphaWithInt=length(truAlphaVec)
  truAlphaVecNoInt=truAlphaVec[-seq(1,nAlphaWithInt,by=nPredicsWithInt)]
  rm(truAlphaVec)
  nAssociations=length(truAlphaVecNoInt)
  truSelection=truAlphaVecNoInt
  results$truModelSize=sum(truAlphaVecNoInt!=0)
  rm(truAlphaVecNoInt)
  
  results$truSelection=truSelection

   # ancom analysis
  if("ancom" %in% compareMethods){
   ancomResul=runAncom(data=data,taxaNames=results$taxaNames,groupVar=EName)
   rm(EName)
   results$ancomResu=ancomResul$TaxaByAncom
   rm(ancomResul)
   indices=PerformInd(results$ancomResu,truSelection)

   results$recallAncom=indices$recall
   results$precisionAncom=indices$precision
   results$F1Ancom=indices$F1
   results$typeIancom=indices$typeI
   results$modelSizeAncom=indices$modelSize
   }

  # zig analysis
  if("zig" %in% compareMethods){
    zigResu=runZig(data=data)
   
    indices=PerformInd(zigResu$zigSelect,truSelection)
    rm(zigResu)

    results$recallZig=indices$recall
    results$precisionZig=indices$precision
    results$F1Zig=indices$F1
    results$typeIZig=indices$typeI
    results$modelSizeZig=indices$modelSize
    }
   
  # IFAA analysis
   if("IFAA" %in% compareMethods){
     IFAAResu=runIFAA(data=data,seed=sample)
     
     results$allIFAAresults=IFAAResu$allIFAAresults
     indices=PerformInd(IFAAResu$IFAASelect,truSelection)
       
     rm(IFAAResu)
     
     results$recallIFAA=indices$recall
     results$precisionIFAA=indices$precision
     results$F1IFAA=indices$F1
     results$typeIIFAA=indices$typeI
     results$modelSizeIFAA=indices$modelSize
     }
   
  # edgeR analysis
  if("edgeR" %in% compareMethods){
    edgeRResu=runEdgeR(data=data)
    results$edgeRfail=edgeRResu$edgeRfail
    indices=PerformInd(edgeRResu$edgerSelect,truSelection)
    rm(edgeRResu)

    results$recallEdgeR=indices$recall
    results$precisionEdgeR=indices$precision
    results$F1EdgeR=indices$F1
    results$typeIEdgeR=indices$typeI
    results$modelSizeEdgeR=indices$modelSize
    }

  # DESeq2 analysis
  if("DESeq2" %in% compareMethods){
    DESeqResu=runDeseq2(data=data)
    results$DESeqFail=DESeqResu$DESeqFail
    indices=PerformInd(DESeqResu$DESeqSelect,truSelection)
    rm(DESeqResu)

    results$recallDESeq2=indices$recall
    results$precisionDESeq2=indices$precision
    results$F1DESeq2=indices$F1
    results$typeIDESeq2=indices$typeI
    results$modelSizeDESeq2=indices$modelSize
    }
  rm(data,truSelection)

  # return results
  
  results$nSub=nSub
  results$nTaxaTru=nTaxaTru
  results$nTaxa=nTaxa
  results$nPredics=nPredics
  return(results)
 }

# Regularization=Regulariz(sample=1,nRef=3,realData=F,screenMethod=c("lasso"),finalPenalMethod=c("MCP"))

#
## run regularization and save the data
#

IFAAsimu=function(
sample,
compareMethods=c("IFAA","ancom","zig","edgeR","DESeq2")
){
  # set.seed(1)
   start.time = proc.time()[3] 
   loadPacks(ancomRun=ancomRun)

   #
   ## get regularized estimate
   #
   Regularization=Regulariz(sample=sample,
                            compareMethods=compareMethods
                            )

   execution.time = (proc.time()[3] - start.time)/60

   #
   ## save the output
   #

   output.data = paste0("simuResu_",sample,".RData")  
   save(Regularization, sample, execution.time, file=output.data)
  }
 

#
## run the simulation on 100 data sets
#

for (i in 1:100){IFAAsimu(sample=i)}

