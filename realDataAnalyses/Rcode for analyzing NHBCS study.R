
#
## R code for analyzing the NHBCS study data
#
realDataSets=get(load("realData.RData"))
data_NHBCS=realDataSets$data_NHBCS

#-------------------------------------------------------------------
## IFAA method to analyze the NHBCS data
#-------------------------------------------------------------------

library(IFAA)

runIFAA=function(
  data,
  Mprefix,
  covsPrefix,
  fwer
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
  
  IFAAresul=IFAA(MicrobData=w,CovData=xData,linkIDname="id",
                 testCov=c("x1","x2","x3"),
                 nRef=40,nPermu=40,reguMethod="mcp",
                 refReadsThresh=0.2,SDThresh=0,
                 SDquantilThresh=0,balanceCut=0,
                 fwerRate=fwer,bootB=500,seed=3)
  
  results$analysisResults=IFAAresul$analysisResults$estByCovList$x1
  rm(IFAAresul)
  
  # return 
  return(results)
}

NHBCSresu=runIFAA(data=data_NHBCS,Mprefix="rawCount",
                  covsPrefix="x",fwer=0.25)

NHBCSresu$analysisResults

#-------------------------------------------------------------------
## ANCOM method to analyze the NHBCS data.
## Instructions for installing the ANCOM R package can be
## found in the supplimental materials.
#-------------------------------------------------------------------

library("ancom.R")

runAncom=function(
  data,
  Mprefix,
  groupVar,
  imputeValue=0.001, 
  sig=0.25, 
  multcorr=1, 
  tau=0.02, 
  theta=0.1
){
  
  results=list()
  
  MVarNamLength=nchar(Mprefix)
  
  # get taxa variable names
  micros = sapply(substr(colnames(data),1,MVarNamLength), function(x) {grep(Mprefix, x)})
  microPositions=which(micros == 1)
  rm(micros)
  
  taxaNames=colnames(data)[microPositions]
  rm(microPositions)
  
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
  
  # impute zero-valued data points
  ancomData[,subsetTaxaNames][ancomData[,subsetTaxaNames]==0]=imputeValue
  
  ancomResu<-ANCOM(OTUdat=ancomData,sig=sig,multcorr=multcorr,tau=tau,theta=theta)
  
  rm(ancomData)

  results$TaxaNamByAncom=ancomResu$detected
  return(results)
}

ancomResul=runAncom(data=data_NHBCS,Mprefix="rawCount",groupVar="x1")
ancomResul$TaxaNamByAncom

#-------------------------------------------------------------------
## Spearman test method to analyze the NHBCS data.
#-------------------------------------------------------------------

SpearmanTest=function(
  data,
  AAprefix,
  xVar,
  adjAlpha
){
  
  results=list()
  
  # get taxa data
  namLength=nchar(AAprefix)
  micros = sapply(substr(colnames(data),1,namLength), function(x) {grep(AAprefix, x)})
  nTaxa = length(which(micros == 1))
  
  taxaName=paste0(AAprefix,1:nTaxa)
  
  taxaData=data[,taxaName]
  
  rm(micros,taxaName)

  # get predictor data
  
  covariates=as.matrix(data[,xVar,drop=F])
  microbData=as.matrix(taxaData)
  
  nCovariates=ncol(covariates)
  nTaxa=ncol(microbData)

  microbName=colnames(microbData)
  
  ## test correlations with AA
  pMatrix=matrix(NA,nrow=nTaxa,ncol=nCovariates)

  for(i in 1:nTaxa){
    for(j in 1:nCovariates){
      corTest.ij=cor.test(x=microbData[,i], y=covariates[,j],
                          alternative = c("two.sided"),
                          method = c("spearman"),
                          exact=F)
      pMatrix[i,j]=corTest.ij$p.value
      rm(corTest.ij)
    }
  }

  raw.pVec=c(t(pMatrix))
  rm(pMatrix)

  fdrP=p.adjust(p=raw.pVec,method="fdr")
  results$sigTaxa.AA=microbName[fdrP<=adjAlpha]

  
  ## test correlations with RA
  pMatrix=matrix(NA,nrow=nTaxa,ncol=nCovariates)
  
  microbData=microbData/rowSums(microbData)

  for(i in 1:nTaxa){
    for(j in 1:nCovariates){
      corTest.ij=cor.test(x=microbData[,i], y=covariates[,j],
                          alternative = c("two.sided"),
                          method = c("spearman"),
                          exact=F)
      pMatrix[i,j]=corTest.ij$p.value
      rm(corTest.ij)
    }
  }

  raw.pVec=c(t(pMatrix))
  rm(pMatrix)
  
  fdrP=p.adjust(p=raw.pVec,method="fdr")

  results$sigTaxa.RA=microbName[fdrP<=adjAlpha]
  
  return(results)
}

SpearmanTest(data=data_NHBCS,AAprefix="rawCount",xVar="x1",adjAlpha=0.25)

  