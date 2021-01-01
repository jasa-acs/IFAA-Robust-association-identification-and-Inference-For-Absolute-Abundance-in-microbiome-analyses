# setwd("C:/Dropbox (UFL)/Papers/IFAA/dataGeneratR")

#
## data generation utility functions
#

source("TruParam_minorR.R")

#
## generate covariates data 
#

getPredics=function(
  nPredics=getNPredics(),
  totRepN,
  quantileCut=0.5
){
  PredicsM=matrix(NA,nrow=totRepN,ncol=nPredics)
  
  x=rnorm(n=totRepN)
  
  # convert the first predictor into a binary covariate
  PredicsM[,1]=(x>quantile(x,probs=quantileCut))+0
  
  # adding intercept to the covariates
  PredicsWithInter=cbind(rep(1,totRepN),PredicsM)
  
  return(PredicsWithInter)
}

#
## generate sequencing depth
#
getSeqDep=function(
  totRepN,
  predicsWithInter=getPredics(totRepN),
  seqDepthDeno=getSeqDepthDeno()
  ){
   seqDepth=rep(NA,totRepN)
   treatGrp=which(predicsWithInter[,2]==1)
   controlGrp=which(predicsWithInter[,2]==0)
   for(i in 1:length(treatGrp)){
     seqDepth[treatGrp[i]]=30
     #seqDepth[treatGrp[i]]=90
     }
   for(i in 1:length(controlGrp)){
     seqDepth[controlGrp[i]]=30
     #seqDepth[controlGrp[i]]=18
     #seqDepth[controlGrp[i]]=9
     #seqDepth[controlGrp[i]]=6
     }
  return(1/seqDepth)
 }

# getSeqDep()

microAncom=function(
 nTaxa,
 totRepN,
 predicsWithInter=getPredics(totRepN),
 truBeta=getTruBetaEcho(),
 seqDep=getSeqDep(totRepN,predicsWithInter=getPredics(totRepN))
 ){

 if(nTaxa!=nrow(truBeta))
   stop("Number of taxa does not match the true number of parameters")
    
 if(sum(truBeta!=0)==0) stop("error: all ANCOM model Poisson means are 0")

 ## generate the null abundance matrix
 rawReads=matrix(NA,nrow=totRepN,ncol=nTaxa)

 # generate the mean 
 muNormM=predicsWithInter%*%t(truBeta)
 #print(predicsWithInter)

 # generate the observed counts
 echoMat=matrix(NA,nrow=totRepN,ncol=nTaxa)
 for (i in 1:totRepN){
    echoMat[i,]=rpois(n=nTaxa,lambda=muNormM[i,])
    rawReads[i,]=floor(seqDep[i]*echoMat[i,])
    }
 rm(echoMat)

 return(rawReads)
 }

# microAncom()

#
### function for data generation in batches with ANCOM paper setting
#
generateDataAncom=function(
  batch
 ){
  #
  ## calculate the total number of subjects and the number of taxa after log-ratio transformation
  #
  nDataset=getNdataPerBatch()# number of data sets in each batch
  nSub=getNSub() # number of subjects in each data set
  nTaxa=getNTaxa() # number of taxa
  nPredics=getNPredics() # number of all predictors
  totRepN=nDataset*nSub

  # load the Poisson means 
  truBeta=getTruBetaEcho()

  #
  ## generate replicate and sample identifiers
  #
  repid=matrix(NA,totRepN,2)
  repid[,1]=rep(seq((batch-1)*nDataset+1,batch*nDataset),each=nSub)
  repid[,2]=rep(seq(1,nSub),nDataset)

  #
  ## generate predictors data 
  #
  PredicsWithInter=getPredics(totRepN=totRepN)

  #
  ## generate sequecing depth
  #
  seqDep=getSeqDep(totRepN=totRepN,predicsWithInter=PredicsWithInter)

  #
  ## generate microbiome data
  #

  microAncomData=microAncom(
    nTaxa=nTaxa,
    totRepN=totRepN,
    predicsWithInter=PredicsWithInter,
    truBeta=truBeta,
    seqDep=seqDep
    )

  #
  ## combine data components
  #

  dataMatrix=cbind(repid,PredicsWithInter[,-1],microAncomData)
  predicsName=paste0("x",1:nPredics)
  rawCountsName=paste0("rawCount",1:nTaxa)

  varNames <-c(c("rep","id"),predicsName,rawCountsName)

  colnames(dataMatrix)<-varNames

  #
  ## save the data in the data folder
  #
 
  dataSetName=paste0("IFAAancom",batch,".csv",sep="")
  write.csv(cbind(dataMatrix),file=dataSetName, row.names=F)
 }



