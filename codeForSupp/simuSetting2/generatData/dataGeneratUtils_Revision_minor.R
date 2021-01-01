# setwd("C:/Dropbox (UFL)/papers/IFAA/Rpack/dataGeneratR")

#
## data generation utility functions
#

library(rockchalk) ## for function mvrnorm() 

source("TruParam_Revision_minor.R")

#
## generate covariates data including the main predictor 
## of interest
#

getPredics=function(
 nPredics=getNPredics(),
 totRepN,
 rhoPredics=getTruPredicsRho(),
 SDpredics=getTruPredicsSigma(),
 outlierBound=c(5,6)
 ){
  PredicsM=matrix(NA,nrow=totRepN,ncol=nPredics)

  meanPredics=rep(0,nPredics)

  rhoPredicsMatrix=matrix(NA,nrow=nPredics,ncol=nPredics)

  for (i in 1:nPredics){
    for (j in 1:nPredics){
     rhoPredicsMatrix[i,j]=rhoPredics^(abs(i-j))
   }
  }
  SDpredicsVec=rep(SDpredics,nPredics)
  varPredics=diag(SDpredicsVec)%*%rhoPredicsMatrix%*%diag(SDpredicsVec)

  PredicsM=mvrnorm(n=totRepN,mu=meanPredics,Sigma=varPredics)
  
  PredicsM[,1]=abs(PredicsM[,1])
  PredicsM[,2]=(PredicsM[,2]>0)+0
  
  #mixtureProb=rbinom(totRepN,1,0.97) # with outlier
  mixtureProb=rep(1,totRepN) # no outlier
  PredicsM[,3]=PredicsM[,3]*mixtureProb+
    (1-mixtureProb)*runif(totRepN,outlierBound[1],outlierBound[2])
  
  # adding intercept to the covariates
  PredicsWithInter=cbind(rep(1,totRepN),PredicsM)

  return(PredicsWithInter)
 }


#
## generate  multivariate Bernoulli distribution 
## for zero-inflation data structure
#

conBern=function(
 nTaxa=getNTaxa(),
 totRepN,
 startFromScratch
 ){
    zeroInflatM=matrix(NA,nrow=totRepN,ncol=nTaxa)

    mvnMatForCut=mvrnorm(n=totRepN,mu=rep(0,nTaxa),Sigma=diag(nTaxa))
    if(!startFromScratch){
     #cutVec=rep(qnorm(0.84),nTaxa)
     #cutVec=rep(qnorm(0.74),nTaxa)
     #cutVec=rep(qnorm(0.64),nTaxa)
     cutVec=rep(qnorm(0.44),nTaxa)
     }
    if(startFromScratch){
      cutVec=rep(-10^6,nTaxa)
      }
    
    zeroInflatM=t(t(mvnMatForCut)>=cutVec)+0
    return(zeroInflatM)
 }

#
## generate Ci to determine sequencing depth
#
getSeqDep=function(
  totRepN,
  predicsWithInter=getPredics(totRepN=totRepN),
  continSeqDepDeno=getContinuSeqDepDeno()
){

  lowBound=continSeqDepDeno[1]
  hiBound=continSeqDepDeno[2]
  
  confoundCoef=10
  mean=min(confoundCoef*predicsWithInter[,2])
  perTurb=mvrnorm(n=1,mu=rep(mean,totRepN),Sigma=diag(rep(0.1,totRepN)))
  
  Ci=exp(-confoundCoef*(predicsWithInter[,2])+perTurb)

  Ci1inver=pmin((1/Ci),hiBound)
  Ci1inver=pmax(Ci1inver,lowBound)
  Ci=1/Ci1inver
  return(Ci)
}

#
## generate multivariate zero-inflated taxa data
#

Microb=function(
 nTaxa=getNTaxa(),
 totRepN,
 startFromScratch,
 predicsWithInter=getPredics(totRepN=totRepN),
 beta=getTruBeta(),
 zeroInflatM=conBern(totRepN=totRepN,startFromScratch=startFromScratch),
 realDataParam=getRealDataParam(),
 Ci=getSeqDep(totRepN=totRepN),
 varianceMax=4,
 varianceRareTaxa=2
 ){
  
 if(sum(beta[-1,]!=0)==0) stop("error: all regression coefficients are 0")

 truAAMat=matrix(NA,nrow=totRepN,ncol=nTaxa)
 epsMat=matrix(NA,nrow=totRepN,ncol=nTaxa)
 
 bi=rnorm(n=totRepN)*sqrt(realDataParam$biVar)

 muNormM=predicsWithInter%*%beta
 
 varMuNorm=diag(var(muNormM))

 # generate multivariate normal residual error
 varianceVec=pmin(realDataParam$sigma2k,varianceMax)
 varianceVec[(realDataParam$rareTaxa)]=varianceRareTaxa
 varianceVec=pmax((varianceVec-varMuNorm),0.001)

 epsMat=mvrnorm(n=totRepN,mu=rep(0,nTaxa),Sigma=diag(varianceVec))

 for (i in 1:totRepN){
    truAA=exp(muNormM[i,]+bi[i]+epsMat[i,])
    truAAMat[i,]=truAA*(zeroInflatM[i,])
    }

 obsAAMat=floor(truAAMat*Ci)

 results=list()
 results$obsAAMat=obsAAMat
 
 return(results)
 }

# micro=Microb()

#
## extract existing multivariate zero-inflated taxa data 
#

useExtMicrob=function(
  extData="IFAA_Rev_origin1.csv",
  nTaxa=getNTaxa(),
  zeroInflatM=conBern(totRepN=totRepN)
){
  
  allData=read.csv(extData,header=T)
  rawCountsName=paste0("rawCount",1:nTaxa)
  extM=allData[,rawCountsName]
  obsAAMat=extM*zeroInflatM
  
  allData[,rawCountsName]=obsAAMat
  allData[,"dataSpars"]=mean(apply(obsAAMat,c(1,2),function(x)x==0))
  return(allData)
}

# micro=useExtMicrob()


#
### function for data generation in batches
#
generateData_Rev=function(
  batch,
  startFromScratch=T
 )
 {
  #
  ## calculate the total number of subjects and the number of taxa
  #
  nDataset=getNdataPerBatch()# number of data sets in each batch
  nSub=getNSub() # number of subjects in each data set
  nTaxa=getNTaxa() # number of taxa
  nPredics=getNPredics() # number of all predictors

  totRepN=nDataset*nSub

  # load the regression coefficients alpha
  beta=getTruBeta()

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
  ## generate zero-inflation data structure
  #
  zeroStruc=conBern(totRepN=totRepN,startFromScratch=startFromScratch)
  Ci=getSeqDep(totRepN=totRepN,predicsWithInter=PredicsWithInter)
  
  #
  ## generate microbiome data
  #
  if(startFromScratch){
   microbData=Microb(
    totRepN=totRepN,
    startFromScratch=startFromScratch,
    predicsWithInter=PredicsWithInter,
    beta=beta,
    zeroInflatM=zeroStruc,
    Ci=Ci
    )
   seqDepLow=getContinuSeqDepDeno()[1]
   seqDepHi=getContinuSeqDepDeno()[2]
   }else{
     extDataSetName=paste0("IFAA_Rev_origin",(((batch-1)%%5)+1),".csv",sep="")
     microbData=useExtMicrob(
       extData=extDataSetName,
       nTaxa=getNTaxa(),
       zeroInflatM=conBern(totRepN=totRepN,startFromScratch=startFromScratch)
       )
     microbData[,"rep"]=microbData[,"rep"]+floor((batch-1)/5)*40
     dataSetName=paste0("IFAA_Rev",batch,".csv",sep="")
     write.csv(microbData,file=dataSetName, row.names=F)
     return("Done")
     }
   cat("microbiome generation done.","\n")
   microbiome=microbData$obsAAMat
   rm(microbData)
   
  #
  ## calculate average microbiome data sparsity
  #
  
  dataSpars=mean(apply(microbiome,c(1,2),function(x){x==0}))

  #
  ## combine data components
  #

  dataMatrix=cbind(repid,PredicsWithInter[,-1],microbiome)
  predicsName=paste0("x",1:nPredics)
  rawCountsName=paste0("rawCount",1:nTaxa)

  cat("data merge done.","\n")
  
  varNames <-c(c("rep","id"),predicsName,rawCountsName)

  colnames(dataMatrix)<-varNames

  #
  ## save the data in the data folder
  #
  dataSetName=paste0("IFAA_Rev_origin",batch,".csv",sep="")

  write.csv(cbind(dataSpars,dataMatrix),file=dataSetName, row.names=F)
  }

# generateData_Rev(batch=1)

