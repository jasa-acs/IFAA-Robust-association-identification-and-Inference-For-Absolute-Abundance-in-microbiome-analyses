#
### R code for summarizing simulation results
#

setwd("C:/Dropbox (UFL)/papers/IFAA/Rpack/R")

source("./generatData/TruParam_Revision_minor.R") # for getTruBeta(), etc
source("simu_Rev_minor.R") # for getMprefix()

analyzeResults=function(
    dataFile="IFAA_Rev1.csv",
    fwer,
    Mprefix=getMprefix()
){

 options(digits=3)
 options(scipen = 999)

 result.files = list.files(pattern="IFAA_simu_sample_*", full.names = F)
 to.process = length(result.files)
 
 output.file = paste(result.files[1], sep="")  
 result.i=get(load(output.file))
 
 print("result.files:")
 print(result.files)
 
 Regularization=result.i$analysisResults
 
 allData=data.matrix(read.csv(file=dataFile,header=T,na.strings=c(NA,"")))
 MVarNamLength=nchar(Mprefix)
 # get taxa variable names
 micros = sapply(substr(colnames(allData),1,MVarNamLength), function(x) {grep(Mprefix, x)})
 microPositions=which(micros == 1)
 rm(micros)
 allTaxaNames=colnames(allData)[microPositions]
 rm(microPositions)
 rm(allData,allTaxaNames)
 
 finalPenalMethod=Regularization$reguMethod
 nTaxa=Regularization$nTaxa
 nRef=result.i$nRef
 nPermu=result.i$nPermu
 nPredics=Regularization$nPredics
 truBetaMat=getTruBeta()[-1,] # without intercept
 biasCalIndex=which(truBetaMat[1,,drop=F]!=0,arr.ind=T)
 
 #biasCalIndex=which(truBetaMat[1,,drop=F]>0.11,arr.ind=T) # for violated assumption
 print("Percentage of weak associations:")
 print(1-mean(truBetaMat[1,,drop=F]>0.11))
 
 print("biasCalIndex:")
 print(biasCalIndex)
 
 truModelSizFWER=sum(colSums(truBetaMat[1,,drop=F]!=0)>0)
 truModelSizDect=sum(truBetaMat[1,]!=0)
 
 bootB=result.i$bootB
 bootLassoAlpha=result.i$bootLassoAlpha

 seqDepthDeno=getContinuSeqDepDeno()
 
 rep.no=NA
 exec.time=NA
 AnalyzedSubs=NA
 nTaxa=NA
 nTaxaTru=NA
 nPredics=NA
 nPredTru=NA
 truModelSize=NA
 finalizedRefTaxon=NA
 final2GoodRefTaxa=vector()
 truGoodIndpRefTax=NA
 finalizedBootRefTaxon=NA
 covForTruNon0=NA
 truIndpRefTax=NA
 bootLassoAlpha=NA
 sdTaxa=NA
 non0perct=NA
 finalBetaHat=NA
 finalRefTaxonQualified=NA
 bias=NA
 unEstiTax=list()
 
 for (i in 1:to.process) {	
   output.file = paste(result.files[i], sep="")  

   #message("Loading saved results from ", data.file, ", for simulation ", sim.number)
   result.i=get(load(output.file))
   Regularization=result.i$analysisResults
   totTime.i=result.i$totalTimeMins
   phase1time.i=Regularization$MCPExecuTime
   phase2time.i=totTime.i-phase1time.i
   dataSpars.i=result.i$dataSpars
   x1permut.i=result.i$x1permut
   AnalyzedSubs.i=Regularization$nSub
   nTaxa.i=Regularization$nTaxa
   nTaxaTru.i=floor(getNTaxa()*getEchoPropDiff())
   nPredics.i=Regularization$nPredics
   nPredTru.i=getNPredTru()
   
   truModelSize.i=nTaxaTru.i
   finalPenalMethod.i=Regularization$reguMethod
   finalizedRefTaxon.i.name=Regularization$finalizedBootRefTaxon
   namLength.i=nchar(finalizedRefTaxon.i.name)
   namPrefLength=nchar(Mprefix)
   
   numPart=c(namPrefLength+1,namLength.i)
   finalizedRefTaxon.i=as.numeric(substr(finalizedRefTaxon.i.name,numPart[1],numPart[2]))
   truIndpRefTax=which(colSums(truBetaMat[1,,drop=F])==0)
   truIndpRefTax.yesNo.i=finalizedRefTaxon.i%in%truIndpRefTax
   bootLassoAlpha.i=result.i$bootLassoAlpha

   beta.i=Regularization$betaMat
   CILow.i=Regularization$CILowMat
   CIUp.i=Regularization$CIUpMat

   beta.LPR.i=Regularization$betaMat.LPR
   CILow.LPR.i=Regularization$CILowMat.LPR
   CIUp.LPR.i=Regularization$CIUpMat.LPR
   
   beta10sum=matrix(0,nrow=nPredics.i,ncol=nTaxa.i)
   CILow10sum=matrix(0,nrow=nPredics.i,ncol=nTaxa.i)
   CIUp10sum=matrix(0,nrow=nPredics.i,ncol=nTaxa.i)
   beta10Prod=matrix(1,nrow=nPredics.i,ncol=nTaxa.i)

   beta10sum.LPR=matrix(0,nrow=nPredics.i,ncol=nTaxa.i)
   CILow10sum.LPR=matrix(0,nrow=nPredics.i,ncol=nTaxa.i)
   CIUp10sum.LPR=matrix(0,nrow=nPredics.i,ncol=nTaxa.i)
   beta10Prod.LPR=matrix(1,nrow=nPredics.i,ncol=nTaxa.i)

   nRefUsedForEst.i=Regularization$nRefUsedForEsti
   totRefTaxForEsti.i=length(Regularization$estiList)
   taxaNam10=names(Regularization$estiList)
   
   finalRefTaxa10.i=as.numeric(substr(taxaNam10,numPart[1],numPart[2]))
   
   truRefTaxa10.i=sum(finalRefTaxa10.i%in%truIndpRefTax)
   
   vec1=rep(NA,(nPredics.i*nTaxa.i))
   vec2=rep(NA,(nPredics.i*nTaxa.i))
   
   for (iii in 1:totRefTaxForEsti.i){
     vec1=rbind(vec1,as.vector(Regularization$estiList[[iii]]$finalBetaEst))
     vec2=rbind(vec2,as.vector(Regularization$estiList[[iii]]$finalBetaEst.LPR))
     }
   rm(iii)
   vec1=vec1[-1,]
   vec2=vec2[-1,]
   
   if(totRefTaxForEsti.i>1){
    beta10Stabi=matrix((diag(var(vec1))),nrow=nPredics.i)
    beta10Stabi.LPR=matrix((diag(var(vec2))),nrow=nPredics.i)
    }else{
     beta10Stabi=matrix(0,nrow=nPredics.i,ncol=nTaxa.i)
     beta10Stabi.LPR=matrix(0,nrow=nPredics.i,ncol=nTaxa.i)
     }
   
   for (ii in 1:totRefTaxForEsti.i){
    beta10sum=beta10sum+matrix(Regularization$estiList[[ii]]$finalBetaEst,nrow=nPredics.i)
    CILow10sum=CILow10sum+matrix(Regularization$estiList[[ii]]$CIvecLow,nrow=nPredics.i)
    CIUp10sum=CIUp10sum+matrix(Regularization$estiList[[ii]]$CIvecUp,nrow=nPredics.i)
    beta10Prod=beta10Prod*((matrix(Regularization$estiList[[ii]]$finalBetaEst,nrow=nPredics.i))>=0)
    
    beta10sum.LPR=beta10sum.LPR+matrix(Regularization$estiList[[ii]]$finalBetaEst.LPR,nrow=nPredics.i)
    CILow10sum.LPR=CILow10sum.LPR+matrix(Regularization$estiList[[ii]]$CIvecLow.LPR,nrow=nPredics.i)
    CIUp10sum.LPR=CIUp10sum.LPR+matrix(Regularization$estiList[[ii]]$CIvecUp.LPR,nrow=nPredics.i)
    beta10Prod.LPR=beta10Prod*((matrix(Regularization$estiList[[ii]]$finalBetaEst.LPR,nrow=nPredics.i))>=0)
    }
   beta10Ave.i=beta10sum/totRefTaxForEsti.i
   CILow10Ave.i=CILow10sum/totRefTaxForEsti.i
   CIUp10Ave.i=CIUp10sum/totRefTaxForEsti.i
   
   beta10Ave.LPR.i=beta10sum.LPR/totRefTaxForEsti.i
   CILow10Ave.LPR.i=CILow10sum.LPR/totRefTaxForEsti.i
   CIUp10Ave.LPR.i=CIUp10sum.LPR/totRefTaxForEsti.i

   biasBeta.i=(beta.i-truBetaMat)[biasCalIndex]
   biasBeta.LPR.i=(beta.LPR.i-truBetaMat)[biasCalIndex]
   biasBeta10Ave.i=(beta10Ave.i-truBetaMat)[biasCalIndex]
   biasBeta10Ave.LPR.i=(beta10Ave.LPR.i-truBetaMat)[biasCalIndex]
   beta10Prod.i=beta10Prod[biasCalIndex]
   beta10Prod.LPR.i=beta10Prod.LPR[biasCalIndex]

   beta10Stabi.i=beta10Stabi[biasCalIndex]
   beta10Stabi.LPR.i=beta10Stabi.LPR[biasCalIndex]
   
   print("i:")
   print(i)
   print("beta10Prod.LPR.i:")
   print(beta10Prod.LPR.i)
   print("beta10Prod.i:")
   print(beta10Prod.i)
   
   cpBeta.i=(truBetaMat>=CILow.i) & (truBetaMat<=CIUp.i)
   cpBeta.LPR.i=(truBetaMat>=CILow.LPR.i) & (truBetaMat<=CIUp.LPR.i)
   cpBeta10Ave.i=(truBetaMat>=CILow10Ave.i) & (truBetaMat<=CIUp10Ave.i)
   cpBeta10Ave.LPR.i=(truBetaMat>=CILow10Ave.LPR.i) & (truBetaMat<=CIUp10Ave.LPR.i)
   
   if(i==1){
     totTime=totTime.i
     phase1time=phase1time.i
     phase2time=phase2time.i
     AnalyzedSubs=AnalyzedSubs.i
     nTaxa=nTaxa.i
     nTaxaTru=nTaxaTru.i
     nPredics=nPredics.i
     nPredTru=nPredTru.i
     dataSpars=dataSpars.i
     x1permut=x1permut.i
     truModelSize=truModelSize.i
     finalPenalMethod=finalPenalMethod.i
     truIndpRefTax.yesNo=truIndpRefTax.yesNo.i
     truRefTaxa10=truRefTaxa10.i
     nRefUsedForEst=nRefUsedForEst.i
     finalizedRefTaxon=finalizedRefTaxon.i
     finalRefTaxa10=finalRefTaxa10.i
     bootLassoAlpha=bootLassoAlpha.i
     totRefTaxForEsti=totRefTaxForEsti.i
     biasBeta=biasBeta.i
     biasBeta.LPR=biasBeta.LPR.i
     biasBeta10Ave=biasBeta10Ave.i
     biasBeta10Ave.LPR=biasBeta10Ave.LPR.i
     beta10signChk=beta10Prod.i
     beta10signChk.LPR=beta10Prod.LPR.i
     beta10Var=beta10Stabi.i
     beta10Var.LPR=beta10Stabi.LPR.i
     
     cpBeta=cpBeta.i
     cpBeta.LPR=cpBeta.LPR.i
     cpBeta10Ave=cpBeta10Ave.i
     cpBeta10Ave.LPR=cpBeta10Ave.LPR.i
     
      } else {
        totTime= rbind(totTime,totTime.i)
        phase1time= rbind(phase1time,phase1time.i)
        phase2time= rbind(phase2time,phase2time.i)
        AnalyzedSubs= rbind(AnalyzedSubs,AnalyzedSubs.i)
        nTaxa= rbind(nTaxa,nTaxa.i)
        nTaxaTru= rbind(nTaxaTru,nTaxaTru.i)
        nPredics= rbind(nPredics,nPredics.i)
        nPredTru= rbind(nPredTru,nPredTru.i)
        truModelSize= rbind(truModelSize,truModelSize.i)
        dataSpars= rbind(dataSpars,dataSpars.i)
        x1permut= rbind(x1permut,x1permut.i)
        finalPenalMethod= rbind(finalPenalMethod,finalPenalMethod.i)
        truIndpRefTax.yesNo= rbind(truIndpRefTax.yesNo,truIndpRefTax.yesNo.i)
        truRefTaxa10= rbind(truRefTaxa10,truRefTaxa10.i)
        finalRefTaxa10= rbind(finalRefTaxa10,finalRefTaxa10.i)
        nRefUsedForEst= rbind(nRefUsedForEst,nRefUsedForEst.i)
        finalizedRefTaxon= rbind(finalizedRefTaxon,finalizedRefTaxon.i)
        bootLassoAlpha= rbind(bootLassoAlpha,bootLassoAlpha.i)
        totRefTaxForEsti=rbind(totRefTaxForEsti,totRefTaxForEsti.i)
        biasBeta=biasBeta+biasBeta.i
        biasBeta.LPR=biasBeta.LPR+biasBeta.LPR.i
        biasBeta10Ave=biasBeta10Ave+biasBeta10Ave.i
        biasBeta10Ave.LPR=biasBeta10Ave.LPR+biasBeta10Ave.LPR.i
        beta10signChk=beta10signChk+beta10Prod.i
        beta10signChk.LPR=beta10signChk.LPR+beta10Prod.LPR.i
        beta10Var=beta10Var+beta10Stabi.i
        beta10Var.LPR=beta10Var.LPR+beta10Stabi.LPR.i
        
        cpBeta=cpBeta+cpBeta.i
        cpBeta.LPR=cpBeta.LPR+cpBeta.LPR.i
        cpBeta10Ave=cpBeta10Ave+cpBeta10Ave.i
        cpBeta10Ave.LPR=cpBeta10Ave.LPR+cpBeta10Ave.LPR.i
        }
    }
 
 biasBeta=biasBeta/to.process
 biasBeta.LPR=biasBeta.LPR/to.process
 biasBeta10Ave=biasBeta10Ave/to.process
 biasBeta10Ave.LPR=biasBeta10Ave.LPR/to.process
 beta10signChk=beta10signChk/to.process
 beta10signChk.LPR=beta10signChk.LPR/to.process
 beta10Var=beta10Var/to.process
 beta10Var.LPR=beta10Var.LPR/to.process

 cpBeta=cpBeta[biasCalIndex]/to.process
 cpBeta.LPR=cpBeta.LPR[biasCalIndex]/to.process
 cpBeta10Ave=cpBeta10Ave[biasCalIndex]/to.process
 cpBeta10Ave.LPR=cpBeta10Ave.LPR[biasCalIndex]/to.process

 rownames(AnalyzedSubs) <- NULL
 colnames(AnalyzedSubs) <- c("nSub")

 rownames(nTaxa) <- NULL
 colnames(nTaxa) <- c("nTaxa")

 rownames(x1permut) <- NULL
 colnames(x1permut) <- c("x1permut")
 
 rownames(nTaxaTru) <- NULL
 colnames(nTaxaTru) <- c("nTaxaTru")

 rownames(nPredics) <- NULL
 colnames(nPredics) <- c("nPredics")

 rownames(nPredTru) <- NULL
 colnames(nPredTru) <- c("nPredTru")

 rownames(finalizedRefTaxon) <- NULL
 colnames(finalizedRefTaxon) <- c("finalizedRefTaxon")

 rownames(truIndpRefTax.yesNo) <- NULL
 colnames(truIndpRefTax.yesNo) <- c("truIndpRefTax.yesNo")

 rownames(truRefTaxa10) <- NULL
 colnames(truRefTaxa10) <- c("truRefTaxa10")
 
 rownames(nRefUsedForEst) <- NULL
 colnames(nRefUsedForEst) <- c("nRefUsedForEst")
 
 meanAbsBiasBeta=mean(abs(biasBeta))
 meanAbsBiasBeta.LPR=mean(abs(biasBeta.LPR))
 print("abs(biasBeta.LPR):")
 print(abs(biasBeta.LPR))
 
 meanAbsBiasBeta10Ave=mean(abs(biasBeta10Ave))
 meanAbsBiasBeta10Ave.LPR=mean(abs(biasBeta10Ave.LPR))

 meanCPBeta=mean(cpBeta)
 meanCPBeta.LPR=mean(cpBeta.LPR)
 meanCPBeta10Ave=mean(cpBeta10Ave)
 meanCPBeta10Ave.LPR=mean(cpBeta10Ave.LPR)

 print("meanAbsBiasBeta:")
 print(meanAbsBiasBeta)
 print("meanAbsBiasBeta.LPR:")
 print(meanAbsBiasBeta.LPR)
 print("meanAbsBiasBeta10Ave:")
 print(meanAbsBiasBeta10Ave)
 print("meanAbsBiasBeta10Ave.LPR:")
 print(meanAbsBiasBeta10Ave.LPR)

 print("meanCPBeta:")
 print(meanCPBeta)
 print("meanCPBeta.LPR:")
 print(meanCPBeta.LPR)
 print("meanCPBeta10Ave:")
 print(meanCPBeta10Ave)
 print("meanCPBeta10Ave.LPR:")
 print(meanCPBeta10Ave.LPR)
 
 print("meanAbsBeta:")
 print((abs(truBetaMat[biasCalIndex])))
 print(mean(abs(truBetaMat[biasCalIndex])))
 print("meanAbsBiasBeta percentage:")
 print(100*(abs(biasBeta/(truBetaMat[biasCalIndex]))))
 print(100*mean(abs(biasBeta/(truBetaMat[biasCalIndex]))))
 print("meanAbsBiasBeta.LPR percentage:")
 print(100*(abs(biasBeta.LPR/(truBetaMat[biasCalIndex]))))
 print(100*mean(abs(biasBeta.LPR/(truBetaMat[biasCalIndex]))))
 print("meanAbsBiasBeta10Ave percentage:")
 print(100*(abs(biasBeta10Ave/(truBetaMat[biasCalIndex]))))
 print(100*mean(abs(biasBeta10Ave/(truBetaMat[biasCalIndex]))))
 print("meanAbsBiasBeta10Ave.LPR percentage:")
 print(100*(abs(biasBeta10Ave.LPR/(truBetaMat[biasCalIndex]))))
 print(100*mean(abs(biasBeta10Ave.LPR/(truBetaMat[biasCalIndex]))))
 
 dataFeatures=cbind(AnalyzedSubs,nTaxa,nTaxaTru,nPredics,
        nPredTru,finalizedRefTaxon,truIndpRefTax.yesNo,truRefTaxa10,
        nRefUsedForEst,
        truModelSizFWER,truModelSizDect,dataSpars,x1permut)

 meanAnalyzedSubs=mean(as.numeric(AnalyzedSubs))
 cat("meanAnalyzedSubs:",meanAnalyzedSubs,"\n")

 print(paste(" ")) 
 print(paste("Number of subjects generated in each data set: ", AnalyzedSubs[1]))
 print(paste("Average analyzed number of subjects per data set: ", meanAnalyzedSubs))
 print(paste("Number of taxa: ", nTaxa[1]))
 print(paste("Number of predictors: ", nPredics[1]))
 print(paste("Number of predictors with non-zero coeffi: ", nPredTru[1])) 
 print(paste("Bootstrap samples: ", bootB))    
 print(paste("Number of analyzed data sets: ", to.process))    
 print(paste("Sequencing depth: ", seqDepthDeno)) 
 
 print("dataFeatures[,]:")
 print(dataFeatures[,])
 
 colnames(truIndpRefTax.yesNo) <- c("truIndpRefTax.yesNo")
 cat("Mean truIndp.yesno",mean(truIndpRefTax.yesNo),"\n")
 
 rownames(truRefTaxa10) <- NULL
 colnames(truRefTaxa10) <- c("truRefTaxa10")
 cat("Mean truRefTaxa10",mean(truRefTaxa10),"\n")
 
 fwerInd=fwerIndex(nDataSet=to.process,
                   truBeta=truBetaMat,
                   fwer=fwer)
 print(fwerInd)

 cat("dataSpars:",mean(dataSpars),"\n")
 
 print("Proportion of no sign change:")
 print(sum(beta10signChk)/length(beta10signChk))
 print("beta10signChk:")
 print(beta10signChk)
 
 print("Proportion of no sign change with LPR:")
 print(sum(beta10signChk.LPR)/length(beta10signChk.LPR))
 print("beta10signChk.LPR:")
 print(beta10signChk.LPR)
 
 print("mean variance of 10 beta:")
 print(mean(beta10Var))
 
 print("mean variance of 10 beta.LPR:")
 print(mean(beta10Var.LPR))

 cat("min totTime:",min(totTime),"minutes.","\n",
     "Mean totTime:", mean(totTime),"minutes.","\n",
     "max totTime:",max(totTime),"minutes.","\n")
 
 cat("min phase 1 time:",min(phase1time),"minutes.","\n",
     "Mean phase 1 time:", mean(phase1time),"minutes.","\n",
     "max phase 1 time:",max(phase1time),"minutes.","\n")
 
 cat("min phase 2 time:",min(phase2time),"minutes.","\n",
     "Mean phase 2 time:", mean(phase2time),"minutes.","\n",
     "max phase 2 time:",max(phase2time),"minutes.","\n")

 print(paste(" ")) 
}


#-------------------------------------------------------------------
## Function to calculate recall, precision and F1
#-------------------------------------------------------------------

PerformInd=function(coefEstVec,truVec){
  results=list()
  
  if(length(coefEstVec)!=length(truVec)){
    stop("Performance indices cannot be calculated due to different vector dimensions")
  }
  
  TruPos=sum((coefEstVec*truVec)>0)
  TruNeg=sum((coefEstVec==0) & (truVec==0))
  FalsePos=sum(((coefEstVec!=0) & (truVec==0))|((coefEstVec*truVec)<0))
  FalseNeg=sum((coefEstVec==0) & (truVec!=0))
  modelSize=sum(coefEstVec!=0)
  rm(coefEstVec,truVec)
  
  recall=TruPos/(TruPos+FalseNeg)
  precision=TruPos/(TruPos+FalsePos)
  F1=2*(recall*precision)/(recall+precision)
  typeIerror=FalsePos/(FalsePos+TruNeg)
  
  if(modelSize==0){
    recall=0
    precision=1
    F1=0
    typeIerror=0
  }
  
  results$recall=recall
  results$precision=precision
  results$F1=F1
  results$typeI=typeIerror
  results$modelSize=modelSize
  
  return(results)
}


#
## function to load family wise error rate (FWER) control performance indices
#

fwerIndex=function(nDataSet,
                   truBeta,
                   fwer
){
  result.files = list.files(pattern="IFAA_simu_sample_*", full.names = F)
  recall=NA
  precision=NA
  F1=NA
  beta=NA
  modelSize=NA
  exeTime=NA
  truPredicsCorr=NA
  covForTruNon0=NA
  rep.no=NA
  for (i in 1:nDataSet) {	
    output.file.i = paste(result.files[i], sep="")  
    result.i=get(load(output.file.i))

    Regularization=result.i$analysisResults
      
    if(length(fwer)>0){
      cut=quantile(Regularization$maxVec,(1-fwer))
      fwerSelect=(Regularization$selecCountOverall>=cut)+0

      indicesFWER=PerformInd(fwerSelect,((colSums(truBeta[1,,drop=F]!=0))+0))
      
      # for violated assumption
      #indicesFWER=PerformInd(fwerSelect,((colSums(truBeta[1,,drop=F]>0.05))+0))

      recall.i=t(indicesFWER$recall)
      precision.i=t(indicesFWER$precision)
      F1.i=t(indicesFWER$F1)
      typeI.i=t(indicesFWER$typeI)
      modelSize.i=indicesFWER$modelSize
      
      if(length(result.i$testCov)>1){
       X1SelectCut=quantile(Regularization$MaxMatTestCovByPermu[4,],probs=(1-fwer))
       X1Select=(Regularization$selecCountMatIndv[4,]>=X1SelectCut)+0
       cat("(X1SelectCut):",X1SelectCut,"\n")
       cat("length(X1Select):",length(X1Select),"\n")
       cat("Regularization$selecCountMatIndv:",Regularization$selecCountMatIndv,"\n")
       cat("Regularization$MaxMatTestCovByPermu:",Regularization$MaxMatTestCovByPermu,"\n")
      
       cat("length((truBeta[4,]!=0)+0):",length((truBeta[4,]!=0)+0),"\n")
      
       detecIndicesFWER=PerformInd(X1Select,((truBeta[4,]!=0)+0))
       recallDect.i=t(detecIndicesFWER$recall)
       precisionDect.i=t(detecIndicesFWER$precision)
       F1Dect.i=t(detecIndicesFWER$F1)
       typeIDect.i=t(detecIndicesFWER$typeI)
       modelSizeDect.i=detecIndicesFWER$modelSize
       X1SelCut.i=X1SelectCut
       }
      fwerRate.i=fwer
      fwerCut.i=cut
    }
    
    if(i==1){
      output.file=output.file.i
      #rep.no=rep.no.i
      recall=recall.i
      precision=precision.i
      F1=F1.i
      typeI=typeI.i
      if(length(result.i$testCov)>1){
        recallDect=recallDect.i
        precisionDect=precisionDect.i
        F1Dect=F1Dect.i
        typeIDect=typeIDect.i
        modelSizeDect=modelSizeDect.i
        X1SelCut=X1SelCut.i
        }
      modelSize=modelSize.i
      fwerCut=fwerCut.i
      fwerRate=fwerRate.i
    }else {
      
      output.file=rbind(output.file,output.file.i)
      #rep.no=rbind(rep.no,rep.no.i)
      recall=rbind(recall,recall.i)
      precision=rbind(precision,precision.i)
      if(length(result.i$testCov)>1){
        recallDect=rbind(recallDect,recallDect.i)
        precisionDect=rbind(precisionDect,precisionDect.i)
        F1Dect=rbind(F1Dect,F1Dect.i)
        typeIDect=rbind(typeIDect,typeIDect.i)
        modelSizeDect=rbind(modelSizeDect,modelSizeDect.i)
        X1SelCut=rbind(X1SelCut,X1SelCut.i)
        }
      F1=rbind(F1,F1.i)
      typeI=rbind(typeI,typeI.i)
      modelSize=rbind(modelSize,modelSize.i)
      fwerCut=rbind(fwerCut,fwerCut.i)
      fwerRate=rbind(fwerRate,fwerRate.i)
    }
  }
  
  if(length(result.i$testCov)>1){
   results=cbind(recall,recallDect,precision,precisionDect,
                F1,F1Dect,typeI,typeIDect,modelSize,modelSizeDect,
                fwerCut,X1SelCut,
                fwerRate)
   
   varNames=c("recallFWER","recallDect","precFWER","precisionDect",
             "F1FWER","F1Dect","typeIFWER","typeIDect",
             "modelSizeFWER","modelSizeDect","fwerCut",
             "X1SelCut","fwerRate"
             )
  }
  
  if(length(result.i$testCov)==1){
    results=cbind(recall,precision,
                  F1,typeI,modelSize,
                  fwerCut,
                  fwerRate)
    
    varNames=c("recallFWER","precFWER",
               "F1FWER","typeIFWER",
               "modelSizeFWER","fwerCut",
               "fwerRate")
   }
  colnames(results)=varNames
  rownames(results) <- NULL
  
  nUsed=nrow(na.omit(results))
  
  results=as.data.frame(results)
  cat("length(output.file):",length(output.file),"\n")
  cat("nrow(results):",nrow(results),"\n")
  
  print(cbind(output.file,results))
  
  meanResults=t(as.matrix(apply(results,2,mean)))
  rownames(meanResults) <- "mean"

  # return results
  results=list()
  results$meanResults=meanResults
  return(results)
}

#
## run summary
#
analyzeResults()
