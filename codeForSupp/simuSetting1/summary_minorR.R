#
### analyze the saved output data sets to generate the results
#

setwd("C:/Dropbox (UFL)/papers/IFAA/R")

source("../generatData/TruParam_minorR.R")
source("utils_minorR.R")

analyzeResults=function(
    fwer
){

 options(digits=3)
 options(scipen = 999)

 result.files = list.files(pattern="simuResu_*", full.names = F)
 to.process = length(result.files)
 
 output.file = paste(result.files[1], sep="")  
 result.i=load(output.file)
 
 print("result.files:")
 print(result.files)
 
 meanLogDiffAncom=getTruLogDiffAnc()
 
 finalBetaHat=NA
 bias=NA

 for (i in 1:to.process) {	
   output.file = paste(result.files[i], sep="")  
   result.i=load(output.file)
   
   finalBetaHat.i=as.vector(Regularization$estiResults$finalBetaEst)
   biasCalIndex.i=(meanLogDiffAncom!=0)
   bias.i=t(finalBetaHat.i-meanLogDiffAncom)[biasCalIndex.i]

   if(i==1){
      finalBetaHat=finalBetaHat.i
      bias=bias.i
      } else {
          finalBetaHat=rbind(finalBetaHat,finalBetaHat.i)
          bias=rbind(bias,bias.i)
          }
   }

 biases=colMeans(bias)
 meanAbsBias=mean(abs(biases))

 meanDiffTru=mean(meanLogDiffAncom[!(meanLogDiffAncom==0)])
 print("meanDiffTru:")
 print(meanDiffTru)
 print("meanAbsBias:")
 print(meanAbsBias)
 print("meanAbsBias percentage:")
 print(100*meanAbsBias/meanDiffTru)
 
 fwerInd=fwerIndex(nDataSet=to.process,fwer=fwer)
 print(fwerInd)

 ancomInd=ancomIndex(nDataSet=to.process)
 print(ancomInd)

 zigInd=zigIndex(nDataSet=to.process)
 print(zigInd)
  
 edgeRInd=edgeRIndex(nDataSet=to.process)
 print(edgeRInd)
  
 DESeq2Ind=DESeqRIndex(nDataSet=to.process)
 print(DESeq2Ind)
 }


#
## function to load family wise error control performance indices
#

fwerIndex=function(nDataSet,fwer
  ){
   result.files = list.files(pattern="simuResu_*", full.names = F)
   recall=NA
   precision=NA
   F1=NA
   typeI=NA
   modelSize=NA

 for (i in 1:nDataSet) {	
   output.file = paste(result.files[i], sep="")  
   result.i=load(output.file)
 
   recall.i=t(Regularization$recallFWER)
   precision.i=t(Regularization$precisionFWER)
   F1.i=t(Regularization$F1FWER)
   typeI.i=t(Regularization$typeIFWER)
   modelSize.i=Regularization$modelSizeFWER

   if(i==1){
      recall=recall.i
      precision=precision.i
      F1=F1.i
      typeI=typeI.i
      modelSize=modelSize.i
      }else {
          recall=rbind(recall,recall.i)
          precision=rbind(precision,precision.i)
          F1=rbind(F1,F1.i)
          typeI=rbind(typeI,typeI.i)
          modelSize=rbind(modelSize,modelSize.i)
         }
   }

  results=cbind(recall,precision,,
     F1,typeI,,modelSize)

  varNames=c("recallFWER","precFWER",
   "F1FWER","typeIFWER",
      "modelSizeFWER")

  colnames(results)=varNames
  rownames(results) <- NULL
  
  nUsed=nrow(na.omit(results))

  results=as.data.frame(results)
  print(results)

  meanResults=t(as.matrix(apply(results,2,mean)))
  rownames(meanResults) <- "mean"

  # return results
  results=list()
  results$meanResults=meanResults
  return(results)
 }


#
## function to load ZIG performance indices
#

zigIndex=function(nDataSet=3
  ){
   result.files = list.files(pattern="simuResu_*", full.names = F)
   recall=NA
   precision=NA
   F1=NA
   beta=NA
   modelSize=NA
   exeTime=NA
   truPredicsCorr=NA

 for (i in 1:nDataSet) {	
   output.file = paste(result.files[i], sep="")  

   result.i=load(output.file)

   recall.i=t(Regularization$recallZig)
   precision.i=t(Regularization$precisionZig)
   F1.i=t(Regularization$F1Zig)
   typeI.i=t(Regularization$typeIZig)
   modelSize.i=Regularization$modelSizeZig
   if(i==1){
      recall=recall.i
      precision=precision.i
      F1=F1.i
      typeI=typeI.i
      modelSize=modelSize.i
      }else {
          recall= rbind(recall,recall.i)
          precision= rbind(precision,precision.i)
          F1= rbind(F1,F1.i)
          typeI= rbind(typeI,typeI.i)
          modelSize= rbind(modelSize,modelSize.i)

         }
   }

  results=cbind(recall,precision,F1,typeI,modelSize)
  varNames=c("recallZIG","precZIG","F1ZIG","typeIZIG",
      "modelSizeZIG")

  colnames(results)=varNames
  rownames(results) <- NULL
  
  results=as.data.frame(results)
  print(results)

  meanResults=t(as.matrix(apply(results,2,mean)))
  rownames(meanResults) <- "mean"

  # return results
  results=list()
  results$meanResults=meanResults
  return(results)
 }

#
## function to load edgeR performance indices
#

edgeRIndex=function(nDataSet=3
  ){
   result.files = list.files(pattern="simuResu_*", full.names = F)
   recall=NA
   precision=NA
   F1=NA
   beta=NA
   modelSize=NA
   exeTime=NA
   truPredicsCorr=NA

 for (i in 1:nDataSet) {	
   output.file = paste(result.files[i], sep="")  

   #message("Loading saved results from ", data.file, ", for simulation ", sim.number)
   result.i=load(output.file)

   recall.i=t(Regularization$recallEdgeR)
   precision.i=t(Regularization$precisionEdgeR)
   F1.i=t(Regularization$F1EdgeR)
   typeI.i=t(Regularization$typeIEdgeR)
   modelSize.i=Regularization$modelSizeEdgeR
   if(i==1){
      recall=recall.i
      precision=precision.i
      F1=F1.i
      typeI=typeI.i
      modelSize=modelSize.i
      }else {
          recall= rbind(recall,recall.i)
          precision= rbind(precision,precision.i)
          F1= rbind(F1,F1.i)
          typeI= rbind(typeI,typeI.i)
          modelSize= rbind(modelSize,modelSize.i)
         }
   }

  results=cbind(recall,precision,F1,typeI,modelSize,seqDepthDeno)
  varNames=c("recallEdgeR","precEdgeR","F1EdgeR","typeIEdgeR",
     "modelSizeEdgeR")

  colnames(results)=varNames
  rownames(results) <- NULL
  
  results=as.data.frame(results)
  print(results)

  meanResults=t(as.matrix(apply(results,2,mean)))
  rownames(meanResults) <- "mean"

  # return results
  results=list()
  results$meanResults=meanResults
  return(results)
 }

#
## function to load DESeq2 performance indices
#

DESeqRIndex=function(nDataSet
){
   result.files = list.files(pattern="simuResu_*", full.names = F)
   
   DESeq2Fail=NA
   recall=NA
   precision=NA
   F1=NA
   beta=NA
   modelSize=NA
   exeTime=NA
   truPredicsCorr=NA
   
   for (i in 1:nDataSet) {	
      output.file = paste(result.files[i], sep="")  
      
      #message("Loading saved results from ", data.file, ", for simulation ", sim.number)
      result.i=load(output.file)
      
      DESeq2Fail.i=t(Regularization$DESeqFail)
      recall.i=t(Regularization$recallDESeq2)
      precision.i=t(Regularization$precisionDESeq2)
      F1.i=t(Regularization$F1DESeq2)
      typeI.i=t(Regularization$typeIDESeq2)
      modelSize.i=Regularization$modelSizeDESeq2
      if(i==1){
         DESeq2Fail=DESeq2Fail.i
         recall=recall.i
         precision=precision.i
         F1=F1.i
         typeI=typeI.i
         modelSize=modelSize.i
      }else {
         DESeq2Fail= rbind(DESeq2Fail,DESeq2Fail.i)
         recall= rbind(recall,recall.i)
         precision= rbind(precision,precision.i)
         F1= rbind(F1,F1.i)
         typeI= rbind(typeI,typeI.i)
         modelSize= rbind(modelSize,modelSize.i)
      }
   }
   
   results=cbind(DESeq2Fail,recall,precision,F1,typeI,modelSize)
   varNames=c("DESeq2Fail","recallDESeq2","precDESeq2","F1DESeq2","typeIDESeq2",
              "modelSizeDESeq2")
   
   colnames(results)=varNames
   rownames(results) <- NULL
   
   results=as.data.frame(results)
   print(results)

   meanResults=t(as.matrix(apply(results,2,mean)))
   rownames(meanResults) <- "mean"

   # return results
   results=list()
   results$meanResults=meanResults
   return(results)
}

#
## function to load ANCOM performance indices
#

ancomIndex=function(nDataSet=3
  ){
   result.files = list.files(pattern="simuResu_*", full.names = F)
   recall=NA
   precision=NA
   F1=NA
   beta=NA
   modelSize=NA
   exeTime=NA
   truPredicsCorr=NA

 for (i in 1:nDataSet) {	
   output.file = paste(result.files[i], sep="")  

   #message("Loading saved results from ", data.file, ", for simulation ", sim.number)
   result.i=load(output.file)

   recall.i=t(Regularization$recallAncom)
   precision.i=t(Regularization$precisionAncom)
   F1.i=t(Regularization$F1Ancom)
   typeI.i=t(Regularization$typeIancom)
   modelSize.i=Regularization$modelSizeAncom
   if(i==1){
      recall=recall.i
      precision=precision.i
      F1=F1.i
      typeI=typeI.i
      modelSize=modelSize.i
      }else {
          recall= rbind(recall,recall.i)
          precision= rbind(precision,precision.i)
          F1= rbind(F1,F1.i)
          typeI= rbind(typeI,typeI.i)
          modelSize= rbind(modelSize,modelSize.i)
         }
   }

  results=cbind(recall,precision,F1,typeI,modelSize)
  varNames=c("recallANCOM","precANCOM","F1ANCOM","typeIANCOM",
     "modelSizeANCOM")

  colnames(results)=varNames
  rownames(results) <- NULL
  
  nUsed=nrow(na.omit(results))

  results=as.data.frame(results)
  print(results)

  meanResults=t(as.matrix(apply(results,2,mean)))
  rownames(meanResults) <- "mean"
  
  # return results
  results=list()
  results$meanResults=meanResults
  return(results)
 }

 analyzeResults(fwer=0.20)



