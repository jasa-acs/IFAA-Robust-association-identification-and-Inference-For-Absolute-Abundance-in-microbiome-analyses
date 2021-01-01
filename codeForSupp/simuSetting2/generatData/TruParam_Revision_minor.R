# setwd("C:/Dropbox (UFL)/papers/IFAA/Rpack/dataGeneratR")

#
## true parameter values and settings
#

truParams=get(load("TruParam.RData"))
getNSub=function(){return(truParams$RealDataParam$nSubReal)}
getNdataPerBatch=function(){return(2)}
getTruPredicsSigma=function(){return(1)}
getTruPredicsRho=function(){return(0.4)}
getNPredics=function(){return(3)}
getNTaxa=function(){return(truParams$RealDataParam$nTaxaReal)}

getTruBeta=function(){
  return(truParams$truBeta)
  #return(truParams$truBetaVioAssump)
  }

getRealDataParam=function(){return(truParams$RealDataParam)}

getEchoAbundanDist=function(){return(c(0.1,0.3,0.6))}

getEchoPropDiff=function(){
    return(0.25)
   }

getEchoEffectsDist=function(){return(c(1,2,3))}

getEchoEffects=function(){
  return(matrix(c(4,4,2.8,2.8,2,2),nrow=2))
}

getContinuSeqDepDeno=function(){return(c(20,80))}


