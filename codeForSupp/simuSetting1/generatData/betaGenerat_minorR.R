# setwd("C:/Dropbox (UFL)/papers/IFAA/dataGeneratR")

source("TruParam_minorR.R")

# define the taxa abundance distribution interm high,medium and low abundance
abundanDist=function(
   n,
   par=getEchoAbundanDist(),
   gammaPar=getEchoGammaPar()
   ){
    return(sample(x=gammaPar,size=n,prob=par,replace=T))
    }
 
# abundanDist(40)

# generate the baseline taxa counts in the echo system
ecoSysBasePar=function(
    nTaxa=getNTaxa()
   ){
    taxaHMLGroup=abundanDist(n=nTaxa)
    basePoiPar=rgamma(n=nTaxa,shape=taxaHMLGroup)
    return(basePoiPar)
    }

# ecoSysBasePar()

# define the distribution of associations
assocDist=function(
   n,
   par=getEchoAbundanDist(),# same distribution as the abundance distribution
   echoEffe=getEchoEffects(),
   echoEffeDist=getEchoEffectsDist()
   ){
    effeCol=sample(x=echoEffeDist,size=n,prob=par,replace=T)
    return(echoEffe[,effeCol])
    }
 
# assocDist(n=4)

# geneate the association matrix between sparsity and covariates
generatEchoTruAssoc=function(
  nTaxa=getNTaxa(),
  propDiff=getEchoPropDiff(),
  nPredics=getNPredics(),
  nPredTru=getNPredTru()
 ){
  # load baseline spars
  loadBasePar=ecoSysBasePar(nTaxa=nTaxa)
  alphaInter=loadBasePar

  # start to generate full regression coefficients matrix
  nColAlpha=nPredics+1
  alphaM=matrix(0,nrow=nTaxa,ncol=nColAlpha)
  alphaM[,1]=alphaInter

  # randomly generate true associated predictors and taxa
  nTaxaTruAncom=floor(nTaxa*propDiff)

 if(nTaxaTruAncom>0){
  truPredInd=sample(nPredics,nPredTru)+1
  truTaxaInd=matrix(NA,ncol=nPredTru,nrow=nTaxaTruAncom)

  for(i in 1:nPredTru){
       truTaxaInd[,i]=sample(nTaxa,nTaxaTruAncom)
       }

  # randomly generate effect size 
  for (j in truPredInd) {
    assocPar=assocDist(n=nTaxaTruAncom)
      for (i in 1:nTaxaTruAncom){
        iTaxon=truTaxaInd[i,which(truPredInd==j)]
        posUni=runif(n=1, min=assocPar[1,i], max=assocPar[2,i])
        alphaM[iTaxon,j]=runif(n=1, min=assocPar[1,i], max=assocPar[2,i])
        }
      }
  }
  return(alphaM)
 }


 TruBetaAncom=generatEchoTruAssoc()
 write.csv(TruBetaAncom,file="TruBetaAncom.csv",row.names=F)
 
