# setwd("C:/Dropbox (UFL)/papers/IFAA/Rpack/dataGeneratR")

source("TruParam_Revision_minor.R")

#
## define the distribution of associations
#
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

#
## generate the matrix of beta^k
#
generatBetaK=function(
  nTaxa=getNTaxa(),
  nPredics=getNPredics(),
  propDiff=getEchoPropDiff(),
  realDataParam=getRealDataParam(),
  CiforGeneratBeta0k=1/40,
  intercepMax=3,
  rareTaxaMeanAdj=2,
  assumption=F
  ){
  
  emptBetaMat=matrix(0,nrow=(nPredics+1),ncol=nTaxa)
  beta0kVec=pmin(realDataParam$MeanLogPosVec,intercepMax)-log(CiforGeneratBeta0k)
  beta0kVec[(realDataParam$rareTaxa)]=beta0kVec[(realDataParam$rareTaxa)]+rareTaxaMeanAdj
  emptBetaMat[1,]=beta0kVec
  
  # randomly generate truly associated taxa
  nTaxaTru=floor(nTaxa*propDiff) 
  truTaxaInd=sample(nTaxa,nTaxaTru)
  assocPar=assocDist(n=nTaxaTru)

  # randomly generate effect size 
  for (j in 1:nTaxaTru){
    iTaxon=truTaxaInd[j]
    beta.i.j=runif(n=1, min=assocPar[1,j], max=assocPar[2,j])
    emptBetaMat[2,iTaxon]=beta.i.j
    }
      
  if(!assumption){emptBetaMat[2,(emptBetaMat[2,]==0)]=0.05}
      
    for (i in 3:4) {
      for (j in 1:nTaxaTru){
        iTaxon=truTaxaInd[j]
        beta.i.j=-runif(n=1,0.1,0.5)
        emptBetaMat[i,iTaxon]=beta.i.j
        }
      }

    write.csv(emptBetaMat,file="TruBetaMat.csv",row.names=F)
    
 }

generatBetaK()

