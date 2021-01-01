
#
## generate true parameter values
#

getNSub=function(){return(50)}
getNdataPerBatch=function(){return(100)}
getNPredics=function(){return(1)}
getNPredTru=function(){return(1)}
getNTaxa=function(){return(500)}


getEchoAbundanDist=function(){return(c(0.1,0.3,0.6))}

getEchoGammaPar=function(){return(c(10000,200,50))}

getEchoPropDiff=function(){
    return(0.25)
   }

getEchoEffectsDist=function(){return(c(1,2,3))}

getEchoEffects=function(){
    return(matrix(c(10000,15000,200,400,100,150),nrow=2))
   }

# getEchoEffects()

# get the true parameters
truParames=get(load("truParams1.RData"))

getTruBetaEcho=function(){
  return(truParames$truBeta)
  }

# get the true treatment effect on log scale

getTruLogDiffAnc=function(){
  return(truParames$truDiffLog)
  }

truDiffLog=getTruLogDiffAnc()


