
#
## this is an hypothetical example of using the IFAA() function, the data
## set dataExampleForPublish.csv is a randomly generated with a similar abundance level
## as the NHBCS data
#

library(IFAA)

data=read.csv("dataExampleForPublish.csv",header=T)
microbMat=data[,c("id",paste0("rawCount",seq(218)))]

x=data[,c("id","x1","x2","x3")]


IFAAresults=IFAA(MicrobData=microbMat,CovData=x,linkIDname="id",
      testCov=c("x1"),ctrlCov=c("x2","x2"),
      nRef=4,
      nPermu=4,
      paraJobs=4,
      fwerRate=0.25,
      bootB=4)

IFAAresults$analysisResults$estByCovList
