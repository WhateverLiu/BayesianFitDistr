






gm2dEM<-function(xCoor, xval, gmWeight, gmMean, gmVar, weightEPS=1e-20, convergeEPS=1e-3, maxit=500L, showProgress=5L, useRelativeDiff=1, ellipseAxisRatioThreshold=14.2){
gm2d(xCoor[[1]], xCoor[[2]], xval, gmWeight, gmMean[[1]], gmMean[[2]],
  gmVar[[1]], gmVar[[3]], gmVar[[2]], weightEPS, convergeEPS, maxit, showProgress, useRelativeDiff,
  ellipseAxisRatioThreshold)
}







gm2dEMparallel<-function(xCoor, xval, gmWeight, gmMean, gmVar, weightEPS=1e-20, convergeEPS=1e-6, maxit=500L, showProgress=5L, NofCore=7L, useRelativeDiff=1, ellipseAxisRatioThreshold=14.2){
gm2dParallel(xCoor[[1]], xCoor[[2]], xval, gmWeight, gmMean[[1]], gmMean[[2]],
  gmVar[[1]], gmVar[[3]], gmVar[[2]], weightEPS, convergeEPS, maxit, showProgress, NofCore, useRelativeDiff,
  ellipseAxisRatioThreshold)
}






gm2dEMcomponentWise<-function(xCoor, xval, gmWeight, gmMean, gmVar, weightEPS=1e-16, convergeEPS=1e-6, maxit=5000L, showProgress=1L, useRelativeDiff=1, ellipseAxisRatio=14.2){
gm2dComponentWise(xCoor[[1]], xCoor[[2]], xval, gmWeight, gmMean[[1]], gmMean[[2]],
  gmVar[[1]], gmVar[[3]], gmVar[[2]], weightEPS, convergeEPS, maxit, showProgress, useRelativeDiff,
  ellipseAxisRatio)
}









# gm2dComponentWise<-function(NumericVector Long, NumericVector Lat, NumericVector val, NumericVector weight,
#               NumericVector miu1, NumericVector miu2, NumericVector var1, NumericVector var2,
#               NumericVector covar, double weightEPS, double convergeEPS,
#               int maxit, int showProgress, double useRelativeDiff)
















GMKtotalVar<-function(location, GMkernel, locationVar, weightEPS=1e-13, outputCorM=F, useMahalanobisDforWeight=F){

if(outputCorM)
{
  corM=matrix(rep(0,length(locationVar)^2L),ncol=length(locationVar))
  tmp=GMKtotalVarCpp(location, GMkernel, locationVar, corM, weightEPS, outputCorM, useMahalanobisDforWeight)
  return(list(tmp,corM))
}
GMKtotalVarCpp(location, GMkernel, locationVar, as.matrix(0), weightEPS, outputCorM, useMahalanobisDforWeight)
}





findEndColumn<-function(N,cpuN){
N=N-1L
C=N/2.0/cpuN*(N-1)
tmp=cumsum(as.numeric((N-1):1))
rst=sapply(cumsum(rep(C,cpuN)),function(x)which(tmp>x)[1]-1)
#rst[length(rst)]=N
rst=as.integer(rst)-1L
#rst[1]=0L
rst[length(rst)]=N-1L
rst
}





GMKtotalVarParallel<-function(location, GMkernel, locationVar, NofCore,
                              weightEPS=1e-12, outputCorM=F, useMahalanobisDforWeight=F,
                              absorbWeightInLocalCovar=T){

colIndex=rep(0, NofCore)
indexVec=0
corM=as.matrix(0)

# if(equalDivide)colIndex=findEndColumn(length(locationVar),NofCore)
# else
  indexVec=sample(0L:(length(locationVar)-1L),length(locationVar))

if(outputCorM)corM=matrix(rep(0,length(locationVar)^2L),ncol=length(locationVar))

tmp=GMKtotalVarCppParallel(location, GMkernel, locationVar, corM,
                          colIndex, indexVec, weightEPS, outputCorM,
                          absorbWeightInLocalCovar, useMahalanobisDforWeight)

if(outputCorM)return(list(tmp,corM))
return(tmp)

# if(outputCorM)
# {
#   corM=matrix(rep(0,length(locationVar)^2L),ncol=length(locationVar))
#   tmp=GMKtotalVarCppParallel(location, GMkernel, locationVar, corM,
#                              colIndex, indexVec, weightEPS, outputCorM)
#   return(list(tmp,corM))
# }
# GMKtotalVarCppParallel(location, GMkernel, locationVar, corM, colIndex,
#                        weightEPS, outputCorM, equalDivide)
}











