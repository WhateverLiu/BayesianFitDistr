


#'
#' Compute empirical distribution.
#' 
#' Compute the empirical PMF of a sample set.
#' 
#' @param x A numeric vector of samples.
#' 
#' @param N An integer. Support size of the result PMF. Default 64.
#'
#' @param rgMethod Regriding method:
#' 
#' \code{"lr"}: linear regrid. Default.
#' 
#' \code{"r4"}: four-point regrid.
#' 
#' \code{"r3"}: three-point regrid. For author's personal use.
#' 
#' \code{"r4g"}: four-point regrid with Gaussian kernel. For author's personal use.
#' 
#' 
#' @inherit rglr return
#' 
#' @details The function creates a PMF with \code{x} as the support and
#' \code{1/length(x)} as the probabilities, and then regrids the PMF onto a 
#' regular grid of \code{N} points within \code{[min(x), max(x)]}.
#'
#' @example inst/examples/emp.R
#'
#'
emp <- function(x, N = 64L, rgMethod = "lr"){
  x = data.table::data.table(
    val = sort(x), P = 1 / length(x))[, .(P = sum(P)), by = val]
  data.table::setDF(x)
  ngd = seq(x$val[1], x$val[nrow(x)], len = N)
  ngd[length(ngd)] = x$val[nrow(x)]
  if      (rgMethod == "lr") rglr(x, ngd)
  else if (rgMethod == "r3") rglr3split(x, ngd)
  else if (rgMethod == "r4") rglr4split(x, ngd)
  else rglr4splitGuassianKernal(x, ngd)
}


# Mean <- function(X){ sum(X[[1]] * X[[2]]) }
# Var <- function(X)
# {
#   sum(X[[1]] ^ 2 * X[[2]]) - sum(X[[1]] * X[[2]]) ^ 2
# }


# como(List X, List Y, NumericVector ngd = NumericVector(0), String rgMethod = "lr") 
# counterComo = function(X, Y, N = 64L, rgMethod = "lr")


mixConvSpatialClaimDistList <- function(
  distList, gridIndex, cors, regridMethod = 'r', headTrunc = 1e-100, 
  tailTrunc = 1e-100, N = 256L, forDeltaCompare = 256L, single = 1L, useFFT = 0)
{
  rst = mixConvSample(list(1L:length(distList)), distList, gridIndex, cors, 
                      regridMethod, headTrunc, tailTrunc, N, forDeltaCompare, 
                      single, maxCore = 1L, useFFT)
  rst[[2]] = rst[[2]][[1]]
  rst
}


mixConvSpatialClaimSampleList = mixConvSample


geocodeToGridIndex <- function(geocode, hierachySize = c(1 / 120, 5 / 120), 
                             startLat = -90, startLong = -180)
{
  rst = lapply(as.list(hierachySize),function(x)
  {
    data.frame(as.integer((geocode[[1]] - startLong) / x),
               as.integer((geocode[[2]] - startLat) / x))
  })
  rst = as.data.frame(rst)
  colnames(rst)[seq(1L,  ncol(rst) - 1L, 2L)] = 
    paste0("long", 1:length(hierachySize))
  colnames(rst)[-seq(1L, ncol(rst) - 1L, 2L)] = 
    paste0("lat",  1:length(hierachySize))
  rst
}


covarMatSumOneSample <- function(variances, gridInd, cors)
{
  covarMatSum(list(1:length(variances)), variances, cors, gridInd, maxCore = 1L)
}


convBtRgd <- function(X, Y, rgMethod = "lr", N = 64L)
{
  rst = convBt(X, Y)
  sp = seq(rst$val[1], rst$val[nrow(rst)], len = N)
  sp[length(sp)] = rst$val[nrow(rst)]
  if(rgMethod == "lr") rglr(rst, sp)
  else rglr4split(rst, sp)
}


rglr2D <- function(dist, xnew, ynew)
{
  indx = as.integer((dist[[1]] - xnew[1]) / (xnew[2] - xnew[1]))
  indy = as.integer((dist[[2]] - ynew[1]) / (ynew[2] - ynew[1]))
  tmpOrder = order(indx, indy)
  rglr2Dcpp(dist[tmpOrder, ], xnew, ynew)
}


normalizeDistributionWithAtoms <- function(X, support = NULL, min = 0, max, mean)
{
  if (min == 0) min = rep(0,length(max))
  if (length(X) == 1L) X = list(X)
  if (is.null(support)){
    rst = mapply(function(p, min, max, m)
    {
      support = seq(min,max,len=length(p))
      midSupport = support[-c(1,length(support))]
      midP  = p[-c(1,length(p))]
      pHead = p[1]
      pTail = p[length(p)]
      head  = support[1]
      tail  = support[length(support)]
      midMean = sum(midSupport * midP)
      midPsum = sum(midP)
      if (midPsum < 1e-16) return(p)
      
      
      # k=solve(matrix(c(pHead+pTail,pTail*tail+head*pHead,midPsum,midMean),ncol=2),c(1,m))
      A = pHead + pTail
      B = pTail * tail + head * pHead
      C = midPsum
      D = midMean
      U=1
      V=m
      k=c((C*V-D*U)/(B*C-A*D),(A*V-B*U)/(A*D-C*B))
      if(min(k)>=0)return(c(k[1]*pHead,k[2]*midP,k[1]*pTail))
      
      
      # k=solve(matrix(c(pTail,pTail*tail,midPsum,midMean),ncol=2),c(1-pHead,m-head*pHead))
      A=pTail
      B=pTail*tail
      C=midPsum
      D=midMean
      U=1-pHead
      V=m-head*pHead
      k=c((C*V-D*U)/(B*C-A*D),(A*V-B*U)/(A*D-C*B))
      if(min(k)>=0)return(c(pHead,k[2]*midP,k[1]*pTail))
      
      
      # k=solve(matrix(c(pHead,head*pHead,midPsum,midMean),ncol=2),c(1-pTail,m-pTail*tail))
      A=pHead
      B=head*pHead
      C=midPsum
      D=midMean
      U=1-pTail
      V=m-pTail*tail
      k=c((C*V-D*U)/(B*C-A*D),(A*V-B*U)/(A*D-C*B))
      c(k[1]*pHead,k[2]*midP,pTail)
      
    },X,min,max,mean,SIMPLIFY = F)
    #print(1.1)
    if(is.data.frame(X))return(as.data.frame(rst))
    return(rst)
  }
  
  rst=mapply(function(p,support,min,max,m){
    midSupport=support[-c(1,length(support))]
    midP=p[-c(1,length(p))]
    pHead=p[1]
    pTail=p[length(p)]
    head=support[1]
    tail=support[length(support)]
    midMean=sum(midSupport*midP)
    midPsum=sum(midP)
    if(midPsum<1e-16)return(p)
    
    # k=solve(matrix(c(pHead+pTail,pTail*tail+head*pHead,midPsum,midMean),ncol=2),c(1,m))
    A=pHead+pTail
    B=pTail*tail+head*pHead
    C=midPsum
    D=midMean
    U=1
    V=m
    k=c((C*V-D*U)/(B*C-A*D),(A*V-B*U)/(A*D-C*B))
    if(min(k)>=0)return(c(k[1]*pHead,k[2]*midP,k[1]*pTail))
    
    
    # k=solve(matrix(c(pTail,pTail*tail,midPsum,midMean),ncol=2),c(1-pHead,m-head*pHead))
    A = pTail
    B = pTail * tail
    C = midPsum
    D = midMean
    U = 1 - pHead
    V = m - head * pHead
    k = c((C * V-D * U) / (B * C - A * D), (A * V - B * U) / (A * D - C * B))
    if (min(k) >= 0) return ( c(pHead, k[2] * midP, k[1] * pTail) )
    
    
    # k=solve(matrix(c(pHead,head*pHead,midPsum,midMean),ncol=2),c(1-pTail,m-pTail*tail))
    A = pHead
    B = head * pHead
    C = midPsum
    D = midMean
    U = 1 - pTail
    V = m - pTail * tail
    k = c((C * V - D * U) / (B * C - A * D),(A * V - B * U) / (A * D - C * B))
    c(k[1] * pHead, k[2] * midP, pTail)
    
  }, X, support, min, max, mean, SIMPLIFY = F)
  if (is.data.frame(X)) return(as.data.frame(rst))
  rst
}




splitList <- function(x, N = 7)
{
  tmp = as.integer(seq(1L, length(x) + 1L, len = N + 1))
  tmp = unique(tmp)
  rst = list()
  for(i in 1L : (length(tmp) - 1L))
  {
    rst[[i]] = x[tmp[i] : (tmp[i + 1] - 1L)]
  }
  rst
}


snowCplyWindows <- function(X, FUN, ..., maxCore = 7)
{
  X = splitList(X, maxCore)
  maxCore = min(maxCore, length(X))
  cl = snow::makeCluster(maxCore, type = "SOCK")
  rst = unlist(snow::clusterApply(cl, X, function(x, FUN)
  {
    lapply(x, function(u)
    {
      FUN(u, ...)
    })
  }, FUN, ...), recursive = F)
  snow::stopCluster(cl)
  rst
}









