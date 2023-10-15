

# makeEmpVStrbPlot = function(condDistrs, trbParams, savePath, )
# {
#   
# }


makeInitialPoints = function(abcLB, abcUB, ndiv = 3L)
{
  abc = apply(rbind(abcLB, abcUB), 2, function(x)
  {
    exp(seq(log(x[1]), log(x[2]), len = ndiv + 2L)[2L:(2L + ndiv - 1L)])
  })
  t(expand.grid(a = abc[,1], b = abc[,2], c = abc[,3]))
}


distTableMean = function(distTable, probsAreColumns = T)
{
  if (!probsAreColumns) distTable = t(distTable)
  apply(distTable, 2, function(x)
  {
    sum(x[-c(1,2)] * seq(0, x[2], len = nrow(distTable) - 2))
  })
}


distTableSD = function(distTable, probsAreColumns = T)
{
  if (!probsAreColumns) distTable = t(distTable)
  apply(distTable, 2, function(x)
  {
    sp = seq(0, x[2], len = nrow(distTable) - 2L)
    p = x[-c(1L,2L)]
    m = sum(sp * p)
    m2 = sum(sp * sp * p)
    sqrt(m2 - m * m)
  })
}


distTableSumProb = function(distTable, probsAreColumns = T)
{
  if (!probsAreColumns) distTable = t(distTable)
  apply(distTable,2,function(x) sum(x[-c(1,2)]))
}


distTableCV = function(distTable, removeNA = T, probsAreColumns = T)
{
  m = distTableMean(distTable, probsAreColumns)
  s = distTableSD(distTable, probsAreColumns)
  rst = s / m
  if (!removeNA) rst
  else rst[is.finite(rst)]
}


checkDistTable = function(rst, threshold = 1e-10, probsAreColumns = T)
{
  if (!probsAreColumns) distTable = t(distTable)
  cvgood = max(diff(distTableCV(rst))) < 1e-10
  meangood = max(abs(distTableMean(rst) - rst[1,])) < 1e-10
  maxgood = min(diff(rst[2,])) > -1e-10
  p0good = max(diff(rst[3,])) < 1e-10
  p1good = min(diff(rst[nrow(rst),])) > -1e-10
  psumgood = max(abs(apply(rst, 2, function(x) sum(x[-c(1,2)])) - 1)) < 1e-10
  pneggood = min(rst) >= 0
  if (all(c(cvgood, meangood, maxgood, p0good, p1good))) return(TRUE)
  if (!pneggood) message("Negative probabilities exist.")
  if (!psumgood) message("Sums of probabilities are not close to 1 enough.")
  if (!cvgood) message("Coefficients of variation are not nonincreasing.")
  if (!meangood) message("Means do not match MDRs.")
  if (!maxgood) message("Maxes are not nondecreasing.")
  if (!p0good) message("P0s are not nonincreasing.")
  if (!p1good) message("P1s are not nondecreasing.")
  FALSE
}



computeWindows = function(xlen, windowSize = 1L, start = 1L, speed = 1L, 
                          returnList = TRUE)
{
  rst = seq(start, xlen + speed + 1L, by = speed)
  rst = rst[rst < xlen - windowSize + 1L]
  rst = c(rst, xlen - windowSize + 1L)
  if (returnList) lapply(rst, function(x) c(x, x + windowSize - 1L))
  else sapply(rst, function(x) c(x, x + windowSize - 1L))
}


levtrbeta = function(abc_lm1, limit = 1, eps = 1e-8, maxit = 100)
{
  d = solve_d(abc_lm1, eps, maxit)
  if (!is.matrix(abc_lm1)) abc_lm1 = as.matrix(abc_lm1)
  actuar::levtrbeta(limit, abc_lm1[1, ] / abc_lm1[3, ], abc_lm1[3, ], 
                    abc_lm1[2, ] / abc_lm1[3, ], scale = d)
}


rtrbeta = function(n, abc_lm1, eps = 1e-8, maxit = 100)
{
  
  d = solve_d(abc_lm1, eps, maxit)
  if (!is.matrix(abc_lm1)) abc_lm1 = as.matrix(abc_lm1)
  actuar::rtrbeta(n, abc_lm1[1, ] / abc_lm1[3, ], abc_lm1[3, ], 
                  abc_lm1[2, ] / abc_lm1[3, ], scale = d)
}


dtrbeta = function(x, abc_lm1, eps = 1e-8, maxit = 100)
{
  d = solve_d(abc_lm1, eps, maxit)
  if (!is.matrix(abc_lm1)) abc_lm1 = as.matrix(abc_lm1)
  actuar::dtrbeta(x, abc_lm1[1,] / abc_lm1[3,], abc_lm1[3,], 
                  abc_lm1[2,] / abc_lm1[3,], scale = d)
}


ptrbeta = function(q, abc_lm1, eps = 1e-8, maxit = 100)
{
  d = solve_d(abc_lm1, eps, maxit)
  if (!is.matrix(abc_lm1)) abc_lm1 = as.matrix(abc_lm1)
  actuar::ptrbeta(q, abc_lm1[1,] / abc_lm1[3,], abc_lm1[3,], 
                  abc_lm1[2,] / abc_lm1[3,], scale = d)
}


qtrbeta = function(p, abc_lm1, eps = 1e-8, maxit = 100)
{
  d = solve_d(abc_lm1, eps, maxit)
  if (!is.matrix(abc_lm1)) abc_lm1 = as.matrix(abc_lm1)
  actuar::qtrbeta(p, abc_lm1[1,] / abc_lm1[3,], abc_lm1[3,], 
                  abc_lm1[2,] / abc_lm1[3,], scale = d)
}




actuar_rtrbeta = function(n, abcd)
{
  if (!is.matrix(abc_lm1)) abc_lm1 = as.matrix(abc_lm1)
  actuar::rtrbeta(n, abcd[1,] / abcd[3,], abcd[3,], 
                  abcd[2,] / abcd[3,], scale = abcd[4,])
}


actuar_levtrbeta = function(abcd)
{
  if (!is.matrix(abcd)) abcd = as.matrix(abcd)
  actuar::levtrbeta(1, abcd[1,] / abcd[3,], abcd[3,], 
                    abcd[2,] / abcd[3,], scale = abcd[4,])
}


actuar_dtrbeta = function(x, abcd, log = F)
{
  if (!is.matrix(abcd)) abcd = as.matrix(abcd)
  actuar::dtrbeta(x, abcd[1,] / abcd[3,], abcd[3,], 
                  abcd[2,] / abcd[3,], scale = abcd[4,], log = log)
}


actuar_ptrbeta = function(q, abcd)
{
  if (!is.matrix(abcd)) abcd = as.matrix(abcd)
  actuar::ptrbeta(q, abcd[1,] / abcd[3,], abcd[3,], 
                  abcd[2,] / abcd[3,], scale = abcd[4,])
}


actuar_qtrbeta = function(p, abcd)
{
  if (!is.matrix(abcd)) abcd = as.matrix(abcd)
  actuar::qtrbeta(p, abcd[1,] / abcd[3,], abcd[3,], 
                  abcd[2,] / abcd[3,], scale = abcd[4,])
}


# sp is some empirical PMF's support.
# Discretize TrB such that the discrete PMF and the empirical are comparable
# on the plot.
actuar_discretize = function(abcd, sp)
{
  d = sp[2] - sp[1]
  spnew = seq(sp[1] - d/2, by = d, len = length(sp))
  pr = actuar_ptrbeta(spnew, abcd)
  rst = data.frame(val = spnew + d/2, P = diff(c(pr, 1)))
  rst$P = rst$P / sum(rst$P)
  rst
}


# Wrap distribution as histogram object.
wrapAsHist = function(X)
{
  val = X$val
  m = (val[-1] + val[-length(val)]) / 2
  sp = c(val[1] - (m[1] - val[1]), m, 
    val[length(val)] + (val[length(val)] - m[length(m)]))
  H = hist(val, breaks = sp, plot = F)
  H$density = X$P; H
}




readDistrTable = function(fpath, majorityMax = 1, 
                          firstColIndex = TRUE,
                          firstDistrIsDegenerate = TRUE,
                          lastDistrIsDegenerate = TRUE)
{
  dat = data.table::setDF(data.table::fread(fpath))
  if (firstColIndex) dat = dat[-1]
  mdr = dat[[1]]
  dat = t(dat[-1])
  dimnames(dat) = NULL
  Npoint = nrow(dat) - 1L
  sp = seq(0, majorityMax, len = Npoint)
  rst = apply(dat, 2, function(x)
  {
    if (x[1] == majorityMax) list(val = sp, P = x[-1])
    else list(val = seq(0, x[1], len = Npoint), P = x[-1])
  })
  if (firstDistrIsDegenerate)
    rst[[1]] = data.frame(val = 0, P = 1)
  if (lastDistrIsDegenerate)
    rst[[length(rst)]] = data.frame(val = 1, P = 1)
  
  
  condDistr = lapply(rst, function(x)
  {
    if (length(x$val) == 1) x
    else
    {
      d = list(val = x$val[-1], P = x$P[-1])
      d$P = d$P / sum(d$P); d
    }
  })
  
  
  condMDR = unlist(lapply(condDistr, function(x) sum(x[[1]] * x[[2]])))
  
  
  list(MDR = mdr, distrs = rst, condMDR = condMDR, condDistr = condDistr)
}




modelP0 = function(mdr, cdr, windows, targetMDRs, 
                   eps = 1e-10, P0modeltype = "hyman", 
                   returnInfo = FALSE, warnUsers = TRUE)
{
  rst = matrix(unlist(lapply(windows, function(x)
  {
    ind = x[1]:x[2]
    r = cdr[ind]
    c(mean(mdr[ind]), sum(r < eps) / length(r))
  })), nrow = 2)
  m = c(0, rst[1,], 1)
  p0 = c(1, rst[2,], 0)
  ind = longestNonDecreasingSubseq(-p0)
  m = m[ind]
  p0 = p0[ind]
  if (P0modeltype == "hyman")
    cbfun = splinefun(x = m, y = p0, method = 'hyman')
  else if (P0modeltype == "linear")
    cbfun = approxfun(x = m, y = p0)
  else stop("P0 model not implemented.")
  P0 = cbfun(targetMDRs)
  if (warnUsers && !all(P0 < 1 - targetMDRs))
    warning(paste0("P0s are too high: P0 >= 1 - targetMDR. ",
            "Only Bernoulli distributions could be fitted to the data. ",
            "Remove some zero claims to make Transformed Beta modeling feasible."))
  if (returnInfo)
  {
    rst = as.data.frame(rst)
    colnames(rst) = c("MDR", "P0")
    list(trainingData = rst, cbfun = cbfun, P0 = P0)  
  }
  else P0
}




# windows is a matrix of 2 rows.
# cdr should all be nonzero.
# P0 responds to targetMDRs.
estimateEmpPMFsOld = function(
    mdr, cdr, windows, P0,
    targetMDRs = c(1:1e4L, seq(1e4L + 10L, 1e5L - 10L, by = 10L)) / 1e5,
    regridMethod = "r4",
    # firstDistr = list(val = 0, P = 1),
    # lastDistr = list(val = 1, P = 1),
    correctBias = FALSE,
    nonDegenerateOldDistrs = NULL,
    sampleWeightOnOldDistrs = 0.5,
    empDistrSupportSize = 20L,
    maxCore = 1000L,
    exponentialScalingForLowerMDRs = TRUE
)
{
  if (length(P0) != length(targetMDRs))
    stop("length(P0) != length(targetMDRs)")
  
  
  if (length(mdr) != length(cdr))
    stop("length(mdr) != length(cdr)")
  
  
  if ( min(diff(mdr)) < 0 ) stop(paste0(
    "mdr has not been sorted. ",
    "Note that cdr should also have been ordered by mdr."))
  
  
  oldDistrs = nonDegenerateOldDistrs
  if (!is.null(oldDistrs))
  {
    if (length(oldDistrs) != length(targetMDRs))
      stop(paste0(
        "List of old distributions is given, ", 
        "but the size of the list is unequal to that of targetMDRs."))
    oldDistrsMeans = unlist(lapply(oldDistrs, function(x) sum(x[[1]] * x[[2]])))
    if (max(abs(oldDistrsMeans - targetMDRs)) > 1e-6)
      stop(paste0("max(abs(oldDistrsMeans - targetMDRs)) > 1e-6. ", 
                  "The old distributions' means should have negligible ",
                  "differences from targetMDRs."))
  }
  
  
  windows = matrix(unlist(windows), nrow = 2)
  
  
  empDistrs = makeEmpDistrList( # Without P0.
    X = cdr, windows = windows, rstSizes = 255L, 
    regridMethod = regridMethod, maxCore = maxCore, 
    fixedMin = 1e300, fixedMax = 1e300)
  # return(empDistrs)
  
  
  m = apply(windows, 2, function(x) mean(mdr[x[1]:x[2]]))
  if (m[1] <= 0) 
    stop("Window has nonpositive mean MDR. Please remove this window of data.")
  
  
  rawP0model = NULL
  targetEmpDistrsP0 = P0
  
  
  mainTargetMDRs = targetMDRs[
    targetMDRs >= m[1] & targetMDRs <= m[length(m)]]
  targetEmpDistrsMain = findEmpDistrGivenMDR( # No degenerate PMFs.
    empDistrs, m, mainTargetMDRs, 255L, maxCore, regridMethod) # mainTargetMDRs is MDRwanted.
  # targetEmpDistrsMain is without P0.
  # return(list(empDistrs = empDistrs, targetEmpDistrsMain = targetEmpDistrsMain,
  #        m = m, mainTargetMDRs = mainTargetMDRs))
  
  
  # ============================================================================
  # Use exponential scaling for extrapolating empirical PMFs whose target 
  #   MDRs are less than the minimum MDR in the given data.
  # ============================================================================
  if (exponentialScalingForLowerMDRs)
  {
    lowerInds = which(targetMDRs < m[1])
    tmpDistr = targetEmpDistrsMain[[1]]
    r = targetMDRs[lowerInds] / targetMDRs[lowerInds[length(lowerInds)] + 1L]
    lowerTargetMeans =  r * sum(tmpDistr[[1]] * tmpDistr[[2]])
    tmpDistr = list(val = tmpDistr[[1]] * r[length(r)], P = tmpDistr[[2]] + 0.0) # Ensure copy.
    # lowerTargetMeans = lowerTargetMeans[-length(lowerTargetMeans)]
    # The seemingly redundant action is to ensure the baseline PMF has
    # less-than-1 max before exponential scaling.
    lowerDistrs = expScalingForMeanMatching(
      distlist = list(tmpDistr),
      lowerTargetMeans[-length(lowerTargetMeans)], 
      regridMethod = regridMethod,
      maxCore = maxCore, eps = 1e-10, 
      maxIter = 100)
    lowerDistrs = c(lowerDistrs, list(tmpDistr))
  } # Without P0
  else
  {
    
    lowerInds = which(targetMDRs < m[1])
    tmpDistr = targetEmpDistrsMain[[1]]
    r = targetMDRs[lowerInds] / targetMDRs[lowerInds[length(lowerInds)] + 1L]
    lowerDistrs = lapply(r, function(mpr)
    {
      list(val = tmpDistr[[1]] * mpr, P = tmpDistr[[2]] + 0.0)
    }) # Without P0
  }
  
  
  # return(list(tmpDistr = tmpDistr, lowerDistrs = lowerDistrs,
  #             lowerTargetMeans = lowerTargetMeans))
  
  
  # ============================================================================
  # Use linear scaling for extrapolating empirical PMFs whose target MDRs are
  #   greater than the maximum MDR in the given data.
  # ============================================================================
  upperInds = which(targetMDRs > m[length(m)])
  tmpDistr = targetEmpDistrsMain[[length(targetEmpDistrsMain)]]
  r = targetMDRs[upperInds] / targetMDRs[upperInds[1] - 1L]
  upperDistrs = upScaleAndBoundPMF(
    pmf = tmpDistr, upperBound = 1.0, upscaler = r, regridMethod = regridMethod)
  
  
  targetEmpDistrs = c(lowerDistrs, targetEmpDistrsMain, upperDistrs) # Without P0.
  
  
  # return(targetEmpDistrs)
  
  
  if (length(targetEmpDistrs) != length(P0))
    stop("length(targetEmpDistrs) != length(P0)")
  
  
  includeP0 = function(distlist)
  {
    for (i in 1:length(distlist))
    {
      distlist[[i]][[1]] = c(0, distlist[[i]][[1]])
      distlist[[i]][[2]] = c(P0[i], distlist[[i]][[2]] * (1 - P0[i]))
    }
    distlist
  }
    
  
  if (is.null(oldDistrs)) # No old distribution for updating.
  {
    # if (empDistrSupportSize == 256L)
    # {
    #   if (correctBias) targetEmpDistrs = correctMeanBias(
    #     distlist = targetEmpDistrs, mdrs = targetMDRs, lim = 1.0, 
    #     eps = 1e-10, maxCore = maxCore)
    #   return(targetEmpDistrs)
    # }
    
    
    targetEmpDistrs = mixDistrList(
      targetEmpDistrs, targetEmpDistrs,
      empDistrSupportSize - 1L, 0, regridMethod, maxCore)
    
    
    # for (i in 1:length(targetEmpDistrsP0))
    #   assignP0(targetEmpDistrs[[i]][[2]], targetEmpDistrsP0[i])
    
    
    if (correctBias) 
    {
      mdr2impose = targetMDRs / (1.0 - P0)
      targetEmpDistrs = correctMeanBias(
        distlist = targetEmpDistrs, mdrs = mdr2impose, lim = 1.0,
        eps = 1e-10, maxCore = maxCore)
    } # Still without P0
    
    
    return(includeP0(targetEmpDistrs))
  }
  
  
  oldDistrsP0 = unlist(lapply(oldDistrs, function(x) x[[2]][1]))
  w = sampleWeightOnOldDistrs
  # targetEmpDistrsP0 = (1 - w) * targetEmpDistrsP0 + w * oldDistrsP0
  P0 = (1 - w) * P0 + w * oldDistrsP0
  oldDistrWithoutP0 = lapply(oldDistrs, function(x)
  {
    list(x[[1]][-1], x[[2]][-1] / (1 - x[[2]][1]))
  })
  
  
  # print("hello")
  # return(oldDistrWithoutP0)
  
  
  targetEmpDistrs = mixDistrList(
    targetEmpDistrs, oldDistrWithoutP0, empDistrSupportSize - 1L,
    sampleWeightOnOldDistrs, regridMethod, maxCore)
  
  
  if (correctBias) 
  {
    mdr2impose = targetMDRs / (1.0 - P0)
    targetEmpDistrs = correctMeanBias(
      distlist = targetEmpDistrs, mdrs = mdr2impose, lim = 1.0,
      eps = 1e-10, maxCore = maxCore)
  }
  
  
  includeP0(targetEmpDistrs)
}




estimateEmpPMFs = function(
    mdr, cdr, windows, P0,
    targetMDRs = c(1:1e4L, seq(1e4L + 10L, 1e5L - 10L, by = 10L)) / 1e5,
    regridMethod = "lmm",
    nonDegenerateOldDistrs = NULL,
    sampleWeightOnOldDistrs = 0.5,
    empDistrSupportSize = 64L,
    maxCore = 1000L
)
{
  if (length(P0) != length(targetMDRs))
    stop("length(P0) != length(targetMDRs)")
  
  
  if (length(mdr) != length(cdr))
    stop("length(mdr) != length(cdr)")
  
  
  if ( min(diff(mdr)) < 0 ) stop(paste0(
    "mdr has not been sorted. ",
    "Note that cdr should also have been ordered by mdr."))
  
  
  oldDistrs = nonDegenerateOldDistrs
  if (!is.null(oldDistrs))
  {
    if (length(oldDistrs) != length(targetMDRs))
      stop(paste0(
        "List of old distributions is given, ", 
        "but the size of the list is unequal to that of targetMDRs."))
    oldDistrsMeans = unlist(lapply(oldDistrs, function(x) sum(x[[1]] * x[[2]])))
    if (max(abs(oldDistrsMeans - targetMDRs)) > 1e-6)
      stop(paste0("max(abs(oldDistrsMeans - targetMDRs)) > 1e-6. ", 
                  "The old distributions' means should have negligible ",
                  "differences from targetMDRs."))
  }
  
  
  windows = matrix(unlist(windows), nrow = 2)
  
  
  empDistrs = makeEmpDistrList( # Without P0.
    X = cdr, windows = windows, rstSizes = 255L, 
    regridMethod = regridMethod, maxCore = maxCore, 
    fixedMin = 1e300, fixedMax = 1e300)
  
  
  m = apply(windows, 2, function(x) mean(mdr[x[1]:x[2]]))
  if (m[1] <= 0) 
    stop("Window has nonpositive mean MDR. Please remove this window of data.")
  
  
  rawP0model = NULL
  targetEmpDistrsP0 = P0
  
  
  mainTargetMDRs = targetMDRs[
    targetMDRs >= m[1] & targetMDRs <= m[length(m)]]
  targetEmpDistrsMain = findEmpDistrGivenMDR( # No degenerate PMFs.
    empDistrs, m, mainTargetMDRs, 255L, maxCore, regridMethod) # mainTargetMDRs is MDRwanted.
  
  
  # ============================================================================
  # Use linear scaling for extrapolating empirical PMFs whose target MDRs are
  #   lower than the minimum MDR in the given data.
  # ============================================================================
  lowerInds = which(targetMDRs < m[1])
  tmpDistr = targetEmpDistrsMain[[1]]
  r = targetMDRs[lowerInds] / targetMDRs[lowerInds[length(lowerInds)] + 1L]
  lowerDistrs = lapply(r, function(mpr)
  {
    list(val = tmpDistr[[1]] * mpr, P = tmpDistr[[2]] + 0.0)
  })
  
  
  # ============================================================================
  # Use linear scaling for extrapolating empirical PMFs whose target MDRs are
  #   greater than the maximum MDR in the given data.
  # ============================================================================
  upperInds = which(targetMDRs > m[length(m)])
  tmpDistr = targetEmpDistrsMain[[length(targetEmpDistrsMain)]]
  r = targetMDRs[upperInds] / targetMDRs[upperInds[1] - 1L]
  upperDistrs = upScaleAndBoundPMF(
    pmf = tmpDistr, upperBound = 1.0, upscaler = r, regridMethod = regridMethod)
  
  
  targetEmpDistrs = c(lowerDistrs, targetEmpDistrsMain, upperDistrs)
  
  
  if (length(targetEmpDistrs) != length(P0))
    stop("length(targetEmpDistrs) != length(P0)")
  
  
  includeP0 = function(distlist)
  {
    for (i in 1:length(distlist))
    {
      distlist[[i]][[1]] = c(0, distlist[[i]][[1]])
      distlist[[i]][[2]] = c(P0[i], distlist[[i]][[2]] * (1 - P0[i]))
    }
    distlist
  }
  
  
  if (is.null(oldDistrs)) # No old distribution for updating.
  {
    targetEmpDistrs = mixDistrList(
      targetEmpDistrs, targetEmpDistrs,
      empDistrSupportSize - 1L, 0, regridMethod, maxCore)
    
    
    mdr2impose = targetMDRs / (1.0 - P0)
    targetEmpDistrsCorrected = correctMeanBias(
      distlist = targetEmpDistrs, mdrs = mdr2impose, lim = 1.0,
      eps = 1e-10, maxCore = maxCore)
    
    
    return(list(biased = includeP0(targetEmpDistrs), 
                biasCorrected = includeP0(
                  targetEmpDistrsCorrected)))
  }
  
  
  oldDistrsP0 = unlist(lapply(oldDistrs, function(x) x[[2]][1]))
  w = sampleWeightOnOldDistrs
  P0 = (1 - w) * P0 + w * oldDistrsP0
  oldDistrWithoutP0 = lapply(oldDistrs, function(x)
  {
    list(x[[1]][-1], x[[2]][-1] / (1 - x[[2]][1]))
  })
  
  
  targetEmpDistrs = mixDistrList(
    targetEmpDistrs, oldDistrWithoutP0, empDistrSupportSize - 1L,
    sampleWeightOnOldDistrs, regridMethod, maxCore)
  
  
  mdr2impose = targetMDRs / (1.0 - P0)
  targetEmpDistrsCorrected = correctMeanBias(
    distlist = targetEmpDistrs, mdrs = mdr2impose, lim = 1.0,
    eps = 1e-10, maxCore = maxCore)
  
  
  list(biased = includeP0(targetEmpDistrs), 
       biasCorrected = includeP0(targetEmpDistrsCorrected))
}








trbLogScore = function(x, p0, abcd, lm = 1, base = exp(1))
{
  nx = length(x)
  p0sum = log(p0, base) * sum(x == 0)
  mul = 1 - p0
  x = x[x != 0]
  lmp = 1 - actuar_ptrbeta(lm, abcd)
  plmSum = sum(x >= lm) * log(lmp * mul, base)
  x = x[x < lm]
  mainSum = actuar_dtrbeta(x, abcd, log = T) / log(base) + 
    length(x) * log(mul, base)
  p0sum + plmSum + mainSum
}


# sp is the support
rawDiscretization = function(sp, abcd, lm = 1, normalize = TRUE)
{
  dx = sp[2] - sp[1]
  lastPoint = sp[length(sp)]
  if (lastPoint < 1) y = c(sp - dx / 2, lastPoint + dx / 2)
  else y = c(sp - dx / 2, Inf)
  p = actuar_ptrbeta(y, abcd)
  rst = list(val = sp, P = diff(p))
  if (normalize) rst$P = rst$P / sum(rst$P)
  rst
}


makePDFplotData = function(min, max, abcd, lm = 1, npoint = 500, normalize = TRUE)
{
  dx = (max - min) / (npoint - 1)
  y = seq(min - dx / 2, max + dx / 2, len = npoint + 1)
  if (max >= lm - 1e-10) y[length(y)] = Inf
  p = actuar_ptrbeta(y, abcd)
  rst = list(val = seq(min, max, len = npoint), P = diff(p))
  if (normalize) rst$P = rst$P / sum(rst$P)
  rst
}


pmftableToDistrList = function(X)
{
  npoint = nrow(X) - 2L
  apply(X, 2, function(x)
  {
    list(val = seq(0, x[2], len = npoint), P = x[-c(1, 2)])
  })
}




# biDirSeqLBFGSBtrbFitList = function(
#     abc,
#     lm1,
#     empDistrList,
#     startingIndex,
#     abcLB,
#     abcUB,
#     scaleEps = 1e-8, 
#     scaleMaxit = 100, 
#     distanceFun = "likelihood",
#     max_iterations = 100, 
#     maxCore = 1000, 
#     RIBlib = "Numerical Recipes", 
#     hgrad = 0, centralDiff = T, m = 6, epsilon = 1e-5, 
#     epsilon_rel = 1e-5, past = 1, delta = 1e-10, max_submin = 10, 
#     max_linesearch = 20, min_step = 1e-20, max_step = 1e+20, 
#     ftol = 1e-4, wolfe = 0.9
#     )
# {
#   startingPmf = startingIndex
#   secondHalf = startingPmf:length(empDistrList)
#   optRstSecondHalf = LBFGSBtrbFitList(
#     abc = abc, 
#     lm1 = lm1[secondHalf],
#     empDistrList = empDistrList[secondHalf], 
#     abcLB = abcLB, 
#     abcUB = abcUB, 
#     scaleEps = scaleEps, 
#     scaleMaxit = scaleMaxit, 
#     distanceFun = distanceFun,
#     max_iterations = max_iterations, 
#     maxCore = maxCore, 
#     RIBlib = RIBlib, 
#     sequentialUpdate = T, 
#     hgrad = hgrad, centralDiff = centralDiff, m = m, epsilon = epsilon, 
#     epsilon_rel = epsilon_rel, past = past, delta = delta, 
#     max_submin = max_submin, max_linesearch = max_linesearch, 
#     min_step = min_step, max_step = max_step, ftol = ftol, wolfe = wolfe)
#   
#   
#   if (startingPmf != 1L)
#   {
#     firstHalf = startingPmf:1  
#     system.time({optRstFirstHalf = LBFGSBtrbFitList(
#       abc = abc, 
#       lm1 = lm1[firstHalf],
#       empDistrList = empDistrList[firstHalf], 
#       abcLB = abcLB, abcUB = abcUB, 
#       scaleEps = 1e-8, scaleMaxit = 100, 
#       distanceFun = distanceFun,
#       max_iterations = 100, maxCore = 1000, 
#       RIBlib = "Numerical Recipes", 
#       sequentialUpdate = sequentialUpdate, 
#       hgrad = 0, centralDiff = T, m = 6, epsilon = 1e-5, 
#       epsilon_rel = 1e-5, past = 1, delta = 1e-10, max_submin = 10, 
#       max_linesearch = 20, min_step = 1e-20, max_step = 1e+20, 
#       ftol = 1e-4, wolfe = 0.9)
#     })
#     optRst = list()
#     optRst$param = cbind(optRstFirstHalf$param[, firstHalf], 
#                          optRstSecondHalf$param[, -1])
#     optRst$fval = c(optRstFirstHalf$fval[firstHalf], 
#                     optRstSecondHalf$fval[-1])
#     optRst$niter = c(optRstFirstHalf$niter[firstHalf], 
#                      optRstSecondHalf$niter[-1])
#   }
#   else optRst = optRstSecondHalf
#   
#   
# }
# 
# 
# ensembleRst[[length(ensembleRst) + 1]] = list(
#   optRst = optRst, empMeans = unlist(
#     lapply(p0andEmpDistrs, function(x) keyALGs::Mean(x))))
# }




# # windows is a matrix of 2 rows.
# modelP0andMakeEmpiricalDistributions = function(
#     mdr, cdr, windows,
#     regridMethod = "r4",
#     givenFixedP0 = NULL,
#     P0modeltype = "hyman", 
#     targetMDRs = c(1:1e4L, seq(1e4L + 10L, 1e5L - 10L, by = 10L)) / 1e5,
#     firstDistr = list(val = 0, P = 1),
#     lastDistr = list(val = 1, P = 1),
#     correctBias = FALSE,
#     nonDegenerateOldDistrs = NULL,
#     sampleWeightOnOldDistrs = 0.5,
#     empDistrSupportSize = 20L,
#     maxCore = 1000L,
#     exponentialScalingForLowerMDRs = TRUE
#     )
# {
#   
#   if (length(mdr) != length(cdr))
#     stop("length(mdr) != length(cdr)")
#   
#   
#   if ( min(diff(mdr)) < 0 ) stop(paste0(
#       "mdr has not been sorted. ",
#       "Note that cdr should also have been ordered by mdr."))
#   
#   
#   if (!is.null(givenFixedP0) && length(givenFixedP0) != length(targetMDRs))
#     stop("P0 is given but its size is unequal to targetMDRs.")
#   
#   
#   oldDistrs = nonDegenerateOldDistrs
#   if (!is.null(oldDistrs))
#   {
#     if (length(oldDistrs) != length(targetMDRs))
#       stop(paste0(
#         "List of old distributions is given, ", 
#         "but the size of the list is unequal to that of targetMDRs."))
#     oldDistrsMeans = unlist(lapply(oldDistrs, function(x) sum(x[[1]] * x[[2]])))
#     if (max(abs(oldDistrsMeans - targetMDRs)) > 1e-6)
#       stop(paste0("max(abs(oldDistrsMeans - targetMDRs)) > 1e-6. ", 
#                   "The old distributions' means should have negligible ",
#                   "differences from targetMDRs."))
#   }
#     
#   
#   windows = matrix(unlist(windows), nrow = 2)
#   
#   
#   empDistrs = makeEmpDistrList(
#     X = cdr, windows = windows, rstSizes = 256L, 
#     regridMethod = regridMethod, maxCore = maxCore, 
#     fixedMin = 0, fixedMax = 1e300)
#   
#   
#   m = apply(windows, 2, function(x) mean(mdr[x[1]:x[2]]))
#   if (m[1] <= 0) stop("Window has nonpositive mean MDR. Please remove this window of data.")
#   
#   
#   if (is.null(givenFixedP0))
#   {
#     p0 = unlist(lapply(empDistrs, function(x) x[[2]][1]))
#     m_  = c(0, m, 1)
#     p0_ = c(1, p0, 0)
#     ind = longestIncreasingSubseq(-p0_)
#     m_ = m_[ind]
#     p0_ = p0_[ind]
#     if (P0modeltype == "hyman")
#       cbfun = splinefun(x = m_, y = p0_, method = 'hyman')
#     else if (P0modeltype == "linear")
#       cbfun = approxfun(x = m_, y = p0_)
#     else stop("P0 model not implemented.")
#     rawP0model = list(mdr = m_, p0 = p0_, p0fun = cbfun)
#     targetEmpDistrsP0 = cbfun(targetMDRs)
#   }
#   else
#   {
#     rawP0model = NULL
#     targetEmpDistrsP0 = givenFixedP0
#   }
#   
#   
#   mainTargetMDRs = targetMDRs[
#     targetMDRs >= m[1] & targetMDRs <= m[length(m)]]
#   targetEmpDistrsMain = findEmpDistrGivenMDR( # No degenerate PMFs.
#     empDistrs, m, mainTargetMDRs, 256L, maxCore, regridMethod)
#   
#   
#   # ============================================================================
#   # Use exponential scaling for extrapolating empirical PMFs whose target 
#   #   MDRs are less than the minimum MDR in the given data.
#   # ============================================================================
#   if (exponentialScalingForLowerMDRs)
#   {
#     lowerInds = which(targetMDRs < m[1])
#     tmpDistr = targetEmpDistrsMain[[1]]
#     r = targetMDRs[lowerInds] / targetMDRs[lowerInds[length(lowerInds)] + 1L]
#     lowerTargetMeans =  r * sum(tmpDistr[[1]] * tmpDistr[[2]])
#     tmpDistr = list(val = tmpDistr[[1]] * r[length(r)], P = tmpDistr[[2]] + 0.0) # Ensure copy.
#     # lowerTargetMeans = lowerTargetMeans[-length(lowerTargetMeans)]
#     # The seemingly redundant action is to ensure the baseline PMF has
#     # less-than-1 max before exponential scaling.
#     lowerDistrs = expScalingForMeanMatching(
#       distlist = list(tmpDistr),
#       lowerTargetMeans[-length(lowerTargetMeans)], 
#       regridMethod = regridMethod,
#       maxCore = maxCore, eps = 1e-10, 
#       maxIter = 100)
#     lowerDistrs = c(lowerDistrs, list(tmpDistr))
#   }
#   else
#   {
#     
#     lowerInds = which(targetMDRs < m[1])
#     tmpDistr = targetEmpDistrsMain[[1]]
#     r = targetMDRs[lowerInds] / targetMDRs[lowerInds[length(lowerInds)] + 1L]
#     lowerDistrs = lapply(r, function(mpr)
#     {
#       list(val = tmpDistr[[1]] * mpr, P = tmpDistr[[2]] + 0.0)
#     })
#   }
#   
#   
#   # return(list(tmpDistr = tmpDistr, lowerDistrs = lowerDistrs,
#   #             lowerTargetMeans = lowerTargetMeans))
#   
#   
#   # ============================================================================
#   # Use linear scaling for extrapolating empirical PMFs whose target MDRs are
#   #   greater than the maximum MDR in the given data.
#   # ============================================================================
#   upperInds = which(targetMDRs > m[length(m)])
#   tmpDistr = targetEmpDistrsMain[[length(targetEmpDistrsMain)]]
#   r = targetMDRs[upperInds] / targetMDRs[upperInds[1] - 1L]
#   upperDistrs = upScaleAndBoundPMF(
#     pmf = tmpDistr, upperBound = 1.0, upscaler = r, regridMethod = regridMethod)
#   
#   
#   # upperDistrs2 = list()
#   # for (i in 1:length(upperInds))
#   # {
#   #   tmp = list(val = r[i] * tmpDistr[[1]], P = tmpDistr[[2]] + 0)
#   #   ind = which(tmp$val > 1)
#   #   if (length(ind) != 0)
#   #   {
#   #     psum = sum(tmp$P[ind])
#   #     tmp$val = tmp$val[-ind]
#   #     tmp$P = tmp$P[-ind]
#   #     sz = length(tmp$val)
#   #     if (tmp$val[sz] == 1) tmp$P[sz] = tmp$P[sz] + psum
#   #     else { tmp$val = c(tmp$val, 1); tmp$P = c(tmp$P, psum) }
#   #   }
#   #   upperDistrs2[[i]] = tmp
#   # }
#   # return (  list(upperDistrs, upperDistrs2) )
#   
#   
#   targetEmpDistrs = c(lowerDistrs, targetEmpDistrsMain, upperDistrs)
#   # return(targetEmpDistrs)
#   
#   
#   for (i in 1:length(targetEmpDistrsP0))
#     assignP0(targetEmpDistrs[[i]][[2]], targetEmpDistrsP0[i])
#   
#   
#   if (is.null(oldDistrs)) # No old distribution for updating.
#   {
#     if (empDistrSupportSize == 256L)
#     {
#       if (correctBias) targetEmpDistrs = correctMeanBias(
#         distlist = targetEmpDistrs, mdrs = targetMDRs, lim = 1.0, 
#         eps = 1e-10, maxCore = maxCore)
#       return(list(empDistrs = targetEmpDistrs, rawP0model = rawP0model,
#                   targetMDRs = targetMDRs, P0 = targetEmpDistrsP0))
#     }
#     
#     
#     targetEmpDistrs = mixDistrList(
#       targetEmpDistrs, targetEmpDistrs,
#       empDistrSupportSize, 0, regridMethod, maxCore)
#     
#     
#     for (i in 1:length(targetEmpDistrsP0))
#       assignP0(targetEmpDistrs[[i]][[2]], targetEmpDistrsP0[i])
#     
#     
#     if (correctBias) targetEmpDistrs = correctMeanBias(
#       distlist = targetEmpDistrs, mdrs = targetMDRs, lim = 1.0, 
#       eps = 1e-10, maxCore = maxCore)
#     
#     
#     return(list(empDistrs = targetEmpDistrs, rawP0model = rawP0model,
#                 targetMDRs = targetMDRs, P0 = targetEmpDistrsP0))
#   }
#   
#   
#   oldDistrsP0 = unlist(lapply(oldDistrs, function(x) x[[2]][1]))
#   w = sampleWeightOnOldDistrs
#   targetEmpDistrsP0 = (1 - w) * targetEmpDistrsP0 + w * oldDistrsP0
#   
#   
#   targetEmpDistrs = mixDistrList(
#     targetEmpDistrs, oldDistrs, empDistrSupportSize, 
#     sampleWeightOnOldDistrs, regridMethod, maxCore)
#   
#   
#   for (i in 1:length(targetEmpDistrsP0))
#     assignP0(targetEmpDistrs[[i]][[2]], targetEmpDistrsP0[i])
#   
#   
#   if (correctBias) targetEmpDistrs = correctMeanBias(
#     distlist = targetEmpDistrs, mdrs = targetMDRs, lim = 1.0, 
#     eps = 1e-10, maxCore = maxCore)
#   
#   
#   list(empDistrs = targetEmpDistrs, 
#        rawP0model = rawP0model, targetMDRs = targetMDRs, 
#        P0 = targetEmpDistrsP0)
# }




















