

# distrTable is a data frame where each column length is either 42 + 2 or 64 + 2,
# the NGFM style.
# repVal is the replacement value vector and cannot be zero.
makeLossDistrFromDtable = function(distrTable, givenMDRs, repVals)
{


  maxes = as.numeric(distrTable[2, ])
  Npoint = nrow(distrTable) - 2L
  distrTable = as.list(distrTable[-c(1, 2), , drop = F])
  sp = seq(0, 1, len = Npoint)
  distrTable[[1]] = data.frame(val = maxes[1], P = 1)
  Ndistr = length(maxes)
  if(length(unique(distrTable[[Ndistr]])) == 2L) # Probabilities of the last distribution.
    distrTable[[Ndistr]] = data.frame(val = maxes[Ndistr], P = 1)
  # print(length(distrTable))
  # print(str(distrTable))
  # print(length(maxes))
  # return()


  distrTable[-c(1, Ndistr)] = mapply(function(x, y)
  {
    data.frame(val = sp * y, P = x)
  }, distrTable[-c(1, Ndistr)], maxes[-c(1, Ndistr)], SIMPLIFY = F)
  # return(distrTable)


  mdrs = unlist(lapply(distrTable, function(x) sum(x[[1]] * x[[2]])))
  if(max(givenMDRs) > mdrs[length(mdrs)]) stop("givenMDRs exceeds maximal mean from distrTable")
  I = lowerBound(givenMDRs, mdrs)
  distrTable = distrTable[I]
  mdrs = mdrs[I]
  multiplier = givenMDRs / mdrs
  multiplier[is.nan(multiplier)] = 1
  multiplier = multiplier * repVals
  mapply(function(x, y)
  {
    x$val = x$val * y; x
  }, distrTable, multiplier, SIMPLIFY = F)


}






















