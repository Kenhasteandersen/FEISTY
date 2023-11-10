#
#Nov 2023
#


#
# F: fishing mortality 1/yr   if F=0 return the original param set
# etaF: the coefficient determining the fish size with 50% fishing selectivity
setFishing = function(p, F=0, etaF=0.05) {
  if(F==0) return(p)
  for (iGroup in 1:p$nGroups) {
    ix = p$ix[[iGroup]]
    mFishing = etaF*max(p$mUpper[ix]) # selectivity at 0.05 of maximum size
    psi = ( 1 + (p$mc[ix]/mFishing)^(-3) )^(-1) # Standard trawl selectivity from Andersen (2019) Fig 5.2
    p$mortF[ix] = psi*F
  }
  return(p)
}


#
# Return the yield of all function groups
#
calcYield = function(
    sim,          # The simulation object to analyse
    etaTime=0.5) {# The fraction of the time series to integrate (default the last half) 
  
  p=sim$p
  
  yieldAllgrid = matrix(nrow=sim$nTime, ncol=length(p$ixFish))
  yield = matrix(nrow=sim$nTime, ncol=p$nGroups)
  yieldMean = rep(data=0, p$nGroups)
  yieldMin = yieldMean
  yieldMax = yieldMean
  
  ixTime = which(sim$t>=(etaTime*sim$t[sim$nTime]))
  
  for (iGroup in 1:p$nGroups) {
    ix = p$ix[[iGroup]]    
    yieldAllgrid[,ix-length(p$ixR)] =t(t(sim$u[, ix]) * p$mortF[ix]) 
    yield[,iGroup]= rowSums(yieldAllgrid[,ix-length(p$ixR)])

    #deltaM = p$mUpper[ix]-p$mLower[ix]
    yieldMean[iGroup] = exp(mean(log(rowSums(yieldAllgrid[ixTime,ix-max(p$ixR)]+1e-10))))
    yieldMin[iGroup] = min(rowSums( yieldAllgrid[ixTime,ix-max(p$ixR)] ))
    yieldMax[iGroup] = max(rowSums( yieldAllgrid[ixTime,ix-max(p$ixR)] ))
  }

  sim$yieldMean=yieldMean
  sim$yieldMin=yieldMin
  sim$yieldMax=yieldMax
  sim$yield2=yield
  return(sim)
}

