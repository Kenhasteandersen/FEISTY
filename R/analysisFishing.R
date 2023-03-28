#
# Analysis of the fishing mortality
#

analyseStages = function(nStages = c(3,6,9,18,36), mInf=20, maxF = 10) {
  y = list()
  for (i in 1:length(nStages)) {
    y[[i]] = plotYield(nStages[i],mInf=mInf,maxF=maxF)
  }
  F = y[[1]][[1]]
  
  defaultplot()
  defaultpanel(xlim=F, ylim=c(0,max(y[[length(nStages)]][[2]])),
               xlab="Fishing mortality (year$^{-1}$)",
               ylab="Yield (gC/vol/yr)")
  for (i in 1:length(nStages)) {
    lines(F, y[[i]][[2]], lwd=i)
  }
}


#
# Make a yield curve
#
plotYield = function(nStages=20, mInf=20, p=setupPelagicSpecies(nStages = nStages, mInf=mInf), 
                     maxF=3, nF=10, F=seq(0,maxF,length.out=nF)) {
  #
  # Fishing with a trawl selectivity
  #
  
  #
  # Run a yield curve:
  #
  yield = rep(0,length(F))
  yieldMin = rep(0,length(F))
  yieldMax = rep(0,length(F))
  SSB = rep(0,length(F))
  
  for (i in 1:length(F)) {
    p = setFishing(p,rep(F[i], p$nGroups))
    sim = simulate(p, tEnd=200, USEdll = FALSE) # Simulate
    yy = calcYield(sim)
    yield[i] = sum(yy[[1]])
    yieldMin[i] = sum(yy[[2]])
    yieldMax[i] = sum(yy[[3]])
    
    SSB[i] = sum(calcSSB(sim)[[1]])
  }
  #
  # Plot
  #
  defaultplot()
  defaultpanel(xlim=F, ylim=c(yieldMax,SSB),
               xlab="Fishing mortality (year$^{-1}$)",
               ylab="Yield (g/$m^2$/yr)")
  lines(F,yield,lwd=3)
  lines(F,yieldMin)
  lines(F,yieldMax)
  lines(F, SSB,col="red", lwd=3)
  
  return(list(F,yield, yieldMin, yieldMax))
}

setFishing = function(p, F, etaF=0.05 ) {
  
  for (iGroup in 1:p$nGroups) {
    ix = p$ix[[iGroup]]
    mFishing = etaF*max(p$mUpper[ix]) # selectivity at 0.05 of maximum size
    psi = ( 1 + (p$mc[ix]/mFishing)^(-3) )^(-1) # Standard trawl selectivity from Andersen (2019) Fig 5.2
    p$mortF[ix] = psi*F[iGroup] 
  }
  return(p)
}
  
#
# Return the yield of all function groups
#
calcYield = function(
    sim,          # The simulation object to analyse
    etaTime=0.5) {# The fraction of the time series to integrate (default the last half) 
  
  yieldMean = rep(data=0, sim$p$nGroups)
  yieldMin = yieldMean
  yieldMax = yieldMean

  for (iGroup in 1:sim$p$nGroups) {
    ix = sim$p$ix[[iGroup]]
    ixTime = which(sim$t>=(etaTime*sim$t[sim$nTime]))
    #deltaM = sim$p$mUpper[ix]-sim$p$mLower[ix]
    B = matrix(0,length(ixTime),length(ix))
    for (i in 1:length(ixTime)) {
      B[i,] = sim$p$mortF[ix] * sim$B[ixTime[i],ix-max(sim$p$ixR)] 
    }
    yieldMean[iGroup] = exp(mean(log(rowSums(B+1e-10 ))))
    yieldMin[iGroup] = min(rowSums( B ))
    yieldMax[iGroup] = max(rowSums( B ))
  }
  
  return(list(yieldMean,yieldMin,yieldMax))
}
#
# Return the SSB of all function groups
#
calcSSB = function(
    sim,          # The simulation object to analyse
    etaTime=0.5) {# The fraction of the time series to integrate (default the last half) 
  
  SSBMean = rep(data=0, sim$p$nGroups)
  SSBMin = SSBMean
  SSBMax = SSBMean
  
  for (iGroup in 1:sim$p$nGroups) {
    ix = sim$p$ix[[iGroup]]
    ixTime = which(sim$t>=(etaTime*sim$t[sim$nTime]))
    #deltaM = sim$p$mUpper[ix]-sim$p$mLower[ix]
    B = matrix(0,length(ixTime),length(ix))
    for (i in 1:length(ixTime)) {
      B[i,] = sim$p$psiMature[ix] * sim$B[ixTime[i],ix-max(sim$p$ixR)] 
    }
    SSBMean[iGroup] = exp(mean(log(rowSums(B+1e-10 ))))
    SSBMin[iGroup] = min(rowSums( B ))
    SSBMax[iGroup] = max(rowSums( B ))
  }
  
  return(list(SSBMean,SSBMin,SSBMax))
}