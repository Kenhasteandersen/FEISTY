#
#Nov 2023
#

#
# Analysis of the fishing mortality
#

analyseStages = function(nStages = c(3,6,9,18,36), maxF = 10) {
  y = list()
  for (i in 1:length(nStages)) {
    y[[i]] = plotYield(p=setupBasic2(szprod = 100,lzprod = 100, bprod  = 5,
                                     depth  = 100,
                                     Tp = 10,
                                     Tb = 8, 
                                     nStages=nStages,
                                     etaMature=0.25,
                                     F=0, # overwritten later
                                     etaF=0.05))
  }
  F = y[[1]][[1]]
  
  defaultplot()
  defaultpanel(xlim=F, ylim=c(0,max(y[[length(nStages)]][[2]])),
               xlab="Fishing mortality (year$^{-1}$)",
               ylab="Yield (g/$m^2$/yr)")
  for (i in 1:length(nStages)) {
    lines(F, y[[i]][[2]], lwd=i)
  }
}


#
# Make a yield curve
#
plotYield = function(p=setupBasic2(szprod = 100,lzprod = 100, bprod  = 5,
                                               depth  = 100,
                                               Tp = 10,
                                               Tb = 8, 
                                               nStages=6,
                                               etaMature=0.25,
                                               F=0, # overwritten later
                                               etaF=0.05), 
                     maxF=3, nF=20, F=seq(0,maxF,length.out=nF)) {
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
    p = setFishing(p,F[i])
    sim = simulateFeisty(bUseRDerivative    = FALSE,
                         p, 
                         #tEnd   = 100,
                         #tStep  = 1,
                         #times  = seq(from=0, to=tEnd, by=tStep),  
                         yini   = p$u0,  
                         USEdll = TRUE,
                         Rmodel = derivativesFeistyR,
                         simpleOutput = TRUE) # Simulate
    yy = calcYield(sim)
    yield[i] = sum(sim$yieldMean)#sum(yy[[1]])
    yieldMin[i] = sum(sim$yieldMin)#sum(yy[[2]])
    yieldMax[i] = sum(sim$yieldMax)#sum(yy[[3]])
    
    SSB[i] = sum(sim$SSBMean)#sum(calcSSB(sim)[[1]])
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

#
# F: fishing mortality 1/yr   if F=0 return the original param set
# etaF: the coefficient determining the fish size with 50% fishing selectivity
setFishing = function(p, F=0, etaF=0.05) {
  p$F=F
  p$etaF=etaF
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
  sim$yield=yield
  return(sim)
}

#
# Return the SSB of all function groups
#
calcSSB = function(
    sim,          # The simulation object to analyse
    etaTime=0.5) {# The fraction of the time series to integrate (default the last half) 
  
  p=sim$p
  
  SSBAllgrid = matrix(nrow=sim$nTime, ncol=length(p$ixFish))
  SSB = matrix(nrow=sim$nTime, ncol=p$nGroups)
  SSBMean = rep(data=0, p$nGroups)
  SSBMin = SSBMean
  SSBMax = SSBMean
  
  ixTime = which(sim$t>=(etaTime*sim$t[sim$nTime]))
      
  for (iGroup in 1:p$nGroups) {
    ix = p$ix[[iGroup]]
    SSBAllgrid[,ix-length(p$ixR)] =t(t(sim$u[, ix]) * p$psiMature[ix]) 
    SSB[,iGroup]= rowSums(SSBAllgrid[,ix-length(p$ixR)])
    #deltaM = p$mUpper[ix]-p$mLower[ix]

    SSBMean[iGroup] = exp(mean(log(rowSums(SSBAllgrid[ixTime,ix-max(p$ixR)]+1e-10))))
    SSBMin[iGroup] = min(rowSums( SSBAllgrid[ixTime,ix-max(p$ixR)] ))
    SSBMax[iGroup] = max(rowSums( SSBAllgrid[ixTime,ix-max(p$ixR)] ))
  }
  
  sim$SSBMean=SSBMean
  sim$SSBMin=SSBMin
  sim$SSBMax=SSBMax
  sim$SSB=SSB

  return(sim)
}
