source("FEISTY.R")
source("plottools.R")

#
# Makes a basic plot of the adult biomass (SSB) of all functional groups
# as a function time.
#
plotSSBtime = function(sim, bPlot=TRUE) {
  if (bPlot)
    defaultplot()
  p = sim$p
  
  semilogypanel(xlim=sim$t, ylim=sim$SSB+1e-10,
                xlab="Time (yr)", ylab = "SSB (gww  m-2)")
  #
  # Plot fish
  #
  for (i in 1:p$nGroups)
    lines(sim$t, sim$SSB[,i])
  #
  # Plot resources:
  #
  for (i in p$ixR)
    lines(sim$t, sim$R[,i], col='blue')
}
#
# Plots the mortalities and feeding levels
#
plotRates = function(p, u=p$u0, bPlot=TRUE) {
  if (p$USEdll) {
    rates = calcDerivativesF(0, u, p, bFullOutput = TRUE)
  }else{
    rates = calcDerivativesR(0, u, p, bFullOutput = TRUE)
  }
    
  if (bPlot)
    defaultplot(mfcol=c(3,1))
  xlim = range(p$mc[p$ixFish])
  
  #
  # Growth rate
  # 
  loglogpanel(xlim = xlim, ylim=rates$g+1e-10,
                ylab="Growth rate (1/year)", xlab="-", xaxis = FALSE)
  for (i in 1:p$nGroups) {
    lines(p$mc[p$ix[[i]]], rates$g[p$ix[[i]]-length(p$ixR)], lwd=i, col='black')
  }
  
  #
  # Mortalities:
  #
  loglogpanel(xlim=xlim, ylim=rates$mortpred+1e-10,
              xlab="-", ylab="mort (1/year)", xaxis = FALSE)
  for (i in 1:p$nGroups) {
    lines(p$mc[p$ix[[i]]], rates$mortpred[p$ix[[i]]], lwd=i, col='red')
    lines(p$mc[p$ix[[i]]], p$mortF[p$ix[[i]]], lwd=i, col='blue')
  }
  hline(p$mort0)
  
  legend(x='bottomleft',
         legend=c('Predation','Background'),
         lty=c(1,dotted),
         col=c('red','black'),
         bty='n')
  #
  # Feeding level
  # 
  semilogxpanel(xlim = xlim, ylim=c(0,1),
                ylab="Feeding level, f", xlab="Mass (gww)")
  for (i in 1:p$nGroups) {
    ix = p$ix[[i]]
    lines(p$mc[ix], rates$f[ix], lwd=i)
    lines(p$mc[ix], p$metabolism[ix]/(p$epsAssim*p$Cmax[ix]), lwd=i, lty=dotted) # Critical feeding level
  }
  
  return(rates)
}
#
# Plot the biomasses of all groups
#
plotSpectra = function(sim, iTime=sim$nTime, bPlot=TRUE) {
  if (bPlot)
    defaultplot()
  p = sim$p
  # mean biomass of last 50% of simulation time
  loglogpanel(xlim=p$mc[p$ixFish], ylim=pmax(1e-10,colMeans(sim$B[round(0.5*iTime):iTime,])),
              xlab = "Mass (gww)", ylab="Biomass (gww m-2)")
  for (i in 1:p$nGroups) {
    lines(p$mc[p$ix[[i]]], sim$B[iTime, p$ix[[i]]-p$ixFish[1]+1], lwd=i)
  }
}
#
# Make 4 panels of simulation
#
plotSimulation = function(sim) {
  defaultplot(mfcol=c(5,1))
  plotSSBtime(sim,bPlot=FALSE)
  plotSpectra(sim, bPlot=FALSE)
  rates = plotRates(sim$p, u=c( sim$R[sim$nTime,], sim$B[sim$nTime,]),bPlot=FALSE)
  
  return(rates)
}
#
# Plot the interaction matrix:
#
plotTheta = function(p) {
  par(mfcol=c(1,1))

  image( t(p$theta[p$ixFish,]), 
         x = seq(1,max(p$ixFish)), 
         y = seq(1, length(p$ixFish)),
         xlab="Prey group", ylab="Predator group" )
  #
  # Add a line for the resource groups
  #
  y = c(0,max(p$ixFish))
  lines(x = 0.5+c(1,1)*max(p$ixR), y = y, col="red")
  #
  # Add a line for each fish group:
  #
  for (iGroup in 1:p$nGroups) {
    lines(x = 0.5+c(1,1)*max(p$ix[[iGroup]]), y-max(p$ixR))
    lines(x=y, y=c(1,1)*max(p$ix[[iGroup]])-max(p$ixR)+0.5)
  }
}

#
# Make a basic run:
#
baserun = function(USEdll=TRUE) {
  p = setupBasic()
  sim = simulate(p,tEnd = 100,USEdll)
  plotSimulation(sim)
  return(sim)
}
