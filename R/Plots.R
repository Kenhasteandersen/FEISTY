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
  
  semilogypanel(xlim=sim$t, ylim=c(1e-2, 100),
                xlab="Time (years)", ylab = "SSB (??)")
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
  rates = calcDerivatives(0, u, p, bFullOutput = TRUE)
  
  if (bPlot)
    defaultplot(mfcol=c(2,1))
  xlim = range(p$mc[p$ixFish])
  #
  # Mortalities:
  #
  loglogpanel(xlim=xlim, ylim=range(rates$mortpred),
              xlab="-", ylab="mort (day^{-1})", xaxis = FALSE)
  for (i in 1:p$nGroups) {
    lines(p$mc[p$ix[[i]]], rates$mortpred[p$ix[[i]]], lwd=i, col='red')
  }
  hline(p$mort0)
  #
  # Feeding level
  # 
  semilogxpanel(xlim = xlim, ylim=c(0,1),
                ylab="Feeding level", xlab="Mass (g)")
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
  
  loglogpanel(xlim=p$mc[p$ixFish], ylim=sim$B[iTime,],
              xlab = "Mass (gram)", ylab="Biomass (gram/??")
  for (i in 1:p$nGroups) {
    lines(p$mc[p$ix[[i]]], sim$B[iTime, p$ix[[i]]-p$ixFish[1]+1], lwd=i)
  }
}
#
# Make 4 panels of simulation
#
plotSimulation = function(sim) {
  defaultplot(mfcol=c(4,1))
  plotSSBtime(sim,bPlot=FALSE)
  plotSpectra(sim, bPlot=FALSE)
  plotRates(sim$p, u=c( sim$R[sim$nTime,], sim$B[sim$nTime,]),bPlot=FALSE)
}
#
# Make a basic run:
#
baserun = function() {
  p = setupBasic()
  sim = simulate(p)
  plotSimulation(sim)
  
  return(sim)
}
