source("R/FEISTY.R")
source("R/plottools.R")
#
# Makes a basic plot of the adult biomass (SSB) of all functional groups
# as a function time.
#
plotSSBtime = function(sim) {
  defaultplot()

  semilogypanel(xlim=sim$t, ylim=c(1e-5, 100),
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
# Make a basic run:
#
baserun = function() {
  p = setupBasic()
  sim = simulate(p)
  plotSSBtime(sim)
  
  return(sim)
}
