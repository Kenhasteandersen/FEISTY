#
# Code for calculating carbon fluxes, carbon injection, and carbon sequestration.
#
# Written by Julie Lemoine and Ken H Andersen.
# Sequestration calculations based on code by Andre W Visser
#
library(FEISTY)
library(data.table)
library(pracma)

#
# Calculate the flux from carcasses, fecal pellets,  reproduction wastes, and respiration.
# The fluxes are calculated at the position of each size class in units of
# gWW/m2/year.
#
# In:
#  A simulation object
#
# Out:
#  A simulation object with additional fields:
#   fluxCarcass, fluxFecal,fluxRepro, and fluxRespiration
#  all in units: gWW/m2/year
#
calcCarbonFluxes <- function(sim) {
  
  p = sim$p # parameters
  
  # Fluxes from grazing
  grazrate = p$Cmax[p$ixFish] * sim$f # grazing rate [yr^-1]
  graz = sim$B * grazrate # grazing flux (before assimilation) [gWW.m^-2.yr^-1]
  feces = graz * 0.15 # feces flux [gWW.m^-2.yr^-1] -> (1 - p$epsAssim)
  
  # Reproduction waste goes to fecal pellets
  rep2feces = sim$Repro
  rep2feces[, sapply(p$ix, tail, n = 1) - p$nResources] =
    sim$Repro[, sapply(p$ix, tail, n = 1) - p$nResources] +
    sim$Fout[, sapply(p$ix, tail, n = 1) - p$nResources]
  rep2feces = rep2feces * (1 - unique(p$epsRepro)) * unique(p$epsRepro / 0.22) # eps_egg = 0.22 from Andersen 2019 p47
  
  feces = feces + rep2feces # feces flux [gWW.m^-2.yr^-1]
  
  # Flux from respiration
  resprate = p$metabolism[p$ixFish] + 0.15 * grazrate # respiration rate [yr-1]
  respiration = sim$B * resprate # respiration flux [gWW.m^-2.yr^-1]
  
  # Get last 40% of timeseries
  etaTime <- 0.4 
  ixTime  <- which(sim$t >= ((1 - etaTime) * sim$t[sim$nTime]))
  
  sim$fluxCarcass   <- colMeans(sim$B[ixTime,] * p$mort0[-c(1:p$nResources)]) # deadfalls flux [gWW.m^-2.yr^-1]
  sim$fluxFecal <- colMeans(feces[ixTime,]) # total feces flux [gWW.m^-2.yr^-1]
  sim$fluxRepro  <- colMeans(rep2feces[ixTime,]) # feces flux from reproduction waste [gWW.m^-2.yr^-1]
  sim$fluxRespiration <- colMeans(respiration[ixTime,]) # respiration flux [gWW.m^-2.yr^-1]
  
  return(sim)
}

#
# Calculate injection fluxes calculated at all z-levels
# Units i gC/m3/year
#
calcCarbonInjection = function(sim) {
  rho_gC_gWW = 9 # Gram carbon per gram wet weight
  #
  # Calculate the injection from sinking POC with velocity v:
  #
  calcPOCinject = function(J, v) {
    
    Solve_Detritus_Euler <- function(z, alpha, zeta_X, v) {
      n <- length(z)
      dz <- c(diff(z), tail(diff(z), 1))
      
      DX <- numeric(n)
      DX[1] <- 0  # Boundary conditions at the surface -> no detritus
      
      for (i in 1:(n - 1)) {
        DX[i + 1] <- (DX[i] + dz[i] * zeta_X[i] / v) / (1 + dz[i] * alpha[i] / v)
      }
      
      return(DX)
    }
    
    # Multiply with the vertical probability distributions
    Jday <-  t(pDay) * J * 0.5 # Units: gWW/m^3/year
    Jnight <- t(pNight) * J * 0.5
    
    # remineralization at each depth
    z <- 0:p$bottom
    alpha <- rep(NA, length(z)) # bacterial degradation
    zeta_X <- rep(0, length(z)) # source = poc production at each depth
    
    # set a depth-dependent bacterial degradation
    alpha[z <= 100] <- rp
    alpha[z > 100 & z <= 1500] <- rm
    alpha[z > 1500] <- rb
    
    # sum total poc production
    Jtotal <- colSums(Jday + Jnight) # Units gWW/m3/year
    
    zeta_X[1:length(Jtotal)] <- Jtotal
    
    # calculate detritus production at each depth
    DX_Euler <- Solve_Detritus_Euler(z, alpha, zeta_X, v)
    # detritus remineralized
    Jinject = alpha * DX_Euler
    
    res = list()
    res$depth <- z
    res$Jbottom <- sum(Jtotal) - sum(Jinject) # gWW/m2/yr
    
    Jinject[length(Jinject)] <- Jinject[length(Jinject)] + res$Jbottom
    #Jinject_fecal[z < p$photic] <- 0 # carbon close to the surface is not sequestered
    
    res$inject <- Jinject
    
    return(res)
  }
  
  
  # parameters
  alpha0 = 0.25 # maximum bacterial degradation rate (from Pinti) [day^-1]
  Q10r = 2 # Q10 remin [-]
  
  p = sim$p
  functy = "all"
  
  # adjust bacterial degradation rate with temperature
  Tr=p$Tp # pelagic water temperature
  rfac = Q10r^((Tr-10)/10)
  rp = rfac*alpha0
  
  Tr=p$Tm # mid-water temperature
  rfac = Q10r^((Tr-10)/10)
  rm = rfac*alpha0
  
  Tr=p$Tb # bottom temperature
  rfac = Q10r^((Tr-10)/10)
  rb = rfac*alpha0
  
  # depth indices
  ix100 = ifelse(test = p$bottom >= 100, yes = 100 + 1, no = 0)
  ixbottom = p$bottom + 1
  
  col_indices <- switch(functy,
                        "smallPel" = 5:10,
                        "mesoPel" = 11:16,
                        "largePel" = 17:25,
                        "mwpred" = 26:34,
                        "dem" = 35:43,
                        "all" = 5:43,
                        stop("Unknown functy value : ", functy))
  #
  # Calculation the injection from respiration:
  #
  pDay = p$depthDay[, col_indices]  
  pNight = p$depthNight[, col_indices]
  respiration = t(pDay+pNight) * sim$fluxRespiration/2 # For each size class
  
  # poc production
  resFecal = calcPOCinject( as.vector(sim$fluxFecal), 500)
  resCarcass = calcPOCinject( as.vector(sim$fluxCarcass), 700)
  resRepro = calcPOCinject( as.vector(sim$fluxRepro), 500)
  
  # Make list with injections as output and convert to carbon units:
  res = list()
  res$z = resFecal$depth
  res$Fecal = resFecal$inject / rho_gC_gWW
  res$Carcass = resCarcass$inject / rho_gC_gWW
  res$Repro = resRepro$inject / rho_gC_gWW
  res$Respiration = colSums( respiration ) / rho_gC_gWW
  
  res$total = res$Fecal + res$Carcass + res$Repro + res$Respiration
  
  return(res)
}

#
# Simulate a given position in the global map:

### PERHAPS MOVE TO MAIN FEISTY (including data file)
#
simulatePosition = function(setup, lat, lon, nStages=9, tEnd=500) {
  # Output from COBALT
  glob <- read.csv("data/Cobalt global data.csv")
  
  ix = which.min( (glob$lat-lat)^2 + (glob$lon-lon)^2 ) # Find the best fitting location
  
  p = setup(
    szprod = glob[ix, "szprod"],        # small zooplankton production
    lzprod = glob[ix, "lzprod"],        # large zooplankton production
    dfbot  = glob[ix, "dfbot"],         # detrital flux reaching the bottom
    photic = glob[ix, "photic"],        # photic zone depth
    depth  = glob[ix, "depth"],         # water column depth
    Tp     = glob[ix, "Tp"],            # pelagic water temperature
    Tm     = glob[ix, "Tm"],            # mid-water temperature
    Tb     = glob[ix, "Tb"],            # bottom water temperature
    nStages = nStages                   # size class number
  )
  
  sim = simulateFEISTY(p = p, tEnd = tEnd)
  return(sim)
}

loadTransportMatrix = function(sFilename="data/CTL.R") {
  #
  # Transport matrix: 
  #
  TM <- readMat(sFilename, sparseMatrixClass="Matrix")
  TM <- CTL.mat$output
  
  ## When loaded the names of every variables were missing in the CTL.mat object
  # missing variable names in "CTL.mat"
  var_names <- dimnames(TM)[[1]]
  TM <- setNames(as.list(TM[,1,1]), var_names)
  
  # missing variable names in "msk"
  msk_names <- dimnames(TM$msk)[[1]]
  TM$msk <- setNames(as.list(TM$msk[,1,1]), msk_names)
  
  # missing variable names in "grid"
  grid_names <- dimnames(TM$grid)[[1]]
  TM$grid <- setNames(as.list(TM$grid[,1,1]), grid_names)
  
  # missing variable names in "MSKS"
  MSKS_names <- dimnames(TM$MSKS)[[1]]
  TM$MSKS <- setNames(as.list(TM$MSKS[,1,1]), MSKS_names)
  
  return(TM)
}
#
# Project the injection calculations onto the TM grid by integrating
# over the entire vertical cell
#
# gC/m2/yr
#
project_injection_to_TM <- function(inject, lat,lon, TM) {
  integral = 0*unique(TM$grid$zt)
  # Find closest grid point:
  ix = list( 
    y=which.min( (lat-TM$grid$yt)^2 ),
    x= which.min( (lon-TM$grid$xt)^2 ))
  # Integrate along the depth:
  for (j in 1:length(TM$grid$zt)) {
      idx = ( (inject$z > TM$grid$zw[j]) 
            & (inject$z <= (TM$grid$zw[j] + TM$grid$dzt[j])))
      integral[j] = trapz( inject$z[idx], inject$total[idx])
      
    }
  return(list(inject=integral, ix=ix))
}




testCarbonCalculations = function() {
  lat = 60
  lon = -15
  # Simulate the position 60, -15 using Cobalt output:
  sim = simulatePosition(setupVertical2, lat, lon)
  
  # Calculate carbon fluxes at the position of the fish:
  sim = calcCarbonFluxes(sim) 
  totalFlux = sim$fluxCarcass + sim$fluxFecal + sim$fluxRepro + sim$fluxRespiration
  barplot(totalFlux, xlab="Size class", ylab="Flux (gWW/m2/yr)")
  
  # Calculate the injection
  inject = calcCarbonInjection(sim)
  
  z = -inject$z
  plot( inject$total, z, type="l", lwd=3, 
        xlim=c(0,max(inject$total[1:length(inject$total)-1])),
        xlab="Injection (gC/m3/yr)",
        ylab="Depth (m)")
  lines( inject$Fecal, z, col="brown" )
  lines( inject$Carcass, z, col="grey" )
  lines( inject$Repro, z, col="darkgreen")
  lines( inject$Respiration, z, col="darkred")
  legend("top",
         c("Total","Fecal pellets","Carcasses","Reproduction","Respiration"),
         lwd=c(3,1,1,1,1),
         col=c("black","brown","grey","darkgreen","darkred")
  )
  
  # Calculate injection on TM grid:
  #TM = loadTransportMatrix()
  injectTM = project_injection_to_TM(inject, lat,lon, TM) 

}

