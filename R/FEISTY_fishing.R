#
#Nov 2023
#

#
# Analysis of the fishing mortality
#

analyseStages = function(nStages = c(3,6,9,12,15,18,21,24,27), maxF = 10) {
  y = list()
  for (i in 1:length(nStages)) {
    y[[i]] = plotYield(p=setupBasic2(szprod = 100,lzprod = 100, bprod  = 5,
                                     depth  = 100,
                                     Tp = 10,
                                     Tb = 8, 
                                     nStages=nStages[i],
                                     etaMature=0.25,
                                     F=0, # overwritten later
                                     etaF=0.05),
                                     maxF=maxF)
  }
  F = y[[1]][[1]]
  
  defaultplot()
  defaultpanel(xlim=F, ylim=c(0,max(sapply(y, function(x) x[[2]]))+5),
               xlab="Fishing mortality (year$^{-1}$)",
               ylab="Yield (g/$m^2$/yr)")
  for (i in 1:length(nStages)) {
    lines(F, y[[i]][[2]], lwd=i)
  }
  
  defaultplot()
  defaultpanel(xlim=F, ylim=c(0,max(sapply(y, function(x) x[[4]]))+5),
               xlab="Fishing mortality (year$^{-1}$)",
               ylab="SSB (g/$m^2$/yr)")
  for (i in 1:length(nStages)) {
    lines(F, y[[i]][[5]], lwd=i)
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
    sim = simulateFEISTY(bCust    = FALSE,
                         p, 
                         #tEnd   = 100,
                         #tStep  = 1,
                         #times  = seq(from=0, to=tEnd, by=tStep),  
                         yini   = p$u0,  
                         USEdll = TRUE,
                         Rmodel = derivativesFEISTYR,
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
  
  return(list(F,yield, yieldMin, yieldMax,SSB))
}

#' Set Fishing Mortality
#'
#' This function sets fishing mortality for all fish based on a baseline fishing mortality and the size-based trawl selectivity.
#'
#' @param p The model parameter list, as created with e.g. \link{setupBasic2}, needed to be updated.
#' @param F The baseline fishing mortality rate [1/year]. Default 0, indicating no fishing.
#' @param etaF A coefficient determining the fish size with 50\% fishing selectivity. The value represents the fraction of the maximum size of a fish functional type. The default value is 0.05.
#' 
#' @return It returns an updated parameter list:
#' \itemize{
#' \item F, baseline fishing mortality, from parameter input.
#' \item etaF, from parameter input.
#' \item mortF, a vector containing fishing mortality of all state variables, including resources (always 0) and fish.
#' }
#' 
#' @details The function sets fishing mortality for all fish.
#' For each group, it calculates the selectivity \code{psi} using the standard trawl selectivity formula from Andersen (2019) \bold{Fig 5.2}. 
#' The fishing mortality \code{mortF} for each group is then updated based on the calculated selectivity \code{psi} and the baseline fishing mortality rate \code{F}.
#'
#' @examples
#' p = setupBasic2(F=0) # No fishing mortality
#' p = setFishing(p, F = 1, etaF = 0.05) # add fishing mortality
#'
#' @references
#' Andersen, K. H. (2019). Fish ecology, evolution, and exploitation: a new theoretical synthesis. Princeton University Press.
#' 
#' @author Yixin Zhao
#' 
#' @aliases setFishing
#' 
#' @seealso 
#' \code{\link{calcYield}} Yield calculation
#' 
#' @export

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

#' Yield Calculation
#'
#' This function calculates the yield for each function type based on a FEISTY simulation result.
#'
#' @usage calcYield(sim, etaTime = 0.4)
#'
#' @param sim The FEISTY simulation result list from \code{\link{simulateFEISTY}}. Note `simpleOutput` must be `TRUE` e.g., \code{sim=simulateFEISTY(simpleOutput=TRUE)}.
#' @param etaTime The last fraction of the simulation period (default the last 40\%).
#'
#' @details
#' This function calculates the yield for each function type based on the FEISTY simulation results within the specified time fraction (default is the last 40\% of the simulation period). \cr
#' Yield is the product of biomass \code{u} [gww/m2] and fishing mortality \code{mortF} [1/year]. Negative values in \code{u} are revised to 0. \cr
#' This function has been integrated into \code{\link{simulateFEISTY}}. It cannot be called independently.
#'
#' @return 
#' Add yield data to the result list:
#' \itemize{
#' \item yieldMean: a vector containing the mean yield data [gww/m2/year] of each functional type of the time range specified.
#' \item yieldMin: a vector containing the minimum yield data [gww/m2/year] of each functional type within the time range specified.
#' \item yieldMax: a vector containing the maximum yield data [gww/m2/year] of each functional type within the time range specified.
#' \item yield: a matrix containing the yield data [gww/m2/year] of each functional type (column) in each time point (row)
#' }
#'
#' @examples
#' no examples
#'
#' @author Yixin Zhao
#' 
#' @seealso 
#' \code{\link{setFishing}} 	Set fishing mortality
#' 
#' @aliases calcYield
#'
# @export

#
# Return the yield of all function groups
#
calcYield = function(
    sim,          # The simulation object to analyse
    etaTime=0.4) {# The last fraction of the simulation period (default the last 40%)
  
  p=sim$p
  
  yieldAllgrid = matrix(nrow=sim$nTime, ncol=length(p$ixFish))
  yield = matrix(nrow=sim$nTime, ncol=p$nGroups)
  yieldMean = rep(data=0, p$nGroups)
  yieldMin = yieldMean
  yieldMax = yieldMean
  
  ixTime = which(sim$t>=((1-etaTime)*sim$t[sim$nTime]))
  
  for (iGroup in 1:p$nGroups) {
    ix = p$ix[[iGroup]]    
    yieldAllgrid[,ix-length(p$ixR)] =t(t(sim$u[, ix]) * p$mortF[ix]) 
    yieldAllgrid[,ix-length(p$ixR)][yieldAllgrid[,ix-length(p$ixR)]<0]=0
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

#' Spawning Stock Biomass Calculation
#'
#' This function calculates the spawning stock biomass (SSB) for each function type based on a FEISTY simulation result.
#'
#' @usage calcSSB(sim, etaTime = 0.4)
#'
#' @param sim The FEISTY simulation result list from \code{\link{simulateFEISTY}}. Note `simpleOutput` must be `TRUE` e.g., \code{sim=simulateFEISTY(simpleOutput=TRUE)}.
#' @param etaTime The last fraction of the simulation period (default the last 40\%).
#'
#' @details
#' This function calculates the spawning stock biomass for each function type based on the FEISTY simulation results within the specified time fraction (default is the last 40\% of the simulation period). \cr
#' Spawning stock biomass is the product of biomass \code{u} [gww/m2] and maturity level \code{psiMature} [1/year]. Negative values in \code{u} are revised to 0. \cr
#' Note the SSB data only represents how much biomass (energy) could be used for reproduction, rather than the amount of offspring. \cr
#' This function has been integrated into \code{\link{simulateFEISTY}}. It cannot be called independently.
#'
#' @return 
#' Add SSB data to the result list:
#' \itemize{
#' \item SSBMean: a vector containing the mean SSB data [gww/m2/year] of each functional type of the time range specified.
#' \item SSBMin: a vector containing the minimum SSB data [gww/m2/year] of each functional type within the time range specified.
#' \item SSBMax: a vector containing the maximum SSB data [gww/m2/year] of each functional type within the time range specified.
#' \item SSB: a matrix containing the SSB data [gww/m2/year] of each functional type (column) in each time point (row)
#' }
#'
#' @examples
#' no examples
#'
#' @author Yixin Zhao
#' 
#' @aliases calcSSB
#' 
#' @seealso 
#' \code{\link{paramAddGroup}} 	Add parameters of one functional type
#'
# @export

#
# Return the SSB of all function groups
#
calcSSB = function(
    sim,          # The simulation object to analyse
    etaTime=0.4) {# The last fraction of the simulation period (default the last 40%) 
  
  p=sim$p
  
  SSBAllgrid = matrix(nrow=sim$nTime, ncol=length(p$ixFish))
  SSB = matrix(nrow=sim$nTime, ncol=p$nGroups)
  SSBMean = rep(data=0, p$nGroups)
  SSBMin = SSBMean
  SSBMax = SSBMean
  
  ixTime = which(sim$t>=((1-etaTime)*sim$t[sim$nTime]))
      
  for (iGroup in 1:p$nGroups) {
    ix = p$ix[[iGroup]]
    SSBAllgrid[,ix-length(p$ixR)] =t(t(sim$u[, ix]) * p$psiMature[ix]) 
    SSBAllgrid[,ix-length(p$ixR)][SSBAllgrid[,ix-length(p$ixR)]<0]=0
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

# Analyse sensitivity for fishing mortality. x: Stage number y: Biomass
analyseStageB = function(nStages = c(3,6,9,12,15,18,21,24,27),maxF=20) {
  y = list()
  F= seq(0,maxF,length.out=21)
 
for (iF in 1:length(F)){  
  
  for (i in 1:length(nStages)) {
     sim = simulateFEISTY(p=setupBasic2(szprod = 100,lzprod = 100, bprod  = 5,
                                     depth  = 100,
                                     Tp = 10,
                                     Tb = 8, 
                                     nStages=nStages[i],
                                     etaMature=0.25,
                                     F=F[iF], # overwritten later
                                     etaF=0.05),
                          simpleOutput = TRUE)
     p=sim$p
     Bpositive=sim$B
     Bpositive[Bpositive<0]=0
# all functional types
     for(j in 1:p$nGroups){
     y[[paste(i)]][[j]]= sum(colMeans(Bpositive[round(0.6*sim$nTime):sim$nTime,p$ix[[j]]-p$nResources])) 
     }
     
# tot B     
     y[[paste(i)]][["totB"]]= sum(colMeans(Bpositive[round(0.6*sim$nTime):sim$nTime,p$ixFish-p$nResources])) 
}
  
  defaultplot()
  defaultpanel(xlim=nStages, ylim=c(0,max(unlist(y))+5),
               xlab="Stage number",
               ylab="Biomass (gww m $^{-2}$)")
  for (i in 1:p$nGroups) {
    lines(nStages, unlist(sapply(y, function(x) x[i])), lwd=i, type = "o", pch = 1)
  }
  
  lines(nStages, unlist(sapply(y, function(x) x[["totB"]])), lwd=2, col="red", type = "o", pch = 1)
  
  title(main=paste("Fishing mortality", F[iF], "(year-1)"),line = -1)
  
 }
  
}



##### Analyse sensitivity for size preference function. x: Stage number y: Biomass

# edit setupBasic2 (,sptype=isp)
# param$theta=sizePrefFeeding(p=param,           # parameter settings
#                             beta = 400,  # preferred predator/prey mass ratio
#                             sigma = 1.3, # width of size preference for feeding
#                             type = sptype)
# 
# # update preference from fish to resources with the simple function
# for (i in param$ixFish) {
#   for(j in param$ixR){
#     param$theta[i,j] = exp( -(log(param$mc[i]/(beta*param$mc[j])))^2 / (2*sigma)^2  )
#     if (param$mc[j] >param$mc[i]) param$theta[i, j] = 0
#   }
# }
#####
analyseSizepref = function(nStages = c(3,6,9,12,15,18,21,24,27),maxF=0) {
  y = list()
  F= seq(0,maxF,length.out=1)

  for (isp in 1:3){  
    for (i in 1:length(nStages)) {
      # sim = simulateFEISTY(bCust=TRUE,p=setupBasic2(szprod = 100,lzprod = 100, bprod  = 5,
      #                                    depth  = 100,
      #                                    Tp = 10,
      #                                    Tb = 8, 
      #                                    nStages=nStages[i],
      #                                    etaMature=0.25,
      #                                    F=F, # overwritten later
      #                                    etaF=0.05  ,sptype=isp),
      #                      simpleOutput = TRUE)
      sim = simulateFEISTY(bCust=TRUE,p=setupVertical2(szprod= 80,lzprod = 80, # Pelagic productivities
                                                       bent = 150, # Detrital flux out of photic zone
                                                       nStages=nStages[i], # No. of size groups
                                                       region = 1, # Temperature profile regions: 1 Tropical, 2 Temperate, 3 Boreal, 4 Default 10 Celsius 
                                                       depth=1500, # Bottom depth
                                                       photic=150, # Photic zone depth
                                                       mesopelagic=250, # mesopelagic depth
                                                       visual=1.5,# >1 visual predation primarily during the day, = 1 equal day and night
                                                       etaMature = 0.25, # Size of matureation relative to
                                                       # asymptotic size. Different from
                                                       # van Denderen (2021), where it is 0.002
                                                       F=0,
                                                       etaF=0.05,sptype=isp),
                           simpleOutput = TRUE)
      p=sim$p
      Bpositive=sim$B
      Bpositive[Bpositive<0]=0
      # all functional types
      for(j in 1:p$nGroups){
        y[[paste(isp)]][[paste(i)]][[j]]= sum(colMeans(Bpositive[round(0.6*sim$nTime):sim$nTime,p$ix[[j]]-p$nResources])) 
      }
      
      # tot B     
      y[[paste(isp)]][[paste(i)]][["totB"]]= sum(colMeans(Bpositive[round(0.6*sim$nTime):sim$nTime,p$ixFish-p$nResources])) 
    }
    

  }
  cc=10+10
  defaultplot()
  defaultpanel(xlim=nStages, ylim=c(0,max(unlist(y))+5),
               xlab="Stage number",
               ylab="Biomass (gww m $^{-2}$)")
  # for (i in 1:p$nGroups) {
  #   lines(nStages, unlist(sapply(y, function(x) x[i])), lwd=i, type = "o", pch = 1)
  # }
  for (isp in 1:3) {
    lines(nStages, unlist(sapply(y[[paste(isp)]], function(x) x[["totB"]])), lwd=isp, col="red", type = "o", pch = 1)
  }
  #title(main=paste("Fishing mortality", F[iF], "(year-1)"),line = -1)
  return(y)
}
