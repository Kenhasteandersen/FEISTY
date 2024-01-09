#
#Nov 2023
#

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
#' @param sim The FEISTY simulation result list from \code{\link{simulateFEISTY}}.
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
#' \item yieldAMean: a vector containing the arithmetic mean yield data [gww/m2/year] of each functional type of the time range specified.
#' \item yieldGMean: a vector containing the natural log-based geometric mean yield data [gww/m2/year] of each functional type of the time range specified.
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
  yieldAMean = rep(data=0, p$nGroups)
  yieldGMean = yieldAMean
  yieldMin = yieldAMean
  yieldMax = yieldAMean
  
  ixTime = which(sim$t>=((1-etaTime)*sim$t[sim$nTime]))
  
  for (iGroup in 1:p$nGroups) {
    ix = p$ix[[iGroup]]    
    yieldAllgrid[,ix-length(p$ixR)] =t(t(sim$u[, ix]) * p$mortF[ix]) 
    yieldAllgrid[,ix-length(p$ixR)][yieldAllgrid[,ix-length(p$ixR)]<0]=0
    yield[,iGroup]= rowSums(yieldAllgrid[,ix-length(p$ixR)])

    #deltaM = p$mUpper[ix]-p$mLower[ix]
    yieldAMean[iGroup] = mean(rowSums(yieldAllgrid[ixTime,ix-max(p$ixR)]))
    yieldGMean[iGroup] = exp(mean(log(rowSums(yieldAllgrid[ixTime,ix-max(p$ixR)]))))
    yieldMin[iGroup] = min(rowSums( yieldAllgrid[ixTime,ix-max(p$ixR)] ))
    yieldMax[iGroup] = max(rowSums( yieldAllgrid[ixTime,ix-max(p$ixR)] ))
  }

  sim$yieldAMean=yieldAMean
  sim$yieldGMean=yieldGMean  
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
#' @param sim The FEISTY simulation result list from \code{\link{simulateFEISTY}}.
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
#' \item SSBAMean: a vector containing the arithmetic mean SSB data [gww/m2/year] of each functional type of the time range specified.
#' \item SSBGMean: a vector containing the natural log-based geometric mean SSB data [gww/m2/year] of each functional type of the time range specified.
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
  SSBAMean = rep(data=0, p$nGroups)
  SSBGMean = SSBAMean
  SSBMin = SSBAMean
  SSBMax = SSBAMean
  
  ixTime = which(sim$t>=((1-etaTime)*sim$t[sim$nTime]))
      
  for (iGroup in 1:p$nGroups) {
    ix = p$ix[[iGroup]]
    SSBAllgrid[,ix-length(p$ixR)] =t(t(sim$u[, ix]) * p$psiMature[ix]) 
    SSBAllgrid[,ix-length(p$ixR)][SSBAllgrid[,ix-length(p$ixR)]<0]=0
    SSB[,iGroup]= rowSums(SSBAllgrid[,ix-length(p$ixR)])
    #deltaM = p$mUpper[ix]-p$mLower[ix]
    
    SSBAMean[iGroup] = mean(rowSums(SSBAllgrid[ixTime,ix-max(p$ixR)]))
    SSBGMean[iGroup] = exp(mean(log(rowSums(SSBAllgrid[ixTime,ix-max(p$ixR)]))))
    SSBMin[iGroup] = min(rowSums( SSBAllgrid[ixTime,ix-max(p$ixR)] ))
    SSBMax[iGroup] = max(rowSums( SSBAllgrid[ixTime,ix-max(p$ixR)] ))
  }
  
  sim$SSBAMean=SSBAMean
  sim$SSBGMean=SSBGMean
  sim$SSBMin=SSBMin
  sim$SSBMax=SSBMax
  sim$SSB=SSB

  return(sim)
}

