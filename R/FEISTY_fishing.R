#' Set Fishing Mortality
#'
#' This function sets fishing mortality for specific functional types based on maximum fishing mortality and size-based trawl selectivity.
#'
#' @param p The model parameter list, as created with e.g. \link{setupBasic2}.
#' @param Fmax The baseline fishing mortality rate [1/year]. Default 0, indicating no fishing.
#'   May be a single value applied to every group in \code{groupidx}, or a vector
#'   of the same length as \code{groupidx} giving a group-specific rate per element
#'   (\code{Fmax[i]} is applied to group \code{groupidx[i]}).
#' @param etaF A coefficient determining the fish size with 50\% fishing selectivity. The value represents the fraction of the maximum size of a fish functional type. The default value is 0.05.
#'   May be a single value applied to every group in \code{groupidx}, or a vector
#'   of the same length as \code{groupidx} giving a group-specific value per element
#'   (\code{etaF[i]} is applied to group \code{groupidx[i]}).
#' @param groupidx An integer vector containing indices of functional types that fishing mortality will be assigned to. Default is all functional types.
#'
#' @return It returns an updated parameter list:
#' \itemize{
#' \item Fmax, group-level maximum fishing mortality.
#' \item etaF, group-level selectivity coefficient.
#' \item mortF, a vector containing fishing mortality of all state variables, including resources (always 0) and fish.
#' }
#'
#' @details The function sets fishing mortality for the functional groups specified by \code{groupidx}.
#' For each specified group, it calculates the selectivity \code{psi} using the standard trawl selectivity formula from Andersen (2019) \bold{Fig 5.2}.
#' The fishing mortality \code{mortF} for specified groups is then updated based on the calculated selectivity \code{psi} and the baseline fishing mortality rate \code{Fmax}.
#' \code{Fmax} and \code{etaF} may each be a single value (applied to every group in
#' \code{groupidx}) or a vector of the same length as \code{groupidx}, in which case
#' element \code{i} is applied to group \code{groupidx[i]}, allowing group-specific
#' fishing mortality and selectivity in a single call. Full group-level vectors
#' of length \code{p$nGroups} are also accepted and subset by \code{groupidx}.
#' A length other than 1, \code{length(groupidx)}, or \code{p$nGroups} raises
#' an error.
#'
#' @examples
#' p = setupVertical2(Fmax=0) # No fishing mortality
#' p = setFishing(p, Fmax = 1, etaF = 0.05, groupidx=c(5)) # add fishing mortality to demersals only.
#' # group-specific fishing: different Fmax and etaF per group
#' p = setFishing(p, Fmax = c(0, 1, 0.2), etaF = c(0.05, 0.01, 0.1), groupidx = c(1, 2, 5))
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
setFishing = function(p, Fmax=0, etaF=0.05, groupidx=seq_len(p$nGroups)) {
  groupidx = as.integer(groupidx)
  if (any(is.na(groupidx)) || any(groupidx < 1) || any(groupidx > p$nGroups))
    stop("groupidx must contain valid functional-group indices.")

  nSel = length(groupidx)
  if (nSel == 0) return(p)

  # Accept full group-level vectors and subset them for this call.
  if (length(Fmax) == p$nGroups && nSel < p$nGroups) Fmax = Fmax[groupidx]
  if (length(etaF) == p$nGroups && nSel < p$nGroups) etaF = etaF[groupidx]

  # Recycle scalars to match groupidx (backwards compatible); validate vector lengths.
  if (length(Fmax) == 1) Fmax = rep(Fmax, nSel)
  if (length(etaF) == 1) etaF = rep(etaF, nSel)
  if (length(Fmax) != nSel)
    stop(sprintf("Fmax must be length 1, length(groupidx) (%d), or p$nGroups (%d); got %d",
                 nSel, p$nGroups, length(Fmax)))
  if (length(etaF) != nSel)
    stop(sprintf("etaF must be length 1, length(groupidx) (%d), or p$nGroups (%d); got %d",
                 nSel, p$nGroups, length(etaF)))

  # Store full group-level vectors so sequential group-specific calls do not
  # leave p$Fmax/p$etaF representing only the last group that was updated.
  oldFmax = p$Fmax
  oldEtaF = p$etaF
  if (is.null(oldFmax) || length(oldFmax) != p$nGroups)
    oldFmax = rep(if (length(oldFmax) == 1) oldFmax else 0, p$nGroups)
  if (is.null(oldEtaF) || length(oldEtaF) != p$nGroups)
    oldEtaF = rep(if (length(oldEtaF) == 1) oldEtaF else 0.05, p$nGroups)

  p$Fmax = oldFmax
  p$etaF = oldEtaF
  p$Fmax[groupidx] = Fmax
  p$etaF[groupidx] = etaF
  for (iGroup in 1:nSel) {
    ix = p$ix[[groupidx[iGroup]]]
    mFishing = p$etaF[groupidx[iGroup]]*max(p$mUpper[ix]) # selectivity at etaF of maximum size
    psi = ( 1 + (p$mc[ix]/mFishing)^(-3) )^(-1) # Standard trawl selectivity from Andersen (2019) Fig 5.2
    p$mortF[ix] = psi*p$Fmax[groupidx[iGroup]]
  }
  return(p)
}

#' Yield Calculation
#'
#' This function calculates the yield for each functional type in units of g/m2/year based on a FEISTY simulation result.
#'
#' @usage calcYield(sim, etaTime = 0.4)
#'
#' @param sim The FEISTY simulation result list from \code{\link{simulateFEISTY}}.
#' @param etaTime The last fraction of the simulation period (default the last 40\%).
#'
#' @details
#' This function calculates the yield for each function type based on the FEISTY simulation results within the specified time fraction (default is the last 40\% of the simulation period). \cr
#' Yield is the product of biomass \code{u} [g/m2] and fishing mortality \code{mortF} [1/year]. Negative values in \code{u} are corrected to 0. 
#' The yield calculation is automatically adaptive to non-time-series simulations or time-series simulations. \cr
#' This function is called automatically inside \code{\link{simulateFEISTY}} (with \code{etaTime = 0.4}),
#' so its outputs are already present in a standard simulation result. It can also be called directly on a
#' simulation result, for example to recompute the yield over a different time fraction \code{etaTime}.
#'
#' @return 
#' Add yield data to the result list:
#' \itemize{
#' \item yieldAMean: a vector containing the arithmetic mean yield data [g/m2/year] of each functional type of the time range specified.
#' \item yieldGMean: a vector containing the natural log-based geometric mean yield data [g/m2/year] of each functional type of the time range specified.
#' \item yieldMin: a vector containing the minimum yield data [g/m2/year] of each functional type within the time range specified.
#' \item yieldMax: a vector containing the maximum yield data [g/m2/year] of each functional type within the time range specified.
#' \item yield: a matrix containing the yield data [g/m2/year] of each functional type (column) in each time point (row)
#' }
#'
# @examples
#'
#' @author Yixin Zhao
#' 
#' @seealso 
#' \code{\link{setFishing}} 	Set fishing mortality
#' 
#' @aliases calcYield
#'
#' @export
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
    if(is.null(dim(p$mortF))){ # if it is not a matrix (non-ts)
    mortF_mat = matrix(data=rep(p$mortF[ix], each = sim$nTime), 
                       nrow = sim$nTime, 
                       ncol = length(p$mortF[ix]), 
                       byrow = FALSE)
  }else{
    mortF_mat = p$mortF[,ix]
    }
    
    #yieldAllgrid[,ix-length(p$ixR)] =t(t(sim$u[, ix]) * p$mortF[ix]) 
    yieldAllgrid[,ix-length(p$ixR)] = sim$u[, ix] * mortF_mat
    yieldAllgrid[,ix-length(p$ixR)][yieldAllgrid[,ix-length(p$ixR)]<0]=0
    yield[,iGroup]= rowSums(yieldAllgrid[,ix-length(p$ixR)])
    
    #deltaM = p$mUpper[ix]-p$mLower[ix]
    yieldAMean[iGroup] = mean(rowSums(yieldAllgrid[ixTime,ix-max(p$ixR),drop=FALSE]))
    yieldGMean[iGroup] = exp(mean(log(rowSums(yieldAllgrid[ixTime,ix-max(p$ixR),drop=FALSE]))))
    yieldMin[iGroup] = min(rowSums( yieldAllgrid[ixTime,ix-max(p$ixR),drop=FALSE] ))
    yieldMax[iGroup] = max(rowSums( yieldAllgrid[ixTime,ix-max(p$ixR),drop=FALSE] ))
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
#' This function calculates the spawning stock biomass (SSB) for each functional 
#' type based on a FEISTY simulation result in units of g/m2.
#'
#' @usage calcSSB(sim, etaTime = 0.4)
#'
#' @param sim The FEISTY simulation result list from \code{\link{simulateFEISTY}}.
#' @param etaTime The last fraction of the simulation period (default the last 40\%).
#'
#' @details
#' This function calculates the spawning stock biomass for each function type based on the FEISTY simulation results within the specified time fraction (default is the last 40\% of the simulation period). \cr
#' Spawning stock biomass is the product of biomass \code{u} [g/m2] and maturity level \code{psiMature} [dimensionless]. Negative values in \code{u} are corrected to 0. \cr
#' Note the SSB data only represents how much biomass (energy) could be used for reproduction, rather than the amount of offspring. \cr
#' This function is called automatically inside \code{\link{simulateFEISTY}} (with \code{etaTime = 0.4}),
#' so its outputs are already present in a standard simulation result. It can also be called directly on a
#' simulation result, for example to recompute the SSB over a different time fraction \code{etaTime}.
#'
#' @return 
#' Add SSB data to the result list:
#' \itemize{
#' \item SSBAMean: a vector containing the arithmetic mean SSB data [g/m2] of each functional type of the time range specified.
#' \item SSBGMean: a vector containing the natural log-based geometric mean SSB data [g/m2] of each functional type of the time range specified.
#' \item SSBMin: a vector containing the minimum SSB data [g/m2] of each functional type within the time range specified.
#' \item SSBMax: a vector containing the maximum SSB data [g/m2] of each functional type within the time range specified.
#' \item SSB: a matrix containing the SSB data [g/m2] of each functional type (column) in each time point (row)
#' }
#'
# @examples
#' 
#'
#' @author Yixin Zhao
#' 
#' @aliases calcSSB
#' 
#' @seealso 
#' \code{\link{paramAddGroup}} 	Add parameters of one functional type
#'
#' @export
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
    
    SSBAMean[iGroup] = mean(rowSums(SSBAllgrid[ixTime,ix-max(p$ixR),drop=FALSE]))
    SSBGMean[iGroup] = exp(mean(log(rowSums(SSBAllgrid[ixTime,ix-max(p$ixR),drop=FALSE]))))
    SSBMin[iGroup] = min(rowSums( SSBAllgrid[ixTime,ix-max(p$ixR),drop=FALSE] ))
    SSBMax[iGroup] = max(rowSums( SSBAllgrid[ixTime,ix-max(p$ixR),drop=FALSE] ))
  }
  
  sim$SSBAMean=SSBAMean
  sim$SSBGMean=SSBGMean
  sim$SSBMin=SSBMin
  sim$SSBMax=SSBMax
  sim$SSB=SSB
  
  return(sim)
}

