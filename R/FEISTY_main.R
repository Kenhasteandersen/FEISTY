#===============================================================================
# The FEISTY model as described in 
# Petrik, C. M., Stock, C. A., Andersen, K. H., van Denderen, P. D., & Watson, J. R. (2019). 
# Bottom-up drivers of global patterns of demersal, forage, and pelagic fishes. 
# Progress in oceanography, 176, 102124.

# van Denderen, P. D., Petrik, C. M., Stock, C. A., & Andersen, K. H. (2021). 
# Emergent global biogeography of marine fish food webs. 
# Global Ecology and Biogeography, 30(9), 1822-1834.

# 
# Main model routines
#      derivativesFEISTYR  : derivative code in R
#      simulateFEISTY      : run the FEISTY model
#
#===============================================================================

#' Calculate Derivatives of State Variables in FEISTY Model in R
#'
#' @description
#' This function calculates the derivatives of state variables in FEISTY model in R. 
#' It is used for time integration of FEISTY simulations or getting rates by computing derivatives of one time step.
#' Note that calculating derivatives in R is much slower than with Fortran.
#'
#' @usage derivativesFEISTYR(t, u, p, FullOutput=TRUE)
#'
#' @param t Current time. Any numeric value, no use currently.
#' @param u A numeric vector containing all state variables, all resources + all size classes.
#' @param p The parameter list.
#' @param FullOutput Logical flag indicating whether to return all calculation results. Default TRUE.
#'
#' @details
#' This function is designed to simulate the dynamics of marine resources and fish populations in FEISTY model in R. 
#' The derivative calculation of resources can be based on different types of dynamics (chemostat or logistic).
#' Negative state variables and NA values are revised to zero to prevent computational errors.
#' ODE solving is conducted by package deSolve. The equations are documented in the vignette.
#' More information can be found in de Roos et al. (2008) and Petrik et al. (2019).
#'
#' @return
#' If \code{FullOutput} is 'TRUE', a list containing the following components:
#' \itemize{
#' \item deriv: a vector of derivatives [g/m2/year] of all resources and all functional type size classes.
#' \item f: a vector containing feeding levels [-] of all resources (0) and all size classes of functional types.
#' \item mortpred: a vector containing predation mortality rate [1/year] of all resources and all size classes of functional types.
#' \item g: Net growth rate (the fraction of available energy invested in growth) [1/year]. It includes all size classes of functional types. Resources not included.
#' \item Repro: Energy used for reproduction of all size classes of functional types [g/m2/year]. Resources not included.
#' \item Fin: Biomass flux into each size class [g/m2/year]. Resources not included.
#' \item Fout: Biomass flux out of each size class [g/m2/year]. Resources not included.
#' \item totMort: a vector containing total mortality [g/m2/year] of each functional type, 
#' which includes predation mortality, background mortality, and fishing mortality.
#' \item totGrazing: a vector containing total grazing (food intake before assimilation) [g/m2/year] of each functional type. Cmax * f * u (maximum consumption rate * feeding level * biomass).
#' \item totLoss: a vector containing total biomass loss [g/m2/year] of each functional type, including unassimilated food, basal metabolism and reproduction cost (1-epsRepro). 
#' The reproduction cost here also include the energy loss (1-epsRepro) from the flux out of the last size class (invested in reproduction) of each functional type.
#' These losses are supposed to be released to environments.
#' \item totRepro: a vector containing total energy used for reproduction [g/m2/year] of each functional type (before multiplying the reproduction efficiency epsRepro).
#' \item totRecruit: a vector containing total recruitment [g/m2/year] of each functional type. 
#' totRecruit = totRepro * epsRepro (reproduction efficiency)
#' \item totBiomass: a vector containing total biomass [g/m2] of each functional type. 
#' }
#' If \code{FullOutput = FALSE}, it only returns derivatives of resources and all size classes.
#'
#' @examples
#' # Parameter list preparation
#' p = setupBasic2()
#' # FEISTY simulation 200 years, use derivativesFEISTYR
#' sim = simulateFEISTY(p=p, tEnd=200, USEdll=FALSE)
#' # get rates of last year, based on biomass of last year 
#' rates = derivativesFEISTYR(t=0, p=p, u=sim$u[sim$nTime,], FullOutput=TRUE)
#' 
#' @references
#' Petrik, C. M., Stock, C. A., Andersen, K. H., van Denderen, P. D., & Watson, J. R. (2019). Bottom-up drivers of global patterns of demersal, forage, and pelagic fishes. Progress in oceanography, 176, 102124.
#' 
#' van Denderen, P. D., Petrik, C. M., Stock, C. A., & Andersen, K. H. (2021). Emergent global biogeography of marine fish food webs. Global Ecology and Biogeography, 30(9), 1822-1834.
#' 
#' de Roos, A. M., Schellekens, T., Van Kooten, T., Van De Wolfshaar, K., Claessen, D., & Persson, L. (2008). Simplifying a physiologically structured population model to a stage-structured biomass model. Theoretical population biology, 73(1), 47-62.
#' 
#' @author Ken H. Andersen, Karline Soetaert, Yixin Zhao
#'
#' @aliases derivativesFEISTYR
#' 
#' @seealso 
#' \code{\link{simulateFEISTY}} The main function for running FEISTY simulations
#'
#' @export
#' 

# ------------------------------------------------------------------------------
#
# Calculate the derivatives of all state variables in R
#
# In:
#  p : a list of parameters
#  u : all state variables
#  bFullOutput : if TRUE returns all internal calculations
#
# Out:
#  deriv : a vector of derivatives
#
# ------------------------------------------------------------------------------

derivativesFEISTYR = function(t,              # current time
                              u,              # all state variables
                              p,              # parameters
                              FullOutput=TRUE) {
  
  # split state variable vector into resource and fish
  u[u<0]=0
  R     = u[p$ixR]       # resource, prey
  iFish = p$ixFish
  B     = u[iFish]       # fish
  
  # update effective temperature for large demersal fish in shallow water
  if(!is.null(p$depth) & !is.null(p$bET))
    if (p$depth<200 & p$bET==TRUE) p=updateET(p=p,u=u)
  
  # ----------------------------------------------
  # Consumption of all fish groups
  # ----------------------------------------------
  
  # V: clearance rate, (m2/g/yr) 
  # theta x u: prey available for consumption
  # Cmax: maximum consumption rate, /yr 
  
  Enc = p$V * (p$theta %*% u)  # /yr
  
  f   = Enc / (p$Cmax + Enc)   # Functional response
  f[is.na(f)] = 0
  
  # net growth rate, /yr
  Eavail  = p$epsAssim * p$Cmax * f - p$metabolism
  
  # ----------------------------------------------
  # Predation mortality, /yr:
  # = t(p$theta) %*% (f*p$Cmax/p$epsAssim*u/p$mc)
  # ----------------------------------------------
  #
  mm = p$Cmax*p$V/(Enc+p$Cmax)*u # temporarily store
  mm[ is.na(mm) ] = 0
  mortpred = t(p$theta) %*% mm
  
  # ----------------------------------------------
  # Total mortality (includes basal and fishing mortality)
  # ----------------------------------------------
  mort = mortpred + p$mort0 + p$mortF   # /year
  
  # ----------------------------------------------
  # Derivative of fish groups
  # ----------------------------------------------
  
  # Flux out of the size group
  #------------------------------
  v     = Eavail[iFish,]   # net growth rate
  vplus = pmax(v,0)
  
  # fraction available for growth
  kappa = 1 - p$psiMature[iFish]   
  g     = kappa*vplus
  
  # growth to the next stage
  gamma = (kappa*vplus - mort[iFish]) /
    (1 - (1/p$z[iFish])^(1-mort[iFish]/(kappa*vplus)) )
  
  gamma[kappa==0] = 0 # No growth of fully mature classes
  
  # goes out of stage (size group)
  Fout = gamma*B
  
  # Energy used for reproduction
  Repro = p$psiMature[iFish]*vplus*B
  
  # Flux into the size group
  #------------------------------
  Fin = Fout*0
  for (i in 1:p$nGroups) {
    ix = p$ix[[i]] - p$ix[[1]][1] +1                  # growth to
    ixPrev  = c(ix[length(ix)], ix[1:(length(ix)-1)]) # growth from
    Fin[ix] = Fout[ixPrev]
    
    # for reproduction: consider the reproduction success
    Fin[ix[1]] = p$epsRepro[i]*(Fin[ix[1]] + sum( Repro[ix] ))
  }
  
  # ----------------------------------------------
  # Assemble derivatives of fish:
  # ----------------------------------------------
  
  dBdt = Fin - Fout + (v - mort[p$ixFish])*B - Repro
  
  # ----------------------------------------------
  # Derivative of resources
  # ----------------------------------------------
  if (p$Rtype == 1)  # chemostat
    dRdt = p$r*(p$K-R) - mortpred[p$ixR]*R
  else               # logistic
    dRdt = p$r*R*(1-R/p$K) - mortpred[p$ixR]*R
  
  # ----------------------------------------------
  # Assemble output:
  # ----------------------------------------------
  if (FullOutput) { # Output everything
    out = list()
    out$deriv = c(dRdt, dBdt)
    out$f     = f[-p$ixR,] # Feeding level only all fish stages, no resources
    out$mortpred = mortpred[,]
    out$g     = g # net growth rate fish stages
    out$Repro = Repro
    out$Fin   = Fin
    out$Fout  = Fout
    
    # for the budget:
    grazing = p$Cmax * f         # grazing rate, /yr
    loss    = (1.-p$epsAssim) * grazing + p$metabolism # Energy loss to environments. Updated below.
    Reprofrac = (1-kappa)*vplus
    
    il <- NULL
    for (i in 1:length(p$ix)){
      # Add the waste energy in reproduction of each stages. Note it does not include the waste energy from last stage energy flux out, added in `totLoss` below.
      loss[p$ix[[i]]] = loss[p$ix[[i]]] + (1-p$epsRepro[i]) * Reprofrac[p$ix[[i]]-p$nResources] 
      il <- c(il, rep(i, times=length(p$ix[[i]])))
    }
    
    out$totMort    = tapply((mort   *u)[p$ixFish], INDEX=il, FUN=sum)
    out$totGrazing = tapply((grazing*u)[p$ixFish], INDEX=il, FUN=sum)
    # Add the waste energy in reproduction from flux out of the last stage of each functional type.
    out$totLoss    = tapply((loss   *u)[p$ixFish], INDEX=il, FUN=sum) +(1-p$epsRepro)*Fout[sapply(p$ix, tail, n = 1)-p$nResources] 
    out$totRepro   = tapply(Repro, INDEX=il, FUN=sum) + Fout[sapply(p$ix, tail, n = 1)-p$nResources] 
    out$totRecruit = out$totRepro* p$epsRepro
    out$totBiomass = tapply(B, INDEX=il, FUN=sum)
    return(out)
  }
  else # Output just the derivatives
    return( list(c(dRdt, dBdt)) )
}

#' Run FEISTY model simulations
#'
#' @description
#' \code{simulateFEISTY} runs simulations of the FEISTY model to resolve the dynamics of marine resources and fish populations over a specified time frame. \cr
#' It provides options for integrating ordinary differential equations in Fortran or R for four prepared setups (\code{setupBasic}, \code{setupBasic2}, \code{setupVertical}, and \code{setupVertical2}).
#'
#' @usage simulateFEISTY (p = setupBasic(), 
#'                        tEnd = 500, tStep  = 1, times = seq(from=0, to=tEnd, by=tStep), 
#'                        yini = p$u0, USEdll = TRUE, Rmodel = derivativesFEISTYR, 
#'                        bCust = TRUE)
#'
#' @param p A complete parameter list. \cr 
#' The parameter of FEISTY setups can be one of the following: \code{\link{setupBasic}}, \code{\link{setupBasic2}}, \code{\link{setupVertical}}, and \code{\link{setupVertical2}}.
#' Moreover, users can modify these four prepared FEISTY setups or customize new setups before calling \code{simulateFEISTY}.
#' @param tEnd The end time for the simulation [year], i.e., simulation period of FEISTY, in years. 
#' @param tStep The time step for ODE solving output [year]. Default is 1. Time series results are in every year.
#' @param times A sequence of time points for FEISTY simulations (ODEs solving), required by function \code{\link{ode}}. For general use, the parameter does not need to be defined, since it will be generated by `tEnd` and `tStep` automatically. \cr
#'              If `NA`, the function returns only the derivatives (one time step running) by the derivative function, which is the \code{\link{derivativesFEISTYR}} in R or the same function in FORTRAN.
#' @param yini A vector containing initial biomass values of all state variables (resources and all size classes). The default is imported from the setup parameter list, `p$u0`. 
#'             If input the whole output list of a simulation (e.g., \code{yini=sim}), the biomass values of all state variables of the last time step (\code{yini=sim$u[sim$nTime,]}) will be applied to the initial values.
#'             It can be used for successive simulations based on previous simulations. 
#' @param USEdll Logical flag, determining whether the ODEs are solved in FORTRAN (`TRUE`) or R (`FALSE`). \cr
#' The \link{deSolve} package is required for both methods. Default is TRUE, \code{USEdll=FALSE} is useful in debugging or model development.
#' `bCust` flag input is ineffective when \code{USEdll=FALSE}.
#' @param Rmodel The R function for computing derivatives, defaults to \code{\link{derivativesFEISTYR}}. Generally, it should not be changed, unless users modify the model profoundly.
#' @param bCust Logical flag, indicates whether to use fixed setups (FALSE) or customized setups (TRUE). \cr 
#' Default is TRUE, which means the core FEISTY parameters generated in R are transmitted to Fortran, and the ode solving is also done by compiled language. 
#' \code{bCust=FALSE} is useful in debugging and model development, e,g., comparing R and FORTRAN results.
#' \code{bCust} flag has a lower priority than \code{USEdll} flag. \code{bCust} flag input is ineffective when \code{USEdll} flag is FALSE.
#' 
#' @details
#' The function runs the FEISTY model simulation over the specified time frame. \cr
#' The simulation supports published FEISTY setups and their revised versions:
#' \code{\link{setupBasic}}, \code{\link{setupBasic2}}, \code{\link{setupVertical}}, and \code{\link{setupVertical2}}, and customized setups by modelers. \cr
#' The simulation can be conducted by a FORTRAN-based approach or an R-based approach. Both methods rely on the \link{deSolve} package for ODE solving. 
#' For efficiency, FORTRAN dll should be used. For model development, the R-version is preferred.
#' Simulations based on customized setups can only be done with the FORTRAN-based approach.
#'
#' @return
#' A list containing the simulation results:
#' \itemize{
#' \item u: a matrix of biomass of each state variable (column) at each time point (row), including resources and all size classes of functional types.
#' \item R: a matrix of biomass of each resource (column) at each time point (row).
#' \item B: a matrix of biomass of each size class (column) at each time point (row).
#' \item t: a vector containing all the simulation time points. From 0 to `tEnd`.
#' \item nTime: the number of time points.
#' \item USEdll: from parameter input.
#' \item p: the parameter list used in the simulation, the same as the input one.
#' \item f: a matrix containing feeding levels [-] of all size classes of functional types over each time point. Resources not included.
#' \item mortpred: a matrix containing a vector containing predation mortality rate [1/year] of all resources and all size classes of functional types over each time point. Resources not included.
#' \item g: a matrix containing the net growth rate [1/year] of all size classes of functional types over each time point. Resources not included.
#' \item Repro: a matrix containing the energy used for reproduction of all size classes of functional types over each time point, rate [g/m2/year]. Resources not included.
#' \item Fin: a matrix containing the biomass flux into each size class over each time point [g/m2/year]. Resources not included.
#' \item Fout: a matrix containing the biomass flux out of each size class over each time point [g/m2/year]. Resources not included.
#' \item totMort: a matrix containing the total mortality [g/m2/year] of each functional type over each time point, 
#' which includes predation mortality, background mortality, and fishing mortality.
#' \item totGrazing: a matrix containing the total grazing (food intake before assimilation) [g/m2/year] of each functional type over each time point.
#' Cmax * f * u (maximum consumption rate * feeding level * biomass).
#' \item totLoss: a matrix containing all biomass loss [g/m2/year] of each functional type over each time point, including unassimilated food, basal metabolism and reproduction cost (1-epsRepro). 
#' The reproduction cost here also include the energy loss (1-epsRepro) from the flux out of the last size class (invested in reproduction) of each functional type. 
#' These losses are supposed to be released to environments.
#' \item totRepro: a matrix containing the total energy used for reproduction [g/m2/year] of each functional type over each time point (before multiplying the reproduction efficiency epsRepro).
#' \item totRecruit: a matrix containing the total recruitment [g/m2/year] of each functional type over each time point. 
#' totRecruit = totRepro * epsRepro (reproduction efficiency)
#' \item totBiomass: a matrix containing the total biomass [g/m2] of each functional type over each time point.
#' \item `SSBAMean`, `SSBGMean`, `SSBMin`, `SSBMax`, and `SSB` can be found in \code{\link{calcSSB}}. \cr
#' `yieldAMean`, `yieldGMean`, `yieldMin`, `yieldMax`, and `yield` can be found in \code{\link{calcYield}}.
#' }
#' 
#' @examples
#' 
#' # run setupBasic with default parameter settings
#' sim <- simulateFEISTY()
#' plotSimulation(sim)
#' 
#' # run FEISTY simulation based on setupVertical2
#' # prepare a parameter list
#' p_V <- setupVertical2(szprod = 100, lzprod = 120, bprodin = 200, depth = 1000)
#' # run the simulation by R
#' sim_Vertical_R <- simulateFEISTY(p = p_V, USEdll = FALSE)
#' plotSimulation(sim_Vertical_R)
#' 
#' # run FEISTY simulation based on setupBasic2 by Fortran
#' sim_Basic2_F <- simulateFEISTY(p = setupBasic2(szprod = 90, lzprod = 100, bprod = 15, 
#'                                                depth = 500, Tp = 11, Tb = 9, 
#'                                                nStages=9, etaMature=0.25, F=0, etaF=0.05),
#'                                 tEnd = 1000, tStep = 1, USEdll = TRUE)
#' plotSimulation(sim_Basic2_F)
#' 
#' # -------------------------------------------------------------------------------
#' 
#' # run FEISTY simulation based on a customized set up. 
#' # The parameter values are provided solely for illustrative purposes.
#' 
#' # Initialize the parameter list.
#' p_cust <- paramInit(szprod=100, lzprod=100, bprod=50, Tp=10,Tb=10,etaMature=0.25, depth=500,
#' mMedium = 0.5, mLarge = 250)
#' 
#' # add three resources
#' p_cust <- paramAddResource(p_cust,
#'           names= c("smallZoo", "largeZoo", "benthos"),
#'           K    = c(100, 120, 80),
#'           r    = c(1, 1, 1),
#'           mLower = c(2e-06,0.001, 0.5e-03),
#'           mUpper = c(0.001, 0.5, 125),
#'           mc   = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000)))
#'           
#' # add two functional types of fish: small pelagic fish and demersal fish           
#' p_cust <- paramAddGroup(p_cust, mMin=0.001, mMax=250, mMature=NA, 
#'                         mortF=0,      nStages=6, name="smallPel")
#' p_cust <- paramAddGroup(p_cust, mMin=0.001, mMax=125000, mMature=NA, 
#'                         mortF=0, nStages=9, name="demersals")
#'                         
#' # add physiological parameters for two functional types
#' p_cust <- paramAddPhysiology(p_cust, 
#'           ac = 20, bc = -0.25,       
#'           am = 0.011*365, bm = -0.175,      
#'           ae = 70, be = -0.2,        
#'           epsRepro = 0.01, 
#'           epsAssim = 0.7)
#' 
#' # Add fishing mortality on demersal fish only. The baseline fishing mortality is 0.2/year.
#' p_cust <- setFishing(p_cust, Fmax=0.2, etaF=0.05,groupidx=c(2))
#' 
#' # Add size preference
#' p_cust$theta <- paramSizepref(p = p_cust,   
#'                         beta = 400,
#'                         sigma = 1.3, 
#'                         type = 1)
#' 
#' # Add temperature effect
#' p_cust=paramTeffect(p=p_cust, Tref=10, Q10=1.88, Q10m=2.35, pelgroupidx=c(1), demgroupidx=c(2))
#' # Turn on the effective temperature effects on large demersals.
#' p_cust$bET=TRUE
#' 
#' # run the simulation for 500 years. 
#' sim_cust <- simulateFEISTY(bCust = TRUE, p = p_cust, tEnd = 500)
#' 
#' @references
#' Petrik, C. M., Stock, C. A., Andersen, K. H., van Denderen, P. D., & Watson, J. R. (2019). Bottom-up drivers of global patterns of demersal, forage, and pelagic fishes. Progress in oceanography, 176, 102124.
#' 
#' van Denderen, P. D., Petrik, C. M., Stock, C. A., & Andersen, K. H. (2021). Emergent global biogeography of marine fish food webs. Global Ecology and Biogeography, 30(9), 1822-1834.
#' 
#' de Roos, A. M., Schellekens, T., Van Kooten, T., Van De Wolfshaar, K., Claessen, D., & Persson, L. (2008). Simplifying a physiologically structured population model to a stage-structured biomass model. Theoretical population biology, 73(1), 47-62.
#' 
#' Soetaert, K., Petzoldt, T., & Setzer, R. W. (2010). Solving differential equations in R: package deSolve. Journal of statistical software, 33, 1-25.
#' 
#' @author Ken H. Andersen, Karline Soetaert, Yixin Zhao
#'
#' @aliases simulateFEISTY
#' 
#' @seealso
#' 
#' \code{\link{setupBasic}} The setup following Petrik et al. (2019) \cr
#' \code{\link{setupBasic2}} A revised setup based on `setupBasic` \cr
#' \code{\link{setupVertical}} The setup following van Denderen et al. (2021) \cr
#' \code{\link{setupVertical2}} A revised setup based on `setupVertical` \cr
#' 
#' \code{\link{calcSSB}} Spawning stock biomass calculation \cr
#' \code{\link{calcYield}} Yield calculation
#' 
#' \code{\link{derivativesFEISTYR}} The derivative function of state variables in FEISTY model
#' 
#' \code{\link{webFEISTY}} A shiny interface for visualizing FEISTY model results
#' 
#' \code{\link{plotSimulation}} Plot simulation results including rates, biomass, and SSB data
#' 
#' @export
#'

simulateFEISTY = function(p      = setupBasic(), 
                          tEnd   = 500,
                          tStep  = 1,
                          times  = seq(from=0, to=tEnd, by=tStep),  
                          yini   = p$u0,  
                          USEdll = TRUE,
                          Rmodel = derivativesFEISTYR,
                          bCust  = TRUE)
{
  
  nR      <- p$nResources[1]  # no of resources. [1] to make sure that this is only one number
  nGroups <- p$nGroups[1] # no of fish groups
  nGrid   <- p$nStages[1] # no of grid points
  nFGrid  <- nGrid-nR # grid points of fish
  
  if (length(yini) != nGrid) 
    stop ("length of 'yini' not ok - should be ", nGrid)  
  
  # Set tolerances for ode solving
  rtol = 1E-8
  atol = 1E-8
  if (max(sapply(p$ix, length))>=21){    
    rtol = 1E-10
    atol = 1E-10}
  if (max(sapply(p$ix, length))>27) stop("The size number cannot be more than 27 due to the low accuracy of integration.")
  
  # prepare output variable names
  Sname <- p$stagenames
  Fname <- p$stagenames[-(1:nR)]
  Gname <- p$groupnames[-(1:nR)]
  outnames <- c(
    paste("f", Fname, sep="."), paste("mortpred", Sname, sep="."),
    paste("g", Fname, sep="."), paste("Repro", Fname, sep="."),
    paste("Fin", Fname, sep="."), paste("Fout", Fname, sep="."),
    paste("totMort", Gname, sep="."), paste("totGrazing", Gname, sep="."),
    paste("totLoss", Gname, sep="."), paste("totRepro", Gname, sep="."),
    paste("totRecruit", Gname, sep="."), paste("totBiomass", Gname, sep="."))    
  
  #
  # calculate in Fortran
  #
  
  if (USEdll==TRUE){

    # names of functions in fortran code to be used
    runfunc  <- "runfeisty"    # the derivative function
    
  if (bCust==TRUE) {    
    # the integers to be passed to the fortran code
    ipar <- c(nGroups,                           # total number of groups
              nR,                                # total number of resources
              unlist(lapply(p$ix, FUN=length)),  # number of stages per fish group
              p$Rtype,                           # type of resource dynamics
              if (is.null(p$pelgrididx)) 1 else length(p$pelgrididx), # length of pelgrididx. 1, if not defined.
              if (is.null(p$pelgrididx)) 1 else p$pelgrididx, # all pelagic fish grid indices. 1, if not defined.
              if (is.null(p$allgrididx)) 1 else length(p$allgrididx), # length of allgrididx. 1, if not defined.
              if (is.null(p$allgrididx)) 1 else p$allgrididx, # all grid indices (resources+fish). 1, if not defined.
              if (is.null(p$lgdemidx))   1 else length(p$lgdemidx), # length of lgdemidx. 1, if not defined.
              if (is.null(p$lgdemidx))   1 else p$lgdemidx, # large demersal fish indices. 1, if not defined.
              if (is.null(p$bET))        0 else as.integer(p$bET)) # effective temperature Boolean flag. 0 (FALSE), if not defined.
    ipar <- as.integer(ipar)
    if (length(c(nGroups,nR,unlist(lapply(p$ix, FUN=length)),p$Rtype)) != 3 + nGroups)
      stop ("length of 'ipar' not ok; check parameters")
    
    if (any(dim(p$theta)-c(nGrid, nGrid) != 0))
      stop ("dimension of 'theta' not ok: should be (", nGrid, ",", nGrid, ")")  
    
    # the double precision numbers to be passed to the fortran code
    rpar   <- c(rep(p$K,  length.out=nR),            # resource parameters
                rep(p$r,  length.out=nR),  
                rep(p$epsRepro, length.out=nGroups), # group-specific parameter
                p$psiMature[-(1:nR)],                # fish-stage parameter
                p$z[-(1:nR)], 
                t(p$theta),                          # check if not transpose
                rep(p$epsAssim,   length.out=nGrid), # all  
                rep(p$V,          length.out=nGrid), 
                rep(p$Cmax,       length.out=nGrid),
                rep(p$metabolism, length.out=nGrid),
                rep(p$mort0,      length.out=nGrid),
                rep(p$mortF,      length.out=nGrid),
                rep(p$Vsave,          length.out=nGrid), 
                rep(p$Cmaxsave,       length.out=nGrid),
                rep(p$metabolismsave, length.out=nGrid),
                p$depth,
                p$Q10,
                p$Q10m,
                p$Tp,
                p$Tb)
    
    # 
    # Run the simulation:
    #
    
      initfunc <- "initfeisty"
      
      if (any(is.na(times)))  # one call and return
        return( DLLfunc(y=yini, times=0, parms=NULL, dllname = "FEISTY",
                        func=runfunc, initfunc=initfunc, outnames=outnames, nout=length(outnames),
                        ipar=ipar, rpar=as.double(rpar)))
      
      u = ode(y=yini, times=times, parms=NULL, dllname = "FEISTY",
              func=runfunc, initfunc=initfunc, outnames=outnames, nout=length(outnames),
              ipar=ipar, rpar=as.double(rpar),
              method = "ode45", rtol = rtol, atol = atol) # Run by dll
    }
    else
    {     # for fixed setups
      # Transmit input file path to Fortran library
      passpath <- function() {
        sys=Sys.info()['sysname']
        
        if (sys=='Darwin') {
          sLibname = system.file("libs", "FEISTY.so", package = "FEISTY")
        }
        if (sys=='Linux') {
          sLibname = system.file("libs", "FEISTY.so", package = "FEISTY")
        }
        if (sys=='Windows'){
          if (Sys.info()['machine']=='x86-64'){
            sLibname = system.file("libs/x64", "FEISTY.dll", package = "FEISTY")
          }else{
            sLibname = system.file("libs/i386", "FEISTY.dll", package = "FEISTY")
          }
        }
        
        # Reload dll to avoid crash?
        if (is.loaded("runfeisty")) { # "runfeisty" is a function name in R_init_feisty.c
          if (sLibname == "") {
            print("sLibname is an empty string")
          } else{
            dyn.unload(sLibname)
            dyn.load(sLibname)
          }
        }
        
        file_path=system.file("extdata", "input.nml", package = "FEISTY")
        dummy=.C("passpath", length=nchar(file_path), file_path_in = charToRaw(file_path))
        file_path_V=system.file("extdata", "tempdata.dat", package = "FEISTY")
        dummy=.C("passpathv", length=nchar(file_path_V), file_path_in=charToRaw(file_path_V))
      }
      
      # Call the Fortran subroutine to pass input file path
      passresult <- passpath()
      
      # Choose the setup:
      if (p$setup=="setupBasic"){
        initfunc <- "initfeistysetupbasic"
        setupinput=c(p$szprod,p$lzprod,p$bprodin,p$dfbot,p$depth,p$Tp,p$Tb)
      }else if(p$setup=="setupBasic2"){
        initfunc <- "initfeistysetupbasic2"
        setupinput=c(p$szprod,p$lzprod,p$bprodin,p$dfbot,length(p$ix[[p$nGroups]]),p$depth,p$Tp,p$Tb,p$etaMature,p$Fmax,p$etaF,as.integer(p$bET))
      }else if(p$setup=="setupVertical"){
        initfunc <- "initfeistysetupvertical"
        setupinput = c(p$szprod,p$lzprod,p$bprodin,p$dfbot,p$dfpho,p$region, p$bottom, p$photic)
      }else if(p$setup=="setupVertical2"){
        initfunc <- "initfeistysetupvertical2"
        setupinput = c(p$szprod,p$lzprod,p$bprodin,p$dfbot,p$dfpho,length(p$ix[[p$nGroups]]), p$Tp, p$Tm, p$Tb, p$bottom,p$photic,p$etaMature,
                       p$shelfdepth,p$visual,p$Fmax,p$etaF)
      }
      
      if (any(is.na(times)))  # one call and return
        return( DLLfunc(y=yini, times=0, parms=as.double(setupinput), dllname = "FEISTY",
                        func=runfunc, initfunc=initfunc, outnames=outnames, nout=length(outnames),
                        ipar=NULL, rpar=NULL))
      # Full simulation:
      u = ode(y=yini, times=times, parms=as.double(setupinput), dllname = "FEISTY",
              func=runfunc, initfunc=initfunc, outnames=outnames, nout=length(outnames),
              ipar=NULL, rpar=NULL,
              method = "ode45", rtol = rtol, atol = atol) # Run by dll
    }   
    #
    # Calculate in R:
    #
  } else if (any(is.na(times))) {  # one call and return
    return (Rmodel(0, yini, p))
  } else {               # R-code
    u = ode(y=yini, times=times, parms=p, func = Rmodel,
            method = "ode45", rtol = rtol, atol = atol) #Run by R
    # assign colnames
    colnames(u)[(1+p$nStages+1):ncol(u)]=outnames
  }
  
  #
  # Assemble output:
  #
  
  sim   = list()
  sim$u = u[,c(p$ixR,p$ixFish)+1]
  sim$R = u[, p$ixR+1]
  sim$B = u[, p$ixFish+1]
  sim$t = times
  sim$nTime = length(times)
  sim$USEdll=USEdll
  sim$p = p
  
  # feeding level
  # "^xx" extracting data starts with "xx"
  col_f=grep("^f.", colnames(u), value = TRUE)
  sim$f=u[,col_f]  
  # predation mortality rate
  col_mortpred=grep("^mortpred", colnames(u), value = TRUE)
  sim$mortpred=u[,col_mortpred]
  # net growth rate
  col_g=grep("^g.", colnames(u), value = TRUE)
  sim$g=u[,col_g]
  # Energy used for reproduction [g/m2/year]
  col_Repro=grep("^Repro", colnames(u), value = TRUE)
  sim$Repro=u[,col_Repro]
  # Biomass flux into each size class
  col_Fin=grep("^Fin", colnames(u), value = TRUE)
  sim$Fin=u[,col_Fin]
  # Biomass flux out of each size class
  col_Fout=grep("^Fout", colnames(u), value = TRUE)
  sim$Fout=u[,col_Fout]
  # total mortality of each functional group [g/m2/year],which includes predation mortality, background mortality, and fishing mortality.
  col_totMort=grep("^totMort", colnames(u), value = TRUE)
  sim$totMort=u[,col_totMort]
  # total grazing of each functional group [g/m2/year], Cmax * f (maximum consumption rate * feeding level), the food intake before assimilation.
  col_totGrazing=grep("^totGrazing", colnames(u), value = TRUE)
  sim$totGrazing=u[,col_totGrazing]
  # total biomass loss of each functional group [g/m2/year], including unassimilated food and metabolism. They are released to environments.
  col_totLoss=grep("^totLoss", colnames(u), value = TRUE)
  sim$totLoss=u[,col_totLoss]
  # total energy used for reproduction of each functional group [g/m2]
  col_totRepro=grep("^totRepro", colnames(u), value = TRUE)
  sim$totRepro=u[,col_totRepro]
  # total recruitment of each functional group [g/m2], TotRecruit = TotRepro * epsRepro (reproduction efficiency)
  col_totRecruit=grep("^totRecruit", colnames(u), value = TRUE)
  sim$totRecruit=u[,col_totRecruit] 
  # total biomass of each functional group [g/m2]
  col_totBiomass=grep("^totBiomass", colnames(u), value = TRUE)
  sim$totBiomass=u[,col_totBiomass]  
  
  #
  # Calculate Spawning Stock Biomass and yield
  #
  sim=calcSSB(sim=sim,etaTime=0.4)
  sim=calcYield(sim=sim,etaTime=0.4)
  
  return(structure(sim, class = 'FEISTY'))
}


derivativesFEISTYR_ts = function(t,              # current time
                                 u,              # all state variables
                                 p,              # parameters
                                 FullOutput=TRUE) {
  p$Tp = p$getts(time=t,y=p$Tp_ts)
  p$Tm = p$getts(time=t,y=p$Tm_ts)
  p$Tb = p$getts(time=t,y=p$Tb_ts)
  p$r[3] = p$getts(time=t,y=p$bprod_ts)
  #p = paramAddPhysiology(p)
  if(p$setup == "setupVertical2"){
    p = paramTeffect_vet(p)
  }else if(p$setup == "setupBasic" | p$setup == "setupBasic2"){
    p = paramTeffect(p=p, # only for setupbasic & 2
                     Tref=p$Tref,
                     Q10=p$Q10,
                     Q10m=p$Q10m,
                     pelgroupidx=c(1:(p$nGroups-1)),
                     demgroupidx=p$nGroups)  
  }
  
  u[1]=p$getts(time=t,y=p$szbio_ts)
  u[2]=p$getts(time=t,y=p$lzbio_ts)
  
  #print(t)
  # split state variable vector into resource and fish
  u[u<0]=0
  R     = u[p$ixR]       # resource, prey
  iFish = p$ixFish
  B     = u[iFish]       # fish
  
  # update effective temperature for large demersal fish in shallow water
  if(!is.null(p$depth) & !is.null(p$bET))
    if (p$depth<200 & p$bET==TRUE) p=updateET(p=p,u=u)
  
  # ----------------------------------------------
  # Consumption of all fish groups
  # ----------------------------------------------
  
  # V: clearance rate, (m2/g/yr) 
  # theta x u: prey available for consumption
  # Encounter # /yr
  Enc = p$V * (p$theta %*% u) 
  
  # ----------------------------------------------
  # Predation mortality, /yr:
  # = t(p$theta) %*% (f*p$Cmax/p$epsAssim*u/p$mc)
  # ----------------------------------------------
  #
  mm = p$Cmax*p$V/(Enc+p$Cmax)*u # temporarily store
  mm[ is.na(mm) ] = 0
  mortpred = t(p$theta) %*% mm
  
  szprod = p$getts(time=t,y=p$szprod_ts)
  lzprod = p$getts(time=t,y=p$lzprod_ts)
  
  dr_fac_theta = matrix(1, nrow = nrow(p$theta), ncol = ncol(p$theta)) 
  if(mortpred[1]*u[1] > szprod | mortpred[2]*u[2] > lzprod){
    
    if (mortpred[1]*u[1] >= szprod) {
      dr_fac_sz = szprod/(mortpred[1]*u[1])
      dr_fac_theta[p$ixFish,1] = dr_fac_sz
      mortpred[1]=dr_fac_sz*mortpred[1]
    }
    if (mortpred[2]*u[2] >= lzprod) {
      dr_fac_lz = lzprod/(mortpred[2]*u[2])
      dr_fac_theta[p$ixFish,2] = dr_fac_lz
      mortpred[2]=dr_fac_lz*mortpred[2]
    }
  }
  
  # f: feeding level
  # Cmax: maximum consumption rate, /yr
  # f = Enc / (p$Cmax + Enc)
  f   = (p$V * ((p$theta*dr_fac_theta) %*% u)) / (p$Cmax + Enc)   # Functional response
  f[is.na(f)] = 0
  
  # net growth rate, /yr
  Eavail  = p$epsAssim * p$Cmax * f - p$metabolism
  
  # ----------------------------------------------
  # Total mortality (includes basal and fishing mortality)
  # ----------------------------------------------
  mort = mortpred + p$mort0 + p$mortF   # /year
  
  # ----------------------------------------------
  # Derivative of fish groups
  # ----------------------------------------------
  
  # Flux out of the size group
  #------------------------------
  v     = Eavail[iFish,]   # net growth rate
  vplus = pmax(v,0)
  
  # fraction available for growth
  kappa = 1 - p$psiMature[iFish]   
  g     = kappa*vplus
  
  # growth to the next stage
  gamma = (kappa*vplus - mort[iFish]) /
    (1 - (1/p$z[iFish])^(1-mort[iFish]/(kappa*vplus)) )
  
  gamma[kappa==0] = 0 # No growth of fully mature classes
  
  # goes out of stage (size group)
  Fout = gamma*B
  
  # Energy used for reproduction
  Repro = p$psiMature[iFish]*vplus*B
  
  # Flux into the size group
  #------------------------------
  Fin = Fout*0
  for (i in 1:p$nGroups) {
    ix = p$ix[[i]] - p$ix[[1]][1] +1                  # growth to
    ixPrev  = c(ix[length(ix)], ix[1:(length(ix)-1)]) # growth from
    Fin[ix] = Fout[ixPrev]
    
    # for reproduction: consider the reproduction success
    Fin[ix[1]] = p$epsRepro[i]*(Fin[ix[1]] + sum( Repro[ix] ))
  }
  
  # ----------------------------------------------
  # Assemble derivatives of fish:
  # ----------------------------------------------
  
  dBdt = Fin - Fout + (v - mort[p$ixFish])*B - Repro
  
  # ----------------------------------------------
  # Derivative of resources
  # ----------------------------------------------
  
  # if (p$Rtype == 1)  # chemostat
  #   dRdt = p$r*(p$K-R) - mortpred[p$ixR]*R
  # else               # logistic
  #   dRdt = p$r*R*(1-R/p$K) - mortpred[p$ixR]*R
  dRdt = c(0,0,0,0)
  dRdt[3] = p$r[3]*R[3]*(1-R[3]/p$K[3]) - mortpred[3]*R[3]
  
  # ----------------------------------------------
  # Assemble output:
  # ----------------------------------------------
  if (FullOutput) { # Output everything
    out = list()
    out$deriv = c(dRdt, dBdt)
    out$f     = f[-p$ixR,] # Feeding level only all fish stages, no resources
    out$mortpred = mortpred[,]
    out$g     = g # net growth rate fish stages
    out$Repro = Repro
    out$Fin   = Fin
    out$Fout  = Fout
    
    # for the budget:
    grazing = p$Cmax * f         # grazing rate, /yr
    loss    = (1.-p$epsAssim) * grazing + p$metabolism # Energy loss to environments. Updated below.
    Reprofrac = (1-kappa)*vplus
    
    il <- NULL
    for (i in 1:length(p$ix)){
      # Add the waste energy in reproduction of each stages. Note it does not include the waste energy from last stage energy flux out, added in `totLoss` below.
      loss[p$ix[[i]]] = loss[p$ix[[i]]] + (1-p$epsRepro[i]) * Reprofrac[p$ix[[i]]-p$nResources] 
      il <- c(il, rep(i, times=length(p$ix[[i]])))
    }
    
    out$totMort    = tapply((mort   *u)[p$ixFish], INDEX=il, FUN=sum)
    out$totGrazing = tapply((grazing*u)[p$ixFish], INDEX=il, FUN=sum)
    # Add the waste energy in reproduction from flux out of the last stage of each functional type.
    out$totLoss    = tapply((loss   *u)[p$ixFish], INDEX=il, FUN=sum) +(1-p$epsRepro)*Fout[sapply(p$ix, tail, n = 1)-p$nResources] 
    out$totRepro   = tapply(Repro, INDEX=il, FUN=sum) + Fout[sapply(p$ix, tail, n = 1)-p$nResources] 
    out$totRecruit = out$totRepro* p$epsRepro
    out$totBiomass = tapply(B, INDEX=il, FUN=sum)
    return(out)
  }
  else # Output just the derivatives
    return( list(c(dRdt, dBdt)) )
}


# only setupBasic setupBasic2 setupVertical2
simulateFEISTY_ts = function(p      = setupBasic(), 
                             tEnd   = 1,
                             tStep  = 1/12,
                             times  = seq(from=0, to=tEnd, by=tStep),  
                             yini   = p$u0,  
                             USEdll = TRUE,
                             Rmodel = derivativesFEISTYR_ts,
                             bCust  = TRUE,
                             spinup = T){
  
  nR      <- p$nResources[1]  # no of resources. [1] to make sure that this is only one number
  nGroups <- p$nGroups[1] # no of fish groups
  nGrid   <- p$nStages[1] # no of grid points
  nFGrid  <- nGrid-nR # grid points of fish
  
  if (length(yini) != nGrid) 
    stop ("length of 'yini' not ok - should be ", nGrid)  
  
  # Set tolerances for ode solving
  rtol = 1E-8
  atol = 1E-8
  if (max(sapply(p$ix, length))>=21){    
    rtol = 1E-10
    atol = 1E-10}
  if (max(sapply(p$ix, length))>27) stop("The size number cannot be more than 27 due to the low accuracy of integration.")
  
  # prepare output variable names
  Sname <- p$stagenames
  Fname <- p$stagenames[-(1:nR)]
  Gname <- p$groupnames[-(1:nR)]
  outnames <- c(
    paste("f", Fname, sep="."), paste("mortpred", Sname, sep="."),
    paste("g", Fname, sep="."), paste("Repro", Fname, sep="."),
    paste("Fin", Fname, sep="."), paste("Fout", Fname, sep="."),
    paste("totMort", Gname, sep="."), paste("totGrazing", Gname, sep="."),
    paste("totLoss", Gname, sep="."), paste("totRepro", Gname, sep="."),
    paste("totRecruit", Gname, sep="."), paste("totBiomass", Gname, sep="."))    
  
  if (spinup == T) {
    loopnum=4
    timesspin=seq(from=0, to=0.1*tEnd, by=tStep)
    pspin = p
    
    pspin$szbio_ts = p$szbio_ts[1:(0.1*length(p$szbio_ts))]
    pspin$szbio_ts[length(pspin$szbio_ts) + 1] = pspin$szbio_ts[length(pspin$szbio_ts)]
    pspin$lzbio_ts = p$lzbio_ts[1:(0.1*length(p$lzbio_ts))]
    pspin$lzbio_ts[length(pspin$lzbio_ts) + 1] = pspin$lzbio_ts[length(pspin$lzbio_ts)]
    pspin$bprod_ts = p$bprod_ts[1:(0.1*length(p$bprod_ts))]
    pspin$bprod_ts[length(pspin$bprod_ts) + 1] = pspin$bprod_ts[length(pspin$bprod_ts)]
    pspin$szprod_ts = p$szprod_ts[1:(0.1*length(p$szprod_ts))]
    pspin$szprod_ts[length(pspin$szprod_ts) + 1] = pspin$szprod_ts[length(pspin$szprod_ts)]
    pspin$lzprod_ts = p$lzprod_ts[1:(0.1*length(p$lzprod_ts))]
    pspin$lzprod_ts[length(pspin$lzprod_ts) + 1] = pspin$lzprod_ts[length(pspin$lzprod_ts)]
    pspin$Tp_ts = p$Tp_ts[1:(0.1*length(p$Tp_ts))]
    pspin$Tp_ts[length(pspin$Tp_ts) + 1] = pspin$Tp_ts[length(pspin$Tp_ts)]
    pspin$Tm_ts = p$Tm_ts[1:(0.1*length(p$Tm_ts))]
    pspin$Tm_ts[length(pspin$Tm_ts) + 1] = pspin$Tm_ts[length(pspin$Tm_ts)]
    pspin$Tb_ts = p$Tb_ts[1:(0.1*length(p$Tb_ts))]
    pspin$Tb_ts[length(pspin$Tb_ts) + 1] = pspin$Tb_ts[length(pspin$Tb_ts)]
  }
  
  #
  # calculate in Fortran
  #
  
  if (USEdll==TRUE){
    
    # names of functions in fortran code to be used
    runfunc  <- "runfeisty_ts"    # the derivative function
    
    if (bCust==TRUE) {    
      # the integers to be passed to the fortran code
      ipar <- c(nGroups,                           # total number of groups
                nR,                                # total number of resources
                unlist(lapply(p$ix, FUN=length)),  # number of stages per fish group
                p$Rtype,                           # type of resource dynamics
                if (is.null(p$pelgrididx)) 1 else length(p$pelgrididx), # length of pelgrididx. 1, if not defined.
                if (is.null(p$pelgrididx)) 1 else p$pelgrididx, # all pelagic fish grid indices. 1, if not defined.
                if (is.null(p$allgrididx)) 1 else length(p$allgrididx), # length of allgrididx. 1, if not defined.
                if (is.null(p$allgrididx)) 1 else p$allgrididx, # all grid indices (resources+fish). 1, if not defined.
                if (is.null(p$lgdemidx))   1 else length(p$lgdemidx), # length of lgdemidx. 1, if not defined.
                if (is.null(p$lgdemidx))   1 else p$lgdemidx, # large demersal fish indices. 1, if not defined.
                if (is.null(p$bET))        0 else as.integer(p$bET)) # effective temperature Boolean flag. 0 (FALSE), if not defined.
      ipar <- as.integer(ipar)
      if (length(c(nGroups,nR,unlist(lapply(p$ix, FUN=length)),p$Rtype)) != 3 + nGroups)
        stop ("length of 'ipar' not ok; check parameters")
      
      if (any(dim(p$theta)-c(nGrid, nGrid) != 0))
        stop ("dimension of 'theta' not ok: should be (", nGrid, ",", nGrid, ")")  
      
      # the double precision numbers to be passed to the fortran code
      rpar   <- c(rep(p$K,  length.out=nR),            # resource parameters
                  rep(p$r,  length.out=nR),  
                  rep(p$epsRepro, length.out=nGroups), # group-specific parameter
                  p$psiMature[-(1:nR)],                # fish-stage parameter
                  p$z[-(1:nR)], 
                  t(p$theta),                          # check if not transpose
                  rep(p$epsAssim,   length.out=nGrid), # all  
                  rep(p$V,          length.out=nGrid), 
                  rep(p$Cmax,       length.out=nGrid),
                  rep(p$metabolism, length.out=nGrid),
                  rep(p$mort0,      length.out=nGrid),
                  rep(p$mortF,      length.out=nGrid),
                  rep(p$Vsave,          length.out=nGrid), 
                  rep(p$Cmaxsave,       length.out=nGrid),
                  rep(p$metabolismsave, length.out=nGrid),
                  p$depth,
                  p$Q10,
                  p$Q10m,
                  p$Tp,
                  p$Tb)
      
      # 
      # Run the simulation:
      #
      
      initfunc <- "initfeisty"
      
      if (any(is.na(times)))  # one call and return
        return( DLLfunc(y=yini, times=0, parms=NULL, dllname = "FEISTY",
                        func=runfunc, initfunc=initfunc, outnames=outnames, nout=length(outnames),
                        ipar=ipar, rpar=as.double(rpar)))
      
      dummy=.Fortran("passnforc", nforcsin = as.integer(nFGrid*3+5) )
      
       #dummy=.C("passnforc", nforcsin = as.integer(nFGrid*3+5))
      
      if(spinup == T){
        pspin=buildforcings(timesspin,p=pspin)
        for (i in 1:loopnum) {
          u <- ode(y		= yini,
                   times		= timesspin,
                   parms		= NULL,
                   ipar = ipar, rpar = as.double(rpar),
                   dllname		= "FEISTY",
                   initfunc	= initfunc,
                   func		= runfunc,
                   initforc	= "initfeistyforc",
                   forcings	= pspin$forcings,
                   fcontrol	= list(method="constant", rule = 2, f = 0, ties = "ordered"),
                   method = "ode45", rtol = rtol, atol = atol,
                   outnames = outnames, nout = length(outnames))
          yini = u[length(timesspin),c(p$ixR,p$ixFish)+1]
          cat(sprintf("spin-up progress: %.2f%%\n", 100*i/loopnum))
        }
        p$u0 = u[length(timesspin),c(p$ixR,p$ixFish)+1]
        yini = p$u0
      }
      
      p=buildforcings(times,p)
      
      runfunc="runfeisty_ts"

      u <- ode(y		= yini,
               times		= times,
               parms		= NULL,
               ipar = ipar, rpar = as.double(rpar),
               dllname		= "FEISTY",
               initfunc	= initfunc,
               func		= runfunc,
               initforc	= "initfeistyforc",
               forcings	= p$forcings,
               fcontrol	= list(method="constant", rule = 2, f = 0, ties = "ordered"),
               method = "ode45", rtol = rtol, atol = atol,
               outnames = outnames, nout = length(outnames))
      
    }
    else
    {     # for fixed setups
      # Transmit input file path to Fortran library
      passpath <- function() {
        sys=Sys.info()['sysname']
        
        if (sys=='Darwin') {
          sLibname = system.file("libs", "FEISTY.so", package = "FEISTY")
        }
        if (sys=='Linux') {
          sLibname = system.file("libs", "FEISTY.so", package = "FEISTY")
        }
        if (sys=='Windows'){
          if (Sys.info()['machine']=='x86-64'){
            sLibname = system.file("libs/x64", "FEISTY.dll", package = "FEISTY")
          }else{
            sLibname = system.file("libs/i386", "FEISTY.dll", package = "FEISTY")
          }
        }
        
        # Reload dll to avoid crash?
        if (is.loaded("runfeisty")) { # "runfeisty" is a function name in R_init_feisty.c
          if (sLibname == "") {
            print("sLibname is an empty string")
          } else{
            dyn.unload(sLibname)
            dyn.load(sLibname)
          }
        }
        
        file_path=system.file("extdata", "input.nml", package = "FEISTY")
        dummy=.C("passpath", length=nchar(file_path), file_path_in = charToRaw(file_path))
        file_path_V=system.file("extdata", "tempdata.dat", package = "FEISTY")
        dummy=.C("passpathv", length=nchar(file_path_V), file_path_in=charToRaw(file_path_V))
      }
      
      # Call the Fortran subroutine to pass input file path
      passresult <- passpath()
      
      # Choose the setup:
      if (p$setup=="setupBasic"){
        initfunc <- "initfeistysetupbasic"
        setupinput=c(p$szprod,p$lzprod,p$bprodin,p$dfbot,p$depth,p$Tp,p$Tb)
      }else if(p$setup=="setupBasic2"){
        initfunc <- "initfeistysetupbasic2"
        setupinput=c(p$szprod,p$lzprod,p$bprodin,p$dfbot,length(p$ix[[p$nGroups]]),p$depth,p$Tp,p$Tb,p$etaMature,p$Fmax,p$etaF,as.integer(p$bET))
      }else if(p$setup=="setupVertical"){
        initfunc <- "initfeistysetupvertical"
        setupinput = c(p$szprod,p$lzprod,p$bprodin,p$dfbot,p$dfpho,p$region, p$bottom, p$photic)
      }else if(p$setup=="setupVertical2"){
        initfunc <- "initfeistysetupvertical2"
        setupinput = c(p$szprod,p$lzprod,p$bprodin,p$dfbot,p$dfpho,length(p$ix[[p$nGroups]]), p$Tp, p$Tm, p$Tb, p$bottom,p$photic,p$etaMature,
                       p$shelfdepth,p$visual,p$Fmax,p$etaF)
      }
      
      if (any(is.na(times)))  # one call and return
        return( DLLfunc(y=yini, times=0, parms=as.double(setupinput), dllname = "FEISTY",
                        func=runfunc, initfunc=initfunc, outnames=outnames, nout=length(outnames),
                        ipar=NULL, rpar=NULL))
      # Full simulation:
      u = ode(y=yini, times=times, parms=as.double(setupinput), dllname = "FEISTY",
              func=runfunc, initfunc=initfunc, outnames=outnames, nout=length(outnames),
              ipar=NULL, rpar=NULL,
              method = "ode45", rtol = rtol, atol = atol) # Run by dll
    }   
    #
    # Calculate in R:
    #
  } else if (any(is.na(times))) {  # one call and return
    return (Rmodel(0, yini, p))
  } else {               # R-code
    
    if (spinup == T) {
      pspin$getts=getts <- function(time, y) {
        approxfun(x = timesspin, y = y, method = "constant", rule = 2, f = 0, ties = "ordered")(time)
      }
      for (i in 1:loopnum) {
        u = ode(y=yini, times=timesspin, parms=pspin, func = Rmodel,
                method = "ode45", rtol = rtol, atol = atol) #Run by R
        yini = u[length(timesspin),c(p$ixR,p$ixFish)+1]
        cat(sprintf("spin-up progress: %.2f%%\n", 100*i/loopnum))
      }
      p$u0 = u[length(timesspin),c(p$ixR,p$ixFish)+1]
      yini = p$u0
    }
    p$getts=getts <- function(time, y) {
      approxfun(x = times, y = y, method = "constant", rule = 2, f = 0, ties = "ordered")(time)
    }
    
    
    u = ode(y=yini, times=times, parms=p, func = Rmodel,
            method = "ode45", rtol = rtol, atol = atol) #Run by R
    # assign colnames
    colnames(u)[(1+p$nStages+1):ncol(u)]=outnames
  }
  
  #
  # Assemble output:
  #
  u[,(1:2)+1]= cbind(p$szbio_ts,p$lzbio_ts)
  sim   = list()
  sim$u = u[,c(p$ixR,p$ixFish)+1]
  sim$R = u[, p$ixR+1]
  sim$B = u[, p$ixFish+1]
  sim$t = times
  sim$nTime = length(times)
  sim$USEdll=USEdll
  sim$p = p
  
  # feeding level
  # "^xx" extracting data starts with "xx"
  col_f=grep("^f.", colnames(u), value = TRUE)
  sim$f=u[,col_f]  
  # predation mortality rate
  col_mortpred=grep("^mortpred", colnames(u), value = TRUE)
  sim$mortpred=u[,col_mortpred]
  # net growth rate
  col_g=grep("^g.", colnames(u), value = TRUE)
  sim$g=u[,col_g]
  # Energy used for reproduction [g/m2/year]
  col_Repro=grep("^Repro", colnames(u), value = TRUE)
  sim$Repro=u[,col_Repro]
  # Biomass flux into each size class
  col_Fin=grep("^Fin", colnames(u), value = TRUE)
  sim$Fin=u[,col_Fin]
  # Biomass flux out of each size class
  col_Fout=grep("^Fout", colnames(u), value = TRUE)
  sim$Fout=u[,col_Fout]
  # total mortality of each functional group [g/m2/year],which includes predation mortality, background mortality, and fishing mortality.
  col_totMort=grep("^totMort", colnames(u), value = TRUE)
  sim$totMort=u[,col_totMort]
  # total grazing of each functional group [g/m2/year], Cmax * f (maximum consumption rate * feeding level), the food intake before assimilation.
  col_totGrazing=grep("^totGrazing", colnames(u), value = TRUE)
  sim$totGrazing=u[,col_totGrazing]
  # total biomass loss of each functional group [g/m2/year], including unassimilated food and metabolism. They are released to environments.
  col_totLoss=grep("^totLoss", colnames(u), value = TRUE)
  sim$totLoss=u[,col_totLoss]
  # total energy used for reproduction of each functional group [g/m2]
  col_totRepro=grep("^totRepro", colnames(u), value = TRUE)
  sim$totRepro=u[,col_totRepro]
  # total recruitment of each functional group [g/m2], TotRecruit = TotRepro * epsRepro (reproduction efficiency)
  col_totRecruit=grep("^totRecruit", colnames(u), value = TRUE)
  sim$totRecruit=u[,col_totRecruit] 
  # total biomass of each functional group [g/m2]
  col_totBiomass=grep("^totBiomass", colnames(u), value = TRUE)
  sim$totBiomass=u[,col_totBiomass]  
  
  #
  # Calculate Spawning Stock Biomass and yield
  #
  sim=calcSSB(sim=sim,etaTime=0.4)
  sim=calcYield(sim=sim,etaTime=0.4)
  
  return(structure(sim, class = 'FEISTY'))
}