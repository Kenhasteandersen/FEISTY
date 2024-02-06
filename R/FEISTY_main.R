#===============================================================================
# The FEISTY model as described in 
# Van Denderen et al. 2020 Global Ecology and Biogeography, DOI: 10.1111/geb.13348 
# Emergent global biogeography of marine fish food webs.

# Petrik, CM, Stock, CA, Andersen, KH, van Denderen, PD, Watson, JR 2019. 
# Bottom-up drivers of global patterns of demersal, forage, and pelagic fishes. 
# Progress in Oceanography, 176, 102124. https://doi.org/10.1016/j.pocean.2019.102124

# Slightly rewritten by Karline Soetaert, based on code from Ken H. Andersen
# 
# Main model routines
#      derivativesFEISTYR  : derivative code in R
#      derivativesFEISTYF  : derivative code in Fortran
#      simulateFEISTY      : run the feisty model
#
#===============================================================================

#' Calculate Derivatives of State Variables in FEISTY Model in R
#'
#' @description
#' This function calculates the derivatives of state variables in FEISTY model in R. 
#' It is used for time integration of FEISTY simulations or getting rates by computing derivatives of one time step.
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
#' ODE solving is conducted by package deSolve.
#' More theories can be found in de Roos et al. (2008) and Petrik et al. (2019).
#'
#' @return
#' If \code{FullOutput} is 'TRUE', a list containing the following components:
#' \itemize{
#' \item deriv: a vector of derivatives [gWW/m2/year] of all resources and all functional type size classes.
#' \item f: a vector containing feeding levels [-] of all resources (0) and all size classes of functional types.
#' \item mortpred: a vector containing predation mortality rate [1/year] of all resources and all size classes of functional types.
#' \item g: Net growth rate [1/year] of all size classes of functional types. Resources not included.
#' \item Repro: Energy used for reproduction of all size classes of functional types, rate [gWW/m2/year]. Resources not included.
#' \item Fin: Biomass flux into each size class [gWW/m2/year]. Resources not included.
#' \item Fout: Biomass flux out of each size class [gWW/m2/year]. Resources not included.
#' \item totMort: a vector containing total mortality [gWW/m2/year] of each functional type, 
#' which includes predation mortality, background mortality, and fishing mortality.
#' \item totGrazing: a vector containing total grazing [gWW/m2/year] of each functional type, Cmax * f (maximum consumption rate * feeding level)
#' To be simply, the food intake before assimilation.
#' \item totLoss: a vector containing all biomass loss [gWW/m2/year] of each functional type, including unassimilated food and metabolism.
#' They are released to environments. where is energy loss from reproduction (1-epsRepro), to be fixed.
#' \item totRepro: a vector containing total energy used for reproduction [gWW/m2] of each functional type.
#' \item totRecruit: a vector containing total recruitment [gWW/m2] of each functional type. 
#' TotRecruit = TotRepro * epsRepro (reproduction efficiency)
#' \item totBiomass: a vector containing total biomass [gWW/m2] of each functional type.
#'   }
#' If \code{FullOutput = FALSE}, only returns \code{deriv}.
#'
#' @examples
#' p = setpBasic2()
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
#' @author Ken H. Andersen, Karline Soetaert <karline.soetaert@nioz.nl>, Yixin Zhao
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
  if(!is.null(p$depth) & !is.null(p$bET) & p$depth<200 & p$bET==TRUE) p=updateET(p=p,u=u)
  
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
  v     = Eavail[iFish]   # net growth rate
  vplus = pmax(0, v)
  
  # fraction available for growth
  kappa = 1 - p$psiMature[iFish]   
  g     = kappa*vplus
  
  # growth to the next stage
  gamma = (kappa*vplus - mort[iFish]) /
    (1 - (1/p$z[iFish])^(1-mort[iFish]/(kappa*vplus)) )
  
  gamma[kappa==0] = 0 # No growth of fully mature classes
  
  # goes out of stage (size group)
  Fout = gamma*B
  
  # fraction to reproduction
  Repro = (1-kappa)*vplus*B
  
  # Flux into the size group
  #------------------------------
  Fin = 0
  RR = 0
  for (i in 1:p$nGroups) {
    ix = p$ix[[i]] - p$ix[[1]][1] +1                  # growth to
    ixPrev  = c(ix[length(ix)], ix[1:(length(ix)-1)]) # growth from
    Fin[ix] = Fout[ixPrev]
    
    # for reproduction: consider the reproduction success
    RR[i] = sum(Repro[ix]) + Fin[ix[1]]
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
    out$f     = f[-p$ixR] # Feeding level only all fish stages, no resources
    out$mortpred = mortpred
    out$g     = g # net growth rate fish stages
    out$Repro = Repro
    out$Fin   = Fin
    out$Fout  = Fout
    
    # for the budget:
    grazing = p$Cmax * f         # grazing rate, /yr
    loss    = (1.-p$epsAssim) * grazing + p$metabolism
    
    il <- NULL
    for (i in 1:length(p$ix))
      il <- c(il, rep(i, times=length(p$ix[[i]])))
    
    out$totMort    = tapply((mort   *u)[p$ixFish], INDEX=il, FUN=sum)
    out$totGrazing = tapply((grazing*u)[p$ixFish], INDEX=il, FUN=sum)
    out$totLoss    = tapply((loss   *u)[p$ixFish], INDEX=il, FUN=sum)
    out$totRepro   = RR    
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
#' It also allows simulations for customized FEISTY setups.
#'
#' @usage simulateFEISTY (bCust = FALSE, p = setupBasic(), tEnd = 100, tStep  = 1, times = seq(from=0, to=tEnd, by=tStep),
#'        yini = p$u0, USEdll = TRUE, Rmodel = derivativesFEISTYR)
#'
#' @param bCust Logical flag, indicates whether to use fixed setups (FALSE) or customized setups (TRUE). \cr 
#' If \code{bCust} is TURE, FEISTY simulations based on customized setups only can be done in FORTRAN, not R. \code{useDLL} input does not work.
#' @param p A complete parameter list. \cr 
#' The parameter of FEISTY setups can be one of the following: \code{\link{setupBasic}}, \code{\link{setupBasic2}}, \code{\link{setVertical}}, and \code{\link{setVertical2}}.
#' Or, modelers can customize new setups before calling \code{simulateFEISTY}.
#' @param tEnd The end time for the simulation [year], i.e., simulation period of FEISTY, in years. Time (in years) over which the model should be simulated.
#' @param tStep The time step for ODE solving output [year]. Default 1
#' @param times A sequence of time points for FEISTY simulations (ODEs solving), required by function \code{\link{ode}}. Generally, it needs nothing, since it will be generated by `tEnd` and `tStep` automatically. \cr
#'              If `NA`, the function returns only the derivative (one time step running) by `Rmodel`. Default \code{\link{derivativesFEISTYR}}.
#' @param yini A vector containing initial biomass values of all state variables (resources and all size classes). The default is from the setup parameter list, `p$u0`.
#' @param USEdll Logical flag, determining whether the ODEs are solved in FORTRAN (`TRUE`) or R (`FALSE`). \cr
#' The \link{deSolve} package is used for both methods.
#' @param Rmodel The R function for computing derivatives, defaults to \code{\link{derivativesFEISTYR}}. Generally, should not be changed, unless modelers modify the model profoundly.
#'
#' @details
#' The function runs the FEISTY model simulation over the specified time frame. \cr
#' The simulation supports published FEISTY setups and their revised versions:
#' \code{\link{setupBasic}}, \code{\link{setupBasic2}}, \code{\link{setVertical}}, and \code{\link{setVertical2}}, and customized setups by modelers. \cr
#' The simulation can be conducted by either a FORTRAN-based approach or an R-based approach. Both methods rely on the \link{desolve} package for ODE solving. 
#' For efficiency, FORTRAN dll should be used. For model development, the R-version is preferred.
#' Simulations based on customized setups only can be done by the FORTRAN-based approach.
#'
#' @return
#' A list containing the simulation results:
#' \itemize{
#' \item u: a matrix of biomass of each state variable (column) at each time point (row), including resources and all size classes of functional types.
#' \item R: a matrix of biomass of each resource (column) at each time point (row).
#' \item B: a matrix of biomass of each size class (column) at each time point (row).
#' \item t: a vector containing all the simulation time points. From 0 to `tEnd`.
#' \item nTime: The number of time points.
#' \item USEdll: from parameter input.
#' \item p: the parameter list used in the simulation, the same as the input one.
#' \item f: A matrix containing feeding levels [-] of all size classes of functional types over each time point. Resources not included.
#' \item mortpred: A matrix containing a vector containing predation mortality rate [1/year] of all resources and all size classes of functional types over each time point. Resources not included.
#' \item g: A matrix containing the net growth rate [1/year] of all size classes of functional types over each time point. Resources not included.
#' \item Repro: A matrix containing the energy used for reproduction of all size classes of functional types over each time point, rate [gWW/m2/year]. Resources not included.
#' \item Fin: A matrix containing the biomass flux into each size class over each time point [gWW/m2/year]. Resources not included.
#' \item Fout: A matrix containing the biomass flux out of each size class over each time point [gWW/m2/year]. Resources not included.
#' \item totMort: a matrix containing the total mortality [gWW/m2/year] of each functional type over each time point, 
#' which includes predation mortality, background mortality, and fishing mortality.
#' \item totGrazing: a matrix containing the total grazing [gWW/m2/year] of each functional type over each time point, Cmax * f (maximum consumption rate * feeding level)
#' To be simply, the food intake before assimilation.
#' \item totLoss: a matrix containing all biomass loss [gWW/m2/year] of each functional type over each time point, including unassimilated food and metabolism.
#' They are released to environments. where is energy loss from reproduction (1-epsRepro), to be fixed.
#' \item totRepro: a matrix containing the total energy used for reproduction [gWW/m2] of each functional type over each time point.
#' \item totRecruit: a matrix containing the total recruitment [gWW/m2] of each functional type over each time point. 
#' TotRecruit = TotRepro * epsRepro (reproduction efficiency)
#' \item totBiomass: a matrix containing the total biomass [gWW/m2] of each functional type over each time point.
#' \item `SSBAMean`, `SSBGMean`, `SSBMin`, `SSBMax`, and `SSB` can be found in \code{\link{calcSSB}}. \cr
#' `yieldAMean`, `yieldGMean`, `yieldMin`, `yieldMax`, and `yield` can be found in \code{\link{calcYield}}.
#' }
#' 
#' @examples
#' # Just some examples, data input and output may not make sense.
#' 
#' 
#' #-----------------------------------------------
#' # run model with default parameter settings
#' #-----------------------------------------------
#' sim <- simulateFEISTY()
#' 
#' simulateFEISTY(sim)
#' 
#' #colnames(sim)
#' 
#' #plot(sim, which = 1:9) # plot first 9 state variables
#' 
#' #plot(sim, which=c( "smallZoo",  "largeZoo", "smallBenthos", "largeBenthos",
#' #                  "totBiomass.smallPel", "totBiomass.largePel", "totBiomass.demersals"))
#'                   
#' #par(mfrow = c(1, 1))
#' #matplot.0D(sim, type= "l", lty=1, ylab="g/m2", log="y", main="Large Pelagics",
#' #         which=c("largePel_1", "largePel_2", "largePel_3", "totBiomass.largePel"))
#' 
#' # -------------------------------------------------------------------------------
#' 
#' # run FEISTY simulation based on setupVertical
#' # prepare a parameter list
#' p_V <- setupVertical(szprod = 100, lzprod = 120, bent = 200, region = 2, depth = 1000, photic = 150)
#' # run the simulation by R and get the simplified output
#' sim_Vertical_R <- simulateFEISTY(bCust = FALSE, p = p_V, tEnd = 1000, tStep = 1,yini = p_V$u0, USEdll = FALSE)
#' plotSimulation(sim_Vertical_R)
#' 
#' # run FEISTY simulation based on setupBasic2 by Fortran and get the simplified output
#' sim_Basic2_F <- simulateFEISTY(bCust = FALSE, p = setupBasic2(szprod = 90, lzprod = 100, bprod = 15, depth = 500, Tp = 11, Tb = 9, 
#' nStages=9, etaMature=0.25, F=0, etaF=0.05), tEnd = 1000, tStep = 1, USEdll = TRUE)
#' 
#' plotSimulation(sim_Basic2_F)
#' 
#' # -------------------------------------------------------------------------------
#' 
#' # run FEISTY simulation based on a customized set up
#' 
#' # Initialize the parameter list.
#' p_cust <- paramInit()
#' 
#' # add three resources
#' p_cust <- paramAddResource(p_cust,
#'           names= c("smallZoo", "largeZoo", "smallBenthos"),
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
#' # Add fishing mortality. The baseline fishing mortality is 2/year.
#' p_cust <- setFishing(p_cust, F=2, etaF=0.05)
#' 
#' # Add size preference
#' p_cust$theta <- paramSizepref(p = p_cust,   
#'                         beta = 400,
#'                         sigma = 1.3, 
#'                         type = 1)
#' 
#' # run the simulation for 500 years and get the detailed output. 
#' sim_cust <- simulateFEISTY(bCust = TRUE, p = p_cust, tEnd = 500)
#' 
#' plotSimulation(sim)
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
#' @author Ken H. Andersen, Karline Soetaert <karline.soetaert@nioz.nl>, Yixin Zhao
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
#' \code{\link{deSolve}} The package for ODEs solving \cr
#' \code{\link{derivativesFEISTYR}} The derivative function of state variables in FEISTY model
#' 
#' \code{\link{webFEISTY}} A shiny interface for visualizing FEISTY model results
#' 
#' \code{\link{plotSimulation}} Plot simulation results including rates, biomass, and SSB data
#' 
#' @export
#'

# ------------------------------------------------------------------------------
#
# Simulate the model.
#
# In: 
#  bCust : FALSE-> Used fixed setups coded in the Fortran library 
#                    (published and revised).
#          TRUE -> Use customized setups. Customize your own setups by
#                     mimicking the setup configuration in FEISTY_setup.R and FEISTY_parms.R.
#  p : fully populated set of parameters
#
#          setup type:  setupBasic: Petrik et al. (2019)
#                       setupBasic2: as (1) but with adjusted to work with arbitrary no of stages
#                       setVertical: van Denderen et al. (2020)
#                       setVertical2: revised based on van Denderen et al. (2020)  
#
#  tEnd : time to simulate (years)
#  times: The times steps to return. If times=NA then it just does one call, essentially
#          used to just get the derivate and not simulate
#  USEdll : TRUE -> ODE solved by Fortran dll / FALSE -> ODE solved by R
#
# Out:
#  A simulation list
# 
# ------------------------------------------------------------------------------

simulateFEISTY = function(bCust    = FALSE,
                          p      = setupBasic(), 
                          tEnd   = 100,
                          tStep  = 1,
                          times  = seq(from=0, to=tEnd, by=tStep),  
                          yini   = p$u0,  
                          USEdll = TRUE,
                          Rmodel = derivativesFEISTYR) 
{
  
  nR      <- p$nResources[1]  # no of resources. [1] to make sure that this is only one number
  nGroups <- p$nGroups[1] # no of fish groups
  nGrid   <- p$nStages[1] # no of grid points
  nFGrid  <- nGrid-nR # grid points of fish
  
  if (length(yini) != nGrid) 
    stop ("length of 'yini' not ok - should be ", nGrid)  
  
  # prepare rtol atol for ode solving
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
  if (USEdll==TRUE || bCust==TRUE) {    
    # the integers to be passed to the fortran code
    ipar <- c(nGroups,                           # total number of groups
              nR,                                # total number of resources
              unlist(lapply(p$ix, FUN=length)),  # number of stages per fish group
              p$Rtype,                           # type of resource dynamics
              if (is.null(p$pelgrididx)) NULL else length(p$pelgrididx),
              p$pelgrididx,
              if (is.null(p$allgrididx)) NULL else length(p$allgrididx),
              p$allgrididx,
              if (is.null(p$lgdemidx)) NULL else length(p$lgdemidx),
              p$lgdemidx,
              if (is.null(p$bET)) NULL else as.integer(p$bET))
    ipar <- as.integer(ipar)
    if (length(c(nGroups,nR,unlist(lapply(p$ix, FUN=length)),p$Rtype) != 3 + nGroups))
      stop ("length of 'ipar' not ok -check parameters")
    
    if (any(dim(p$theta)-c(nGrid, nGrid) != 0))
      stop ("dimension of 'theta' not ok: should be (", nGrid, ",", nGrid, ")")  
    
    # the double precisions to be passed to the fortran code
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
    
    # names of functions in fortran code to be used
    runfunc  <- "runfeisty"    # the derivative function
    
    # 
    # Run the simulation:
    #
    if (bCust==TRUE){     # for customized setups
      
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
    
    if (bCust==FALSE){     # for fixed setups
      # Oct 2023 Transmit input file path to Fortran library
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

        file_path=system.file("data", "input.nml", package = "FEISTY")
        dummy=.C("passpath", length=nchar(file_path), file_path_in = charToRaw(file_path))
        file_path_V=system.file("data", "tempdata.dat", package = "FEISTY")
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
        setupinput=c(p$szprod,p$lzprod,p$bprodin,p$dfbot,length(p$ix[[p$nGroups]]),p$depth,p$Tp,p$Tb,p$etaMature,p$F,p$etaF,as.integer(p$bET))
      }else if(p$setup=="setupVertical"){
        initfunc <- "initfeistysetupvertical"
        setupinput = c(p$szprod,p$lzprod,p$bprodin,p$dfbot,p$dfpho,p$region, p$bottom, p$photic)
      }else if(p$setup=="setupVertical2"){
        initfunc <- "initfeistysetupvertical2"
        setupinput = c(p$szprod,p$lzprod,p$bprodin,p$dfbot,p$dfpho,length(p$ix[[p$nGroups]]),p$region,p$bottom,p$photic,p$etaMature,p$F,p$etaF)
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
  # Energy used for reproduction [gWW/m2/year]
  col_Repro=grep("^Repro", colnames(u), value = TRUE)
  sim$Repro=u[,col_Repro]
  # Biomass flux into each size class
  col_Fin=grep("^Fin", colnames(u), value = TRUE)
  sim$Fin=u[,col_Fin]
  # Biomass flux out of each size class
  col_Fout=grep("^Fout", colnames(u), value = TRUE)
  sim$Fout=u[,col_Fout]
  # total mortality of each functional group [gWW/m2/year],which includes predation mortality, background mortality, and fishing mortality.
  col_totMort=grep("^totMort", colnames(u), value = TRUE)
  sim$totMort=u[,col_totMort]
  # total grazing of each functional group [gWW/m2/year], Cmax * f (maximum consumption rate * feeding level), the food intake before assimilation.
  col_totGrazing=grep("^totGrazing", colnames(u), value = TRUE)
  sim$totGrazing=u[,col_totGrazing]
  # total biomass loss of each functional group [gWW/m2/year], including unassimilated food and metabolism. They are released to environments. where is energy loss from reproduction (1-epsRepro), to be fixed.
  col_totLoss=grep("^totLoss", colnames(u), value = TRUE)
  sim$totLoss=u[,col_totLoss]
  # total energy used for reproduction of each functional group [gWW/m2]
  col_totRepro=grep("^totRepro", colnames(u), value = TRUE)
  sim$totRepro=u[,col_totRepro]
  # total recruitment of each functional group [gWW/m2], TotRecruit = TotRepro * epsRepro (reproduction efficiency)
  col_totRecruit=grep("^totRecruit", colnames(u), value = TRUE)
  sim$totRecruit=u[,col_totRecruit] 
  # total biomass of each functional group [gWW/m2]
  col_totBiomass=grep("^totBiomass", colnames(u), value = TRUE)
  sim$totBiomass=u[,col_totBiomass]  
  
  #
  # Calculate Spawning Stock Biomass and yield
  #
  sim=calcSSB(sim=sim,etaTime=0.4)
  sim=calcYield(sim=sim,etaTime=0.4)
  
  return(sim)
}
