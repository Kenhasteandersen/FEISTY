#===============================================================================
# The Feisty model as described in 
# Van Denderen et al. 2020 Global Ecology and Biogeography, DOI: 10.1111/geb.13348 
# Emergent global biogeography of marine fish food webs.

# Petrik, CM, Stock, CA, Andersen, KH, van Denderen, PD, Watson, JR 2019. 
# Bottom-up drivers of global patterns of demersal, forage, and pelagic fishes. 
# Progress in Oceanography, 176, 102124. https://doi.org/10.1016/j.pocean.2019.102124

# Slightly rewritten by Karline Soetaert, based on code from Ken H. Andersen
# 
# Main model routines
#      derivativesFeistyR  : derivative code in R
#      derivativesFeistyF  : derivative code in Fortran
#      simulateFeisty      : run the feisty model
#
#===============================================================================

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

derivativesFeistyR = function(t,              # current time
                              u,              # all state variables
                              p,              # parameters
                              FullOutput=TRUE) {
  
  # split state variable vector into resource and fish
  u[u<0]=0
  R     = u[p$ixR]       # resource, prey
  iFish = p$ixFish
  B     = u[iFish]       # fish, 
  
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
    out$f     = f # Feeding level all stages
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
    
    out$TotMort    = tapply((mort   *u)[p$ixFish], INDEX=il, FUN=sum)
    out$TotGrazing = tapply((grazing*u)[p$ixFish], INDEX=il, FUN=sum)
    out$TotLoss    = tapply((loss   *u)[p$ixFish], INDEX=il, FUN=sum)
    out$TotRepro   = RR    
    out$TotRecruit = out$TotRepro* p$epsRepro
    out$Fishes     = tapply(B, INDEX=il, FUN=sum)
    return(out)
  }
  else # Output just the derivatives
    return( list(c(dRdt, dBdt)) )
}


# ------------------------------------------------------------------------------
#
# Simulate the model.
#
# In: 
#  bCust : FALSE-> Used fixed setups coded in the Fortran library 
#                    (published and revised).
#          TRUE -> Use customized setups. Customize your own setups by
#                     mimicking the setup configuration in Feisty_setup.R and Feisty_parms.R.
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
#  simpleOutput: boolean that states whether the full output is needed or not
#
# Out:
#  A simulation list
# 
# ------------------------------------------------------------------------------

simulateFeisty = function(bCust    = FALSE,
                          p      = setupBasic(), 
                          tEnd   = 100,
                          tStep  = 1,
                          times  = seq(from=0, to=tEnd, by=tStep),  
                          yini   = p$u0,  
                          USEdll = TRUE,
                          Rmodel = derivativesFeistyR,
                          simpleOutput = FALSE) 
{
  
  nR      <- p$nResources[1]  # no of resources. [1] to make sure that this is only one number
  nGroups <- p$nGroups[1] # no of fish groups
  nGrid   <- p$nStages[1] # no of grid points
  nFGrid  <- nGrid-nR # grid points of fish
  
  if (length(yini) != nGrid) 
    stop ("length of 'yini' not ok - should be ", nGrid)  
  
  # prepare rtol atol for ode
    rtol = 1E-8
    atol = 1E-8
  if (max(sapply(p$ix, length))>=21){    
    rtol = 1E-10
    atol = 1E-10}
  if (max(sapply(p$ix, length))>27) stop("The size number cannot be more than 27 due to the low accuracy of integration.")
  
  #
  # calculate in Fortran
  #
  if (USEdll) {    
    # the integers to be passed to the fortran code
    ipar <- c(nGroups,                           # total number of groups
              nR,                                # total number of resources
              unlist(lapply(p$ix, FUN=length)),  # number of stages per fish group
              p$Rtype)                           # type of resource dynamics
    ipar <- as.integer(ipar)
    if (length(ipar) != 3 + nGroups) 
      stop ("length of 'ipar' not ok -check parameters")  
    
    if (any(dim(p$theta)-c(nGrid, nGrid) != 0))
      stop ("dimension of 'theta' not ok: should be (", nGrid, ",", nGrid, ")")  
    
    # the double precisions to be passed to the fortran code
    rpar   <- c(rep(p$K,  length.out=nR),            # resource parameters
                rep(p$r,  length.out=nR),  
                rep(p$epsRepro, length.out=nGroups), # group-specific parameter
                p$psiMature[-(1:nR)],                # fish-stage parameter
                p$z[-(1:nR)], 
                t(p$theta),                              # check if not transpose
                rep(p$epsAssim,   length.out=nGrid), # all  
                rep(p$V,          length.out=nGrid), 
                rep(p$Cmax,       length.out=nGrid),
                rep(p$metabolism, length.out=nGrid),
                rep(p$mort0,      length.out=nGrid),
                rep(p$mortF,      length.out=nGrid))  
    
    # names of functions in fortran code to be used
    runfunc  <- "runfeisty"    # the derivative function
    
    # output variables
    Sname <- p$stagenames
    Fname <- p$stagenames[-(1:nR)]
    Gname <- p$groupnames[-(1:nR)]
    outnames <- c(
      paste("f", Sname, sep="."), paste("mortpred", Sname, sep="."),
      paste("g", Fname, sep="."), paste("Repro", Fname, sep="."),
      paste("Fin", Fname, sep="."), paste("Fout", Fname, sep="."),
      paste("totMort", Gname, sep="."), paste("totGrazing", Gname, sep="."),
      paste("totLoss", Gname, sep="."), paste("totRepro", Gname, sep="."),
      paste("totRecruit", Gname, sep="."), paste("totBiomass", Gname, sep="."))
    
    # 
    # Run the simulation:
    #
    if (bCust==TRUE){     # If using the R code for calculating the derivative (slow):
      
      initfunc <- "initfeisty"
      
      if (any(is.na(times)))  # one call and return
        return( DLLfunc(y=yini, times=0, parms=NULL, dllname = "FeistyR",
                        func=runfunc, initfunc=initfunc, outnames=outnames, nout=length(outnames),
                        ipar=ipar, rpar=as.double(rpar)))
      
      u = ode(y=yini, times=times, parms=NULL, dllname = "FeistyR",
              func=runfunc, initfunc=initfunc, outnames=outnames, nout=length(outnames),
              ipar=ipar, rpar=as.double(rpar),
              method = "ode45", rtol = rtol, atol = atol) # Run by dll
    }
    
    if (bCust==FALSE){    # If using the derivative from Fortran (fast):
      # Oct 2023 Transmit input file path to Fortran library
      passpath <- function() {
        sys=Sys.info()['sysname']
        
        if (sys=='Darwin') {
          sLibname = system.file("libs", "FeistyR.so", package = "FeistyR")
        }
        if (sys=='Linux') {
          sLibname = system.file("libs", "FeistyR.so", package = "FeistyR")
        }
        if (sys=='Windows'){
          if (Sys.info()['machine']=='x86-64'){
            sLibname = system.file("libs/x64", "FeistyR.dll", package = "FeistyR")
          }else{
            sLibname = system.file("libs/i386", "FeistyR.dll", package = "FeistyR")
          }
        }
        
        if(!is.loaded(sLibname)) dyn.load(sLibname)
        
        file_path=system.file("data", "input.nml", package = "FeistyR")
        dummy=.C("passpath", length=nchar(file_path), file_path_in = charToRaw(file_path))
        file_path_V=system.file("data", "tempdata.dat", package = "FeistyR")
        dummy=.C("passpathv", length=nchar(file_path_V), file_path_in=charToRaw(file_path_V))
      }
      
      # Call the Fortran subroutine to pass input file path
      passresult <- passpath()
      
      # Choose the setup:
      if (p$setup=="setupBasic"){
        initfunc <- "initfeistysetupbasic"
        setupinput=c(p$szprod,p$lzprod,p$bprod,p$depth,p$Tp,p$Tb)
      }else if(p$setup=="setupBasic2"){
        initfunc <- "initfeistysetupbasic2"
        setupinput=c(p$szprod,p$lzprod,p$bprod,length(p$ix[[p$nGroups]]),p$depth,p$Tp,p$Tb,p$etaMature,p$F,p$etaF)
      }else if(p$setup=="setupVertical"){
        initfunc <- "initfeistysetupvertical"
        setupinput = c(p$szprod,p$lzprod,p$bent,p$region, p$bottom, p$photic)
      }else if(p$setup=="setupVertical2"){
        initfunc <- "initfeistysetupvertical2"
        setupinput = c(p$szprod,p$lzprod,p$bent,length(p$ix[[p$nGroups]]),p$region,p$bottom,p$photic,p$etaMature,p$F,p$etaF)
      }
      
      if (any(is.na(times)))  # one call and return
        return( DLLfunc(y=yini, times=0, parms=as.double(setupinput), dllname = "FeistyR",
                        func=runfunc, initfunc=initfunc, outnames=outnames, nout=length(outnames),
                        ipar=NULL, rpar=NULL))
      # Full simulation:
      u = ode(y=yini, times=times, parms=as.double(setupinput), dllname = "FeistyR",
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
  }
  
  # if (any(u[,c(p$ixR,p$ixFish)+1]<0))
  #   stop (paste("Negative biomass occurs! The current timestep length is tStep = ", tStep,
  #               ". Try shortening the timestep length to increase accuracy, e.g., sim = simulateFeisty(..., tStep = ", tStep/10, ").",sep=""))
  
  #
  # Assemble output:
  #
  if (! simpleOutput) return(u)
 
  sim   = list()
  sim$u = u[,c(p$ixR,p$ixFish)+1]
  sim$R = u[, p$ixR+1]
  sim$B = u[, p$ixFish+1]
  sim$t = times
  sim$nTime = length(times)
  p$USEdll=USEdll
  sim$p = p
  
  #
  # Calculate Spawning Stock Biomass and yield
  #
  sim=calcSSB(sim=sim,etaTime=0.6)
  sim=calcYield(sim=sim,etaTime=0.6)
  
  return(sim)
}
