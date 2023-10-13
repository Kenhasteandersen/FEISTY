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

# Load deSolve:
library(deSolve)

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

  if (FullOutput) {
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
  else
    return( list(c(dRdt, dBdt)) )
}


# ------------------------------------------------------------------------------
#
# Simulate the model.
#
# In: 
#  cus : FALSE-> Fixed setups (published and revised).
#        TRUE -> Customize your own setup. 'setup' and 'setupini' are revoked.
#  setup : It only works when 'cus' is TRUE.
#          setup type: 1 = Petrik et al. (2019)
#                      2 = 
#                      3 = van Denderen et al. (2020)
#                      4 =           
#  setupini : It only works when 'cus' is TRUE.......
#
#  p : fully populated set of parameters
#  tEnd : time to simulate (years)
#  USEdll : TRUE -> ODE solved by dll / FALSE -> ODE solved by R
#
# Out:
#  A simulation list
# 
# ------------------------------------------------------------------------------

simulateFeisty = function(cus    = FALSE,
                          setup  = 1,
                          setupini = c(100,100,5,10,8),# setupbasic(smzprod,lgzprod,bprod,Ts,Tb)
                          p      = setupBasic(), 
                          tEnd   = 100,
                          times  = seq(from=0, to=tEnd, by=1),
                          yini   = p$u0,  
                          USEdll = TRUE,
                          Rmodel = derivativesFeistyR,
                          simpleOutput = FALSE) 
{

  nR      <- p$nResources[1]  # [1] to make sure that this is only one number
  nGroups <- p$nGroups[1]
  nGrid   <- p$nStages[1]
  nFGrid  <- nGrid-nR
  
  if (length(yini) != nGrid) 
      stop ("length of 'yini' not ok - should be ", nGrid)  

  if (USEdll) {    # calculate in Fortran

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
    #initfunc <- "initfeisty"   # the initialisation function
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

if (cus==TRUE){
     initfunc <- "initfeisty"
  
     if (any(is.na(times)))  # one call and return
     return( DLLfunc(y=yini, times=0, parms=NULL, dllname = "FeistyR",
                  func=runfunc, initfunc=initfunc, outnames=outnames, nout=length(outnames),
                  ipar=ipar, rpar=as.double(rpar)))
  
   u = ode(y=yini, times=times, parms=NULL, dllname = "FeistyR",
          func=runfunc, initfunc=initfunc, outnames=outnames, nout=length(outnames),
          ipar=ipar, rpar=as.double(rpar)) # Run by dll
}
    
if (cus==FALSE){
  
  file_path=system.file("data", "input.nml", package = "FeistyR")
  # Call the Fortran subroutine to pass input file path
  result <- passpath(file_path)
  
  if (setup==1){
    initfunc <- "initfeistysetupbasic"
    if (any(is.na(times)))  # one call and return
     return( DLLfunc(y=yini, times=0, parms=NULL, dllname = "FeistyR",
            func=runfunc, initfunc=initfunc, outnames=outnames, nout=length(outnames),
            ipar=NULL, rpar=NULL))
    
    
    u = ode(y=yini, times=times, parms=as.double(setupini), dllname = "FeistyR",
            func=runfunc, initfunc=initfunc, outnames=outnames, nout=length(outnames),
            ipar=NULL, rpar=NULL) # Run by dll
  }
  
  if (setup==2){ #setupbasic2
    
  }
}    
    
    
  } else if (any(is.na(times))) {  # one call and return
    return (Rmodel(0, yini, p))
  } else {               # R-code
    
    u = ode(y=yini, times=times, parms=p, func = Rmodel) #Run by R

  }
  if (! simpleOutput) return(u)
  #
  # Assemble output:
  #
  sim   = list()
  sim$R = u[, p$ixR+1]
  sim$B = u[, p$ixFish+1]
  sim$t = times
  sim$nTime = length(times)
  p$USEdll=USEdll
  sim$p = p
  
  #
  # Calculate Spawning Stock Biomass
  #
  SSB = matrix(nrow=length(times), ncol=p$nGroups)
  yield = SSB
  for (i in 1:p$nGroups){
      SSB[,i]   = colSums( t(u[, p$ix[[i]]]) * p$psiMature[p$ix[[i]]] )
      yield[,i] = colSums( t(u[, p$ix[[i]]]) * p$mortF[p$ix[[i]]] )
  }
  sim$SSB = SSB
  sim$yield = yield
  
  return(sim)
  
}

# Oct 2023 Transmit input file path to Fortran library
passpath <- function(file_path) {
  
   sys=Sys.info()['sysname']
   
   if (sys=='Darwin') {
    # sLibname = system.file("libs", "FeistyR.dll", package = "FeistyR")
   }
   if (sys=='Linux') {
    # sLibname = system.file("libs", "FeistyR.so", package = "FeistyR")
   }
   if (sys=='Windows'){
       if (Sys.info()['machine']=='x86-64'){
         sLibname = system.file("libs/x64", "FeistyR.dll", package = "FeistyR")
       }else{
         sLibname = system.file("libs/i386", "FeistyR.dll", package = "FeistyR")
       }
   }
   
  dyn.load(sLibname)
  
  .Fortran("passpath", file_path_in = as.character(file_path))
}
