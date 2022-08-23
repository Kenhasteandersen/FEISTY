#
# Core code for the library
#
library('deSolve')
#
# Initialize the parameter structure
#
# In:
#  depth : depth (meters)
#  pprod : Primary productivity (gww m-3 yr-1)
#

# Out:
#  A parameters list
#
parametersInit = function(depth, pprod) {
  param = list()
  param$depth = depth
  param$pprod = pprod
  
  #
  # Prepare for groups to be added:
  # 
  param$nGroups = 0
  param$ix = list() # Indexes of each fish group
  param$ixFish = list() # Indices to all fish groups
  
  return(param)
}


#new
loadFEISTYmodel = function() {
  sys=Sys.info()['sysname']
  
  if (sys=='Darwin') 
    sLibname = '../lib/libFEISTY.dylib'
  if (sys=='Linux') 
    sLibname = '../lib/libFEISTY.so'
  if (sys=='Windows')
    sLibname = '../lib/libFEISTY.dll'
  
  dyn.load(sLibname)
}





#
# Makes a grid.
#
# In:
#  mMin, mMax : minimum and maximum sizes (gram).
#  nStages : number of stages
#
# Out:
#  A list with lower and upper sizes, the ratio between the two (z) and
#  the center mass (mc).
#
makeGrid = function(mMin, mMax, nStages) {
  #
  # Setup mass grid:
  #
  m = exp(seq(log(mMin), log(mMax), length.out=(nStages+1)) ) # grid of lower sizes, with the last being the upper size of the last cell
  mLower = m[1:nStages]
  mUpper = m[2:(nStages+1)]
  z = mUpper / mLower # The ratio between upper and lower sizes
  mc = exp( log(mLower) + 0.5*(log(z)) ) # Geometric mean center mass
  
  return(list(mLower=mLower, mUpper=mUpper, z=z, mc=mc))
}


# Adds another group to the system. The group is defined by its min and
# max sizes, it's size at maturation, and the number of stages.
#
# In:
#  p - a list of parameters                                                      (It is not clear what is suppose to get into the p, REMY)
#  mMin, mMax - min and max masses (lower mass of the first stage
#               and upper mass of the last stage)
#  nStages - number of stages
#
# Out:
#  And updated parameters list
#
paramAddGroup = function(p ,mMin, mMax, mMature, nStages) {
  p$nGroups = p$nGroups + 1                                                     
  n = p$nGroups
  #
  # Setup mass grid:
  #
  grid = makeGrid(mMin, mMax, nStages)
  if (n==1)
    ixStart = max(p$ixR)+1
  else
    ixStart = max(p$ix[[n-1]])+1
  ix = ixStart:(ixStart+nStages-1)
  
  p$ix[[n]] = ix
  p$mLower[ix] = grid$mLower
  p$mUpper[ix] = grid$mUpper
  p$z[ix] = grid$z
  p$mc[ix] = grid$mc
  #
  # Setup maturation schedule:
  #
  p$mMature[n] = mMature
  p$psiMature[ix] = ( 1 + (p$mc[ix]/mMature)^(-5) )^(-1)
  #p$psiMature[ix[length(ix)]] = 1 # Always fully mature in last stage
  #
  # Zero fishing:
  #
  p$mortF[ix] = 0
  
  p$u0[ix] = 1 # Initial conditiona
  
  p$ixFish = min(p$ix[[1]]) : max(p$ix[[n]])
  p$nStates = max(ix)
  return(p)
}


#
# Make a basic three-species setup as described in Petrik et al (2019): Bottom-up 
# drivers of global patterns of demersal, forage, and pelagic fishes. Progress 
# in Oceanography 176, 102124. doi 10.1016/j.pocean.2019.102124.
#
# Out:
#  An updated parameter list. The list contains:
#   ixResources - the indices to the resources
#   ixGroup - array of nGroups with indices for each group
#   mMature(nGroups) - mass of maturation of each group
#

setupBasic = function(pprod = 100, bprod=5) {
  
  # Initialize the parameters:
  param = parametersInit(0, pprod)
  #
  # Setup resource groups:
  #
  param$ixR = 1:4 # 4 resources: Small - large zoo small - large benthos 
  param$r = c(1, 1, 1, 0) # Resource growth rates g ww/m2/yr
  param$K = c(pprod, pprod, bprod, 0)  # g ww/m2
  param$mc = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500)) # weight central size)
  param$mLower = c(2e-06,0.001, 0.5e-03, 0.25) # weight lower limit)  
  param$mUpper = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500)) # weight central size
  param$u0[param$ixR] = param$K # Initial conditions at carrying capacity
  
  #
  # Add fish groups:
  #
  param = paramAddGroup(param, 0.001, 250, 0.5, 2) # Small pelagics
  param = paramAddGroup(param, 0.001, 125000, 250, 3) # Large pelagics
  param = paramAddGroup(param, 0.001, 125000, 250, 3) # Demersals
  
  #
  # Override the generic psiMature and make only adult classes 50% mature
  #
  param$psiMature = 0*param$psiMature
  param$psiMature[6] = 0.5
  param$psiMature[9] = 0.5
  param$psiMature[12] = 0.5
  #
  # Setup physiology:
  #
  h = 20; # Max. consumption coefficient
  n = -0.25 # Max. consumption exponent
  k = 0.011*365 # Metabolism coefficient
  p = -0.175 # Metabolism exponent
  gamma = 70 # Coef. for clearance rate
  q = -0.2 # Clearance rate exponent
  ix = param$ixFish
  m = param$mc[ix]
  param$Cmax[ix] = h*m^n # maximum consumption rate 
  param$V[ix] = gamma*m^q # clearance rate 
  param$metabolism[ix] = k*m^p # 0.2*param$Cmax[ix] # standard metabolism 
  param$epsRepro = rep(0.01, param$nGroups) # reproduction * recruitment efficiency 
  param$epsAssim = 0.7 # Assimilation efficiency
  param$Cmax[is.na(param$Cmax)] = 0
  param$V[is.na(param$V)] = 0 
  
  #
  # Setup size interaction matrix:
  #
  thetaS = 0.25 # Medium fish pref for small zooplankton 
  thetaA = 0.5  # Large fish pref for medium forage fish
  thetaD = 0.75 # Pref of large demersal on pelagic prey
  
  param$theta = matrix(nrow=param$nStates, ncol=param$nStates, data=0)
  
  # Small pelagics:
  param$theta[5,1] = 1 # Small ones eat only small zooplankton
  param$theta[6,1:10] = c(thetaS, 1, 0, 0, 1, 0, 1, 0,0, 1) 

  # Large pelagics:
  param$theta[7,1] = 1
  param$theta[8,1:10] = c(thetaS, 1, 0, 0, 1, 0, 1, 0,0, 1) 
  param$theta[9,6] = thetaA # medium forage fish
  param$theta[9,8] = 1 # medium large pelagics
  
  # Demersals:
  param$theta[10,1] = 1
  param$theta[11,3] = 1 # Medium demersals eat benthos
  param$theta[12,6] = thetaA*thetaD  # medium forage fish
  param$theta[12,8] = thetaD # medium large pelagics
  param$theta[12,11] = 1 # medium demersals
  
  
  
  # beta = 400
  # sigma = 1.3
  # param$theta = matrix(nrow=param$nStates, ncol=param$nStates)
  # for (i in param$ixFish) {
  #   param$theta[i,] = exp( -(log(param$mc[i]/(beta*param$mc)))^2 / (2*sigma)^2  )
  #   param$theta[i,param$mc>param$mc[i]] = 0
  # }
  # param$theta[is.na(param$theta)] = 0
  # 
  
  # # Setup interactions between groups and resources:
  # #
  # ixSmall = param$ix[[1]]
  # ixLarge = param$ix[[2]]
  # ixDem = param$ix[[3]]
  # 
  # mMedium = 10
  # mLarge = 5000
  # # Pelagic/demersal indices:
  # ixR = param$ixR
  # #ixDemPelagic = ixDem[(param$mc[ixDem] < 10) | (param$mc[ixDem] > 5000)]
  # #ixDemDemersal = ixDem[param$mc[ixDem] >= 10]
  # #ixLargeDemersal = ixLarge[param$mc[ixLarge] > 5000]
  # #ixPelagic = c(param$ixR[1],param$ixR[2], ixSmall, ixLarge, ixDemPelagic)
  # #ixDemersal = c(param$ixR[3],param$ixR[4], ixDemDemersal, ixLargeDemersal)
  # # Small pelagics feed on pelagic resources:
  # 
  # for (i in ixSmall) {
  #   param$theta[i,ixR[3:4]] = 0 
  #   param$theta[i,ixDem] = 0
  # }
  # # Large pelagics feed on pelagic resources:
  # for (i in ixLarge) {
  #   param$theta[i,ixR[3:4]] = 0 
  #   param$theta[i,ixDem[ (param$mc[ixDem]>mMedium) & 
  #                          (param$mc[ixDem]<mLarge) ]] = 0
  # }
  # # Small demersals feed on pelagic resources
  # for (i in ixDem[param$mc[ixDem]<mMedium])
  #   param$theta[i,ixR[3:4]] = 0
  # for (i in ixDem[ (param$mc[ixDem]>mMedium) & 
  #                    (param$mc[ixDem]<mLarge) ]) {
  #   param$theta[i,ixR[1:2]] = 0
  #   param$theta[i,param$ixFish] = 0
  # }
  
# LArge demersals feed on benthic resources and all fish:
#  for (i in ixDemersal[ param$mc[ixDemersal]>mLarge ]) {
#    param$
#  } 

  #
  # Mortality
  #
  param$mort0 = 0.1
  param$mortF[param$ixFish] = 0.3#*c(0,1,0,0,1,0,0,1) # Fishing only on mature stages
  

  param$metabolism[is.na(param$metabolism)]=0
  param$mortF[is.na(param$mortF)]=0
  param$psiMature[is.na(param$psiMature)]=0
  param$z[is.na(param$z)]=0

  return(param)
}

#
# Make a basic three-species setup based up setupBasic(), but generalised to:
# - more realistic sizes
# - Generalized size-based feeding
# - possiblity of more then 3 size groups in each group
#
# Out:
#  An updated parameter list. The list contains:
#   ixResources - the indices to the resources
#   ixGroup - array of nGroups with indices for each group
#   mMature(nGroups) - mass of maturation of each group
#

setupBasic2 = function(pprod = 100, bprod=5, nSizeGroups=9) {
  
  # Initialize the parameters:
  param = parametersInit(0, pprod)
  #
  # Setup resource groups:
  #
  param$ixR = 1:4 # 4 resources: Small - large zoo small - large benthos 
  param$r = c(1, 1, 1, 0) # Resource growth rates g ww/m2/yr
  param$K = c(pprod, pprod, bprod, 0)  # g ww/m2
  param$mc = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500)) # weight central size)
  param$mLower = c(2e-06,0.001, 0.5e-03, 0.25) # weight lower limit)  
  param$mUpper = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500)) # weight central size
  param$u0[param$ixR] = param$K # Initial conditions at carrying capacity
  
  #
  # Add fish groups:
  #
  param = paramAddGroup(param, 0.001, 250, 0.5, round(0.66*nSizeGroups)) # Small pelagics
  param = paramAddGroup(param, 0.001, 125000, 250, nSizeGroups) # Large pelagics
  param = paramAddGroup(param, 0.001, 125000, 250, nSizeGroups) # Demersals
  #
  # Setup physiology:
  #
  h = 20; # Max. consumption coefficient
  n = -0.25 # Max. consumption exponent
  k = 0.011*365 # Metabolism coefficient
  p = -0.175 # Metabolism exponent
  gamma = 70 # Coef. for clearance rate
  q = -0.2 # Clearance rate exponent
  ix = param$ixFish
  m = param$mc[ix]
  param$Cmax[ix] = h*m^n # maximum consumption rate 
  param$V[ix] = gamma*m^q # clearance rate 
  param$metabolism[ix] = k*m^p # 0.2*param$Cmax[ix] # standard metabolism 
  param$epsRepro = rep(0.01, param$nGroups) # reproduction * recruitment efficiency 
  param$epsAssim = 0.7 # Assimilation efficiency
  param$Cmax[is.na(param$Cmax)] = 0
  param$V[is.na(param$V)] = 0 
  
  #
  # Setup size interaction matrix:
  #
  thetaA = 0.5  # Large fish pref for medium forage fish
  thetaD = 0.75 # Pref of large demersal on pelagic prey
  
  param$theta = matrix(nrow=param$nStates, ncol=param$nStates, data=0)
  
  # # Small pelagics:
  # param$theta[5,1] = 1 # Small ones eat only small zooplankton
  # param$theta[6,1:10] = c(thetaS, 1, 0, 0, 1, 0, 1, 0,0, 1) 
  # 
  # # Large pelagics:
  # param$theta[7,1] = 1
  # param$theta[8,1:10] = c(thetaS, 1, 0, 0, 1, 0, 1, 0,0, 1) 
  # param$theta[9,6] = thetaA # medium forage fish
  # param$theta[9,8] = 1 # medium large pelagics
  # 
  # # Demersals:
  # param$theta[10,1] = 1
  # param$theta[11,3] = 1 # Medium demersals eat benthos
  # param$theta[12,6] = thetaA*thetaD  # medium forage fish
  # param$theta[12,8] = thetaD # medium large pelagics
  # param$theta[12,11] = 1 # medium demersals
  
   #
   # Size-based interactions:  
   #
   beta = 400
   sigma = 1.3
   param$theta = matrix(nrow=param$nStates, ncol=param$nStates)
   for (i in param$ixFish) {
     param$theta[i,] = exp( -(log(param$mc[i]/(beta*param$mc)))^2 / (2*sigma)^2  )
     param$theta[i,param$mc>param$mc[i]] = 0
   }
   param$theta[is.na(param$theta)] = 0
  
   #
   # Setup interactions between groups and resources:
   #
   ixSmall = param$ix[[1]]
   ixLarge = param$ix[[2]]
   ixDem = param$ix[[3]]
  
   mMedium = 10
   mLarge = 5000
   ixSmallSizeDem = ixDem[ (param$mc[ixDem]<=mMedium) ]
   ixMediumSizeDem = ixDem[ (param$mc[ixDem]>mMedium) &
     (param$mc[ixDem]<mLarge) ]
   # Pelagic/demersal indices:
   ixR = param$ixR
   #ixDemPelagic = ixDem[(param$mc[ixDem] < 10) | (param$mc[ixDem] > 5000)]
   #ixDemDemersal = ixDem[param$mc[ixDem] >= 10]
   #ixLargeDemersal = ixLarge[param$mc[ixLarge] > 5000]
   #ixPelagic = c(param$ixR[1],param$ixR[2], ixSmall, ixLarge, ixDemPelagic)
   #ixDemersal = c(param$ixR[3],param$ixR[4], ixDemDemersal, ixLargeDemersal)
   
   # Pelagic fish do not feed on benthic resources
   param$theta[ixSmall, 3:4] = 0
   param$theta[ixLarge, 3:4] = 0
   # ... or on medium-sized demersal fish:
   param$theta[ixSmall, ixMediumSizeDem] = 0
   param$theta[ixLarge, ixMediumSizeDem] = 0
   
   # Large pelagics have reduced feeding efficiency on small pelagics:
   param$theta[ixLarge,ixSmall] = thetaA * param$theta[ixLarge,ixSmall] 
   # ... and do not feed on medium-sized demersal:
   param$theta[ixLarge, ixMediumSizeDem ] = 0
   
   # Medium-sized large demersals feed only on benthos:
   param$theta[ixMediumSizeDem, ixSmall] = 0 
   param$theta[ixMediumSizeDem, ixLarge] = 0 
   
   # Large demersals feed have reduced feeding effiiency on pelagic species:
   param$theta[ixDem, ixSmall] = thetaA * thetaD * param$theta[ixDem,ixSmall] 
   param$theta[ixDem, ixLarge] = thetaD * param$theta[ixDem, ixLarge] 

  #
  # Mortality
  #
  param$mort0 = 0.1
  param$mortF[param$ixFish] = 0.#*c(0,1,0,0,1,0,0,1) # Fishing only on mature stages
  
  
  param$metabolism[is.na(param$metabolism)]=0
  param$mortF[is.na(param$mortF)]=0
  param$psiMature[is.na(param$psiMature)]=0
  param$z[is.na(param$z)]=0
  
  return(param)
}
#
# Make a basic setup with just pelagic fish
#
# Out:
#  An updated parameter list. The list contains:
#   ixResources - the indices to the resources
#   ixGroup - array of nGroups with indices for each group
#   mMature(nGroups) - mass of maturation of each group
#

setupPelagicSpecies = function(depth = 500, pprod = 100, nStages=6, mInf) {
  # Initialize the parameters:x
  param = parametersInit(depth, pprod)
  
  #
  # Setup resource groups:
  #
  param$ixR = 1:4 # 4 resources: Small - large zoo small - large benthos 
  param$r = c(1, 1, 1, 0) # Resource growth rates g ww/m2/yr
  param$K = c(pprod, pprod, 5, 0)  # g ww/m2
  param$mc = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500)) # weight central size)
  param$mLower = c(2e-06,0.001, 0.5e-03, 0.25) # weight lower limit)  
  param$mUpper = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500)) # weight central size
  param$u0[param$ixR] = param$K # Initial conditions at carrying capacity
  param$mortF[param$ixR] = 0 # No fishing on resources
  #
  # Add fish groups:
  #
  for (iGroup in 1:length(mInf))
    param = paramAddGroup(param, 0.001, mInf[[iGroup]], 0.25*mInf[[iGroup]], nStages)
  #
  # Setup physiology:
  #
  h = 20; # Max. consumption coefficient
  n = -0.25 # Metabolic exponent
  q = -0.2 # Clearance rate exponent
  gamma = 70 # Coef. for clearance rate
  ix = param$ixFish
  m = param$mc[ix]
  param$Cmax[ix] = h*m^n # maximum consumption rate 
  param$V[ix] = gamma*m^q # clearance rate 
  param$metabolism[ix] = 0.2*param$Cmax[ix] # standard metabolism 
  param$epsRepro = rep(0.01, param$nGroups) # reproduction * recruitment efficiency 
  param$epsAssim = 0.7 # Assimilation efficiency
  param$Cmax[is.na(param$Cmax)] = 0
  param$V[is.na(param$V)] = 0 
  
  #
  # Setup size interaction matrix:
  #
  beta = 400
  sigma = 1.3
  param$theta = matrix(nrow=param$nStates, ncol=param$nStates)
  for (i in param$ixFish) {
    param$theta[i,] = exp( -(log(param$mc[i]/(beta*param$mc)))^2 / (2*sigma)^2  )
    param$theta[i,param$mc>param$mc[i]] = 0
  }
  param$theta[is.na(param$theta)] = 0
  
  #
  # Setup interactions between groups and resources:
  #
  
  mMedium = 10
  mLarge = 5000
  # Pelagic/demersal indices:
  ixR = param$ixR
  #ixDemPelagic = ixDem[(param$mc[ixDem] < 10) | (param$mc[ixDem] > 5000)]
  #ixDemDemersal = ixDem[param$mc[ixDem] >= 10]
  #ixLargeDemersal = ixLarge[param$mc[ixLarge] > 5000]
  #ixPelagic = c(param$ixR[1],param$ixR[2], ixSmall, ixLarge, ixDemPelagic)
  #ixDemersal = c(param$ixR[3],param$ixR[4], ixDemDemersal, ixLargeDemersal)
  # Small pelagics feed on pelagic resources:
  

  param$theta[,ixR[3:4]] = 0  # No demersal feeding
  #
  # Mortality
  #
  param$mort0 = .5 # NOTE: set pretty high to give a stable population
  
  return(param)
}


#
# Calculate the derivatives of all state variables by R
#
# In:
#  p : a list of parameters
#  u : all state variables
#  bFullOutput : if TRUE returns all internal calculations
#
# Out:
#  deriv : a vector of derivatives
#
calcDerivativesR = function(t, u, p, bFullOutput=FALSE) {
  ix = p$ixFish
  B = u[ix]
  #
  # Calc consumption of all fish groups
  #

  # Calc Encounter:
  Enc = p$V * (p$theta %*% u)
  f = Enc / (p$Cmax + Enc) # Functional response
  f[is.na(f)] = 0
  Eavail = p$epsAssim * p$Cmax * f - p$metabolism

  # Predation mortality:
  #mortpred = t(p$theta) %*% (f*p$Cmax/p$epsAssim*u/p$mc)
  mm = p$Cmax*p$V/(Enc+p$Cmax)*u
  mm[ is.na(mm) ] = 0
  mortpred = t(p$theta) %*% mm

  # Total mortality
  mort = mortpred + p$mort0 + p$mortF
  #
  # Derivative of fish groups
  #

  # Flux out of the size group:
  v = Eavail[ix]
  vplus = pmax(0,v)
  kappa = 1-p$psiMature[ix]
  g=kappa*vplus
  gamma = (kappa*vplus - mort[ix]) /
    (1 - (1/p$z[ix])^(1-mort[ix]/(kappa*vplus)) )
  gamma[kappa==0] = 0 # No growth of fully mature classes
  Fout = gamma*B
  Repro = (1-kappa)*vplus*B

  # Flux into the size group
  Fin = 0
  for (i in 1:p$nGroups) {
    ix = p$ix[[i]] - p$ix[[1]][1] +1
    ixPrev = c(ix[length(ix)], ix[1:(length(ix)-1)])
    Fin[ix] = Fout[ixPrev]
    # Reproduction:
    Fin[ix[1]] = p$epsRepro[i]*(Fin[ix[1]] + sum( Repro[ix] ))
  }

  # Assemble derivatives of fish:
  dBdt = Fin - Fout + (v - mort[p$ixFish])*B - Repro

  #
  # Resources:
  #
  R = u[p$ixR]
  dRdt = p$r*(p$K-R) - mortpred[p$ixR]*R

  if (bFullOutput) {
    out = list()
    out$deriv = c(dRdt, dBdt)
    out$f = f # Feeding level
    out$Repro = Repro
    out$Fin = Fin
    out$Fout = Fout
    out$mortpred = mortpred
    out$mort = mort
    out$g=g
    return(out)
  }
  else
    return( list(c(dRdt, dBdt)) )
}

#
# Calculate the derivatives of all state variables by Fortran dll
#
calcDerivativesF = function(t, y, p, bFullOutput=FALSE) {
  derivF = .Fortran("f_calcderivatives", 
                    u= as.numeric(y),
                    #dt=as.numeric(0),
                    dudt=as.numeric(y))
  
  dRdt = derivF$dudt[1:4]
  dBdt= derivF$dudt[p$ixFish]

   if (bFullOutput) {
     nGrid = p$nStates
 out = .Fortran("f_getrates", 
                 u= as.numeric(y), # y is u last step, defined from plotSimulation
                dudt=as.numeric(y), # in dll initially set as 0 and then runs
                flvl_r=vector(length=nGrid, mode="numeric"),
                mortpred_r = vector(length=nGrid, mode="numeric"),
                g_r = vector(length=nGrid-length(p$ixR), mode="numeric")
                )                  
                 
     out$f = out$flvl_r # Feeding level
     out$mortpred = out$mortpred_r
     out$g = out$g_r
     return(out)
   }
   else
  return(derivF$dudt)
}


#
# Simulate the model.
#
# In: 
#  p : fully populated set of parameters
#  tEnd : time to simulate (years)
#  USEdll : TRUE -> ODE solved by dll / FALSE -> ODE solved by R
#
# Out:
#  A simulation list
# 
simulate= function(p = setupBasic(), tEnd = 100,USEdll=TRUE) {
  #
  # Integrate the equations:
  #
  start_time = Sys.time()
  t = seq(0, tEnd,length.out=tEnd+1) 
  #Calc by R or dll
  if (USEdll) {
    loadFEISTYmodel()
    dummy = .Fortran("f_setupFEISTY", pprod=as.double(p$pprod))
    #dudt = assign("dudt", rep(as.double(0),12), envir = .GlobalEnv) 
    u = ode(y=p$u0,
            times = t,
            func = function(t,y,parms) list(calcDerivativesF(t,y,parms)),# Run by dll
            parms = p)
  }else{
    u = ode(y = p$u0,
            times = t,
            func = calcDerivativesR, #Run by R
            parms = p)
}
  
  #
  # Assemble output:
  #
  sim = list()
  p$USEdll=USEdll
  sim$p = p
  sim$R = u[, p$ixR+1]
  sim$B = u[, p$ixFish+1]
  sim$t = t
  sim$nTime = length(t)
  
  #
  # Calculate Spawning Stock Biomass
  #
  SSB = matrix(nrow=length(t), ncol=p$nGroups)
  yield = SSB
  for (i in 1:p$nGroups)
    for (j in 1:length(t)) {
      SSB[j,i] = sum( u[j, p$ix[[i]]] * p$psiMature[p$ix[[i]]] )
      yield[j,i] = sum( u[j, p$ix[[i]]] * p$mortF[p$ix[[i]]] )
    }
  sim$SSB = SSB
  sim$yield = yield
  
  end_time = Sys.time()
  sim$tictoc=end_time-start_time
  print(sim$tictoc)
  return(sim)
  
}

