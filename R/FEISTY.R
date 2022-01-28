#
# Core code for the library
#
library('deSolve')
#
# Initialize the parameter structure
#
# In:
#  depth : depth (meters)
#  pprod : Primary productivity (UNITS?)
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
  m = exp( seq(log(mMin), log(mMax), length.out=(nStages+1)) ) # grid of lower sizes, with the last being the upper size of the last cell
  mLower = m[1:nStages]
  mUpper = m[2:(nStages+1)]
  z = mUpper / mLower # The ratio between upper and lower sizes
  mc = exp( log(mLower) + 0.5*(log(z)) ) # Geometric mean center mass
  
  return( list(mLower=mLower, mUpper=mUpper, z=z, mc=mc))
}
#
# Adds another group to the system. The group is defined by its min and
# max sizes, it's size at maturation, and the number of stages.
#
# In:
#  p - a list of parameters
#  mMin, mMax - min and max masses (lower mass of the first stage
#               and upper mass of the last stage)
#  nStages - number of stages
#
# Out:
#  And updated parameters list
#
paramAddGroup = function(p,mMin, mMax, mMature, nStages) {
  p$nGroups = p$nGroups + 1
  n = p$nGroups
  #
  # Setup mass grid:
  #
  grid = makeGrid(mMin,mMax,nStages)
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
  
  p$u0[ix] = 1 # Initial conditiona
  
  p$ixFish = min(p$ix[[1]]) : max(p$ix[[n]])
  p$nStates = max(ix)
  return(p)
}
#
# Make a basic three-species setup based upon Petrik et al (2019): Bottom-up drivers of global patterns of demersal, forage, and pelagic fishes. Progress in Oceanography 176, 102124. doi 10.1016/j.pocean.2019.102124.
#
# Out:
#  An updated parameter list. The list contains:
#   ixResources - the indices to the resources
#   ixGroup - array of nGroups with indices for each group
#   mMature(nGroups) - mass of maturation of each group
#
setupBasic = function(depth = 500, pprod = 100) {
  # Initialize the parameters:
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
  #
  # Add fish groups:
  #
  param = paramAddGroup(param, 0.001, 250, 0.5, 2) # Small pelagics
  param = paramAddGroup(param, 0.001, 125000, 250, 3) # Large pelagics
  param = paramAddGroup(param, 0.001, 20000, 250, 3) # Demersals
  #
  # Setup physiology:
  #
  h = 20; # Max. consumption coefficient
  n = -0.25 # Metabolic exponent
  q = -0.2 # Clearance rate exponent
  gamma = 70 # Coef. for clearance rate
  ix = param$ixFish
  m = param$mc[ix]
  param$Cmax[ix] = h*m^n
  param$V[ix] = gamma*m^q
  param$metabolism[ix] = 0.2*param$Cmax[ix]
  param$epsRepro = rep(0.01, param$nGroups)
  param$epsAssim = 0.7 # Assimilation efficiency
  param$Cmax[ is.na(param$Cmax) ] = 0
  param$V[ is.na(param$V) ] = 0
  #
  # Setup interaction matrix:
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
  # Mortality
  #
  param$mort0 = 0.1
  
  return(param)
}

#
# Calculate the derivatives of all state variables
#
# In: 
#  p : a list of parameters
#  u : all state variables
#  bFullOutput : if TRUE returns all internal calculations
#
# Out:
#  deriv : a vector of derivatives
#
calcDerivatives = function(t, u, p, bFullOutput=FALSE) {
  ix = p$ixFish
  B = u[ix]
  #
  # Calc consumption of all fish groups
  #
  
  # Calc Encounter:
  Enc = p$V * (p$theta %*% u)
  f = Enc / (p$Cmax + Enc) # Functional response
  f[is.na(f)] = 0
  Eavail = p$epsAssim * p$Cmax*f - p$metabolism
  
  # Predation mortality:
  mortpred = t(p$theta) %*% (f*p$Cmax/p$epsAssim*u/p$mc)

  # Total mortality
  mort = mortpred + p$mort0
  
  #
  # Derivate of fish groups
  #
  
  # Flux out of the size group:
  v = Eavail[ix]
  vplus = pmax(0,v)
  kappa = 1-p$psiMature[ix]
  gamma = (kappa*vplus - mort[ix]) /
    (1 - (1/p$z[ix])^(1-mort[ix]/(kappa*vplus)) )
  gamma[kappa==0] = 0 # No growth of fully mature classes
  Fout = gamma*B
  Repro = (1-kappa)*vplus*B
  
  # Flux into the size group
  Fin = 0
  for (i in 1:p$nGroups) {
    ix = p$ix[[i]]-p$ix[[1]][1]+1
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
    return(out)
  }
  else
    return( list(c(dRdt, dBdt)) )
}

#
# Simulate the model.
#
# In: 
#  p : fully populated set of parameters
#  tEnd : time to simulate (years)
#
# Out:
#  A simulation list
# 
simulate= function(p = setupBasic(), tEnd = 100) {
  #
  # Integrate the equations:
  #
  t = seq(0, tEnd,length.out=tEnd+1)
  u = ode(y = p$u0, 
            times = t, 
            func = calcDerivatives, 
            parms=p)
  #
  # Assemble output:
  #
  sim = list()
  sim$p = p
  sim$R = u[, p$ixR+1]
  sim$B = u[, p$ixFish+1]
  sim$t = t
  sim$nTime = length(t)
  #
  # Calculate Spawning Stock Biomass
  #
  SSB = matrix(nrow=length(t), ncol=p$nGroups)
  for (i in 1:p$nGroups)
    for (j in 1:length(t))
      SSB[j,i] = sum( u[j, p$ix[[i]]] * p$psiMature[p$ix[[i]]] )
  sim$SSB = SSB
  
  return(sim)
}
