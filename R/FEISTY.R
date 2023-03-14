#
# Core code for the library
#
library('deSolve')
library('pracma') # needed for erf, linspace
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

setupBasic = function(pprod = 100, bprod=5, temps=10, tempb=8) {
  
  # Initialize the parameters:
  param = parametersInit(0, pprod)
  param$bprod = bprod
  param$temps = temps
  param$tempb = tempb
  param$setup=1
  
  #
  #update temperature
  #
  Tref=10
  Q10=1.88
  Q10m=2.35 #Petrik et al.,2019
  
  fTemp=Q10^((temps - Tref)/10)
  fTempm=Q10m^((temps - Tref)/10)
  fTempdem=Q10^((tempb - Tref)/10)
  fTempmdem=Q10m^((tempb - Tref)/10)
  
  #
  # Setup resource groups:
  #
  param$ixR = 1:4 # 4 resources: Small - large zoo small - large benthos 
  param$r = c(1, 1, 1, 0) # Resource growth rates g ww/m2/yr
  param$K = c(pprod, pprod, bprod, 0)  # g ww/m2
  param$mc = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500)) # weight central size)
 # param$mLower = c(2e-06,0.001, 0.5e-03, 0.25) # weight lower limit)  
  #param$mUpper = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500)) # weight central size
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
  #ix = param$ixFish
  ix = c(param$ix[[1]],param$ix[[2]])
  m = param$mc[ix]
  param$Cmax[ix] = fTemp* h*m^n # maximum consumption rate 
  param$V[ix] = fTemp* gamma*m^q # clearance rate 
  param$metabolism[ix] = fTempm* k*m^p # 0.2*param$Cmax[ix] # standard metabolism
  ix = param$ix[[3]]
  m = param$mc[ix]
  param$Cmax[ix] = fTempdem* h*m^n # maximum consumption rate 
  param$V[ix] = fTempdem* gamma*m^q # clearance rate 
  param$metabolism[ix] = fTempmdem* k*m^p # 0.2*param$Cmax[ix] # standard metabolism 
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
  param$mortF[param$ixFish] = 0.3*c(0,1,0,0.1,1,0,0.1,1) # Fishing only on mature stages

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

setupBasic2 = function(pprod = 100, bprod=5, nSizeGroups=9, temps=10, tempb=8) {
  
  # Initialize the parameters:
  param = parametersInit(0, pprod)
  param$bprod = bprod
  param$nSizeGroups = nSizeGroups
  param$temps = temps
  param$tempb = tempb
  param$setup=2
  
  #
  #update temperature
  #
  Tref=10
  Q10=1.88
  Q10m=2.35 #Petrik et al.,2019
  
  fTemp=Q10^((temps - Tref)/10)
  fTempm=Q10m^((temps - Tref)/10)
  fTempdem=Q10^((tempb - Tref)/10)
  fTempmdem=Q10m^((tempb - Tref)/10)
  
  #
  # Setup resource groups:
  #
  param$ixR = 1:4 # 4 resources: Small - large zoo small - large benthos 
  param$r = c(1, 1, 1, 0) # Resource growth rates g ww/m2/yr
  param$K = c(pprod, pprod, bprod, 0)  # g ww/m2
  param$mc = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500)) # weight central size)
  #param$mLower = c(2e-06,0.001, 0.5e-03, 0.25) # weight lower limit)  
  #param$mUpper = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500)) # weight central size
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
  #ix = param$ixFish
  ix = c(param$ix[[1]],param$ix[[2]])
  m = param$mc[ix]
  param$Cmax[ix] = fTemp* h*m^n # maximum consumption rate 
  param$V[ix] = fTemp* gamma*m^q # clearance rate 
  param$metabolism[ix] = fTempm* k*m^p # 0.2*param$Cmax[ix] # standard metabolism
  ix = param$ix[[3]]
  m = param$mc[ix]
  param$Cmax[ix] = fTempdem* h*m^n # maximum consumption rate 
  param$V[ix] = fTempdem* gamma*m^q # clearance rate 
  param$metabolism[ix] = fTempmdem* k*m^p # 0.2*param$Cmax[ix] # standard metabolism 
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
   param$theta[ixMediumSizeDem, 1:2] = 0 
   param$theta[ixMediumSizeDem, param$ixFish] = 0 
   
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

setupVertical = function(pprod = 80, nSizeGroups=6,region = 4,
                         bottom=1500,photic=150) {
  
  # Initialize the parameters:
  param = parametersInit(0, pprod)
  
  # habitat and small benthos
  param$bottom=bottom # water depth default 1500m
  param$photic=photic # photic zone depth default 150m
  param$mesop = 250 # ? depth
  param$visual = 1.5 # scalar; >1 visual predation primarily during the day, = 1 equal day and night
  param$bent = 150
  bprod=0.1*(param$bent*(param$bottom/param$photic)^-0.86)
  param$bprod = bprod
  param$nSizeGroups = nSizeGroups
  param$region = region
  param$setup = 3
  
  #
  # Setup resource groups:
  #
  param$ixR = 1:4 # 4 resources: Small - large zoo small - large benthos 
  param$r = c(1, 1, 1, 0) # Resource growth rates g ww/m2/yr
  param$K = c(pprod, pprod, bprod, 0)  # g ww/m2 
  param$mc = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500)) # weight central size)
  param$mLower = c(2e-06,0.001, 0.5e-03, 0.25) # weight lower limit
  param$mUpper = c(0.001, 0.5, 125, 125) #upper limit

  #
  # Add fish groups:
  #  paramAddGroup = function(p ,mMin, mMax, mMature, nStages) 
  param = paramAddGroup(param, 0.001, 250, 0.5, round(0.66*nSizeGroups))    # Small pelagics
  param = paramAddGroup(param, 0.001, 250, 0.5, round(0.66*nSizeGroups))    # Mesopelagics
  param = paramAddGroup(param, 0.001, 125000, 250, nSizeGroups) # Large pelagics
  param = paramAddGroup(param, 0.001, 125000, 250, nSizeGroups) # Bathypelagics
  param = paramAddGroup(param, 0.001, 125000, 250, nSizeGroups) # Large demersal
 
  # initial conditions
  param$u0[param$ixR] = c(0.5,0.5,0.5,0)
  param$u0[param$ixFish]= 0.0001*param$u0[param$ixFish]
  
  if (param$bottom <= param$mesop) {
    param$u0(param$ix[[2]][1]:param$ix[[2]][length(param$ix[[2]])])=0 #mesopelagics to zero
    param$u0(param$ix[[4]][1]:param$ix[[4]][length(param$ix[[4]])])=0 #mid-water pred to zero
  }
  
  #
  # Override the generic psiMature and make only adult classes 50% mature
  #
   # param$psiMature = 0*param$psiMature
   # for (iGroup in 1:length(param$ix)){
   # param$psiMature[max(param$ix[[iGroup]])] = 0.5
   # }
  
  param$psiMature = 0*param$psiMature
  #overwrite psiMature    from matlab simple run
  nsize=nSizeGroups+1
  sizes = logspace(log10(0.001), log10(125000), nsize) # mMin=0.001     mMax=1.25d5 predatory fish
  matstageS = which.min(abs(sizes-0.5))
  matstageL = which.min(abs(sizes-250))
  param$psiMature[param$ix[[1]]][matstageS:length(param$ix[[1]])] = 0.5 # fishSmall
  param$psiMature[param$ix[[2]]][matstageS:length(param$ix[[2]])] = 0.5 # fishMeso
  param$psiMature[param$ix[[3]]][matstageL:length(param$ix[[3]])] = 0.5 # fishLarge
  param$psiMature[param$ix[[4]]][matstageL:length(param$ix[[4]])] = 0.5 # fishBathy
  param$psiMature[param$ix[[5]]][matstageL:length(param$ix[[5]])] = 0.5 # fishDemersal
   
  #
  # Setup physiology:
  #
  h = 20; # Max. consumption coefficient
  n = -0.25 # Max. consumption exponent
  k = 0.2*h # 0.011*365 Metabolism coefficient  # maintenance costs, 20% of h
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
  # theta calc:
  #
  #
  beta = 400
  sigma = 1.3
  param$theta = matrix(nrow=param$nStates, ncol=param$nStates, data=0)
  param$sizeprefer = matrix(nrow=param$nStates, ncol=param$nStates, data=0)
  param$vertover = matrix(nrow=param$nStates, ncol=param$nStates, data=0)
  
# calculate size-preference matrix
  for (i in param$ixFish[[1]]: param$nStates){
       for (j in 1: param$nStates){
            param$sizeprefer[i, j] = sqrt(pi/2)*sigma*(
              erf((log(param$mUpper[j]) - log(param$mc[i]/beta))/(sqrt(2)*sigma))
                  - erf((log(param$mLower[j]) - log(param$mc[i]/beta))/(sqrt(2)*sigma)))
  param$sizeprefer[i, j] = param$sizeprefer[i, j]/(log(param$mUpper[j]) - log(param$mLower[j]))
  }
}
# calculate overlap from depth distribution
  ssigma = 10 # width of initial distribution
  tau = 10    # increase in width
  
  sigmap = ssigma + tau*log10(param$mc/param$mc[1]) # width for each size class
  xrange = linspace(0, param$bottom, param$bottom + 1)
  param$dvm = param$photic + 500 # 650
  #  from matlab
  if (param$bottom < (param$photic + 500)) {
    param$dvm = param$bottom   # migration to bottom in intermediate habitats
  }
  if (param$bottom <= param$mesop) {
    param$dvm = 0              # no migration in shallow habitats
  }
                
  # ixjuv = 2     #minloc(abs(sizes-smat)); from matlab
  # ixadult = 3   #minloc(abs(sizes-lmat));
  
  ixjuv = which.min(abs(sizes-0.5))
  ixadult = which.min(abs(sizes-250))
  
  # zooplankton night
  xloc = 0 
  zp_n = matrix(nrow=length(xrange), ncol=2, data=0) # ncol=2: small zoo & large zoo
  for (i in 1: 2){      
   zp_n[,i] = (1/(sqrt(2*pi*sigmap[i]^2)))* 
       exp(-(((xrange - xloc)^2)/(2*sigmap[i]^2)))
  }
  zp_n = zp_n %*% diag(1/colSums(zp_n))
  
  # zooplankton day (half at surface, half at dvm depth
  zp_d=matrix(nrow=length(xrange), ncol=2, data=0)  #small zoo & large zoo      
  xloc = param$dvm
  for (i in 1: 2){
    zp_d[, i] = (1/(sqrt(2*pi*sigmap[i]^2)))* 
        exp(-(((xrange - xloc)^2)/(2*sigmap[i]^2)))
  }
  zp_d = zp_d %*% diag(1/colSums(zp_d))
  zp_d = (zp_n + zp_d)/2
  
 # benthos small and large (at bottom with width sigma)
  bent_dn= matrix(nrow=length(xrange), ncol=2, data=0)  #small bent & large bent   
  xloc = param$bottom
  for (i in 1: 2) {  # small benthos & large benthos
  bent_dn[, i] = (1/(sqrt(2*pi*ssigma^2)))*exp(-((xrange - xloc)^2/(2*ssigma^2)))
  bent_dn[, i] = bent_dn[, i]/sum(bent_dn[, i])
}
  
 # small pelagic fish (day + night) always at surface
  #allocate (ix(ixEnd(1) - ixStart(1) + 1))
  spel_dn= matrix(nrow=length(xrange), ncol=length(param$ix[[1]]), data=0)
  xloc = 0
  ix = param$ix[[1]]
  for ( i in 1: length(ix)) {
   spel_dn[, i] = (1/(sqrt(2*pi*sigmap[ix[i]]^2)))* 
     exp(-((xrange - xloc)^2/(2*sigmap[ix[i]]^2)))
  }
  spel_dn = spel_dn %*% diag(1/colSums(spel_dn))
 
  # meso pelagic night   at surface  
  mpel_n = spel_dn
  
  # meso pelagic day (all at dvm)
  mpel_d = matrix(nrow=length(xrange), ncol=length(param$ix[[2]]), data=0)
  xloc = param$dvm
  ix = param$ix[[2]]
  for (i in 1: length(ix)){
    mpel_d[, i] = (1/(sqrt(2*pi*sigmap[ix[i]]^2)))*
     exp(-((xrange - xloc)^2/(2*sigmap[ix[i]]^2)))
  }
  mpel_d = mpel_d %*% diag(1/colSums(mpel_d))
  
  # large pelagic fish night (all at surface)
  lpel_n= matrix(nrow=length(xrange), ncol=length(param$ix[[3]]), data=0)
  xloc = 0
  ix = param$ix[[3]]
  for (i in 1:length(ix)){
  lpel_n[, i] = (1/(sqrt(2*pi*sigmap[ix[i]]^2)))* 
    exp(-((xrange - xloc)^2/(2*sigmap[ix[i]]^2)))
  }
  lpel_n = lpel_n %*% diag(1/colSums(lpel_n))
  
  #large pelagic fish day (non-adult at surface   adult at dvm)
  lpel_d= matrix(nrow=length(xrange), ncol=length(param$ix[[3]]), data=0)
  xlocvec = rep(0,length(ix)) # initialization #ix same as above
  xlocvec[ixadult:length(xlocvec)] = param$dvm # non-adult at surface   adult at dvm
  for (i in 1: length(ix)){ #ix same as above
  lpel_d[, i] = (1/(sqrt(2*pi*sigmap[ix[i]]^2)))* 
    exp(-((xrange - xlocvec[i])^2/(2*sigmap[ix[i]]^2)))
  }
  lpel_d = lpel_d   %*%  diag(1/colSums(lpel_d))
  lpel_d = (lpel_d + lpel_n)/2
  
  # bathypelagic night (adults in midwater, others at surface)
  bpel_n = matrix(nrow=length(xrange), ncol=length(param$ix[[4]]), data=0)
  ix = param$ix[[4]]
  xlocvec = rep(0,length(ix)) # initialization
  xlocvec[ixadult:length(xlocvec)] = param$dvm # non-adult at surface   adult at dvm
  for (i in 1: length(ix)) {
   bpel_n[,i] = (1/(sqrt(2*pi*sigmap[ix[i]]^2)))* 
    exp(-((xrange - xlocvec[i])^2/(2*sigmap[ix[i]]^2)))
  }
  bpel_n = bpel_n %*% diag(1/colSums(bpel_n))
  
  # bathypelagic day (all at dvm)
  bpel_d = matrix(nrow=length(xrange), ncol=length(param$ix[[4]]), data=0)
  xlocvec = rep(param$dvm,length(ix)) # overwrite all elements by dvm
  for (i in 1:length(ix)) {
  bpel_d[, i] = (1/(sqrt(2*pi*sigmap[ix[i]]^2)))*
    exp(-((xrange - xlocvec[i])^2/(2*sigmap[ix[i]]^2)))
}
  bpel_d = bpel_d %*% diag(1/colSums(bpel_d))
  
  # demersal fish night
  dem_n = matrix(nrow=length(xrange), ncol=length(param$ix[[5]]), data=0)
  ix = param$ix[[5]]
  xlocvec = rep(0,length(ix)) # initialization
  xlocvec[ixjuv:length(xlocvec)] = param$bottom # larvae at surface   juvenile and adult at bottom
  for (i in 1: length(ix)) {
  dem_n[, i] = (1/(sqrt(2*pi*sigmap[ix[i]]^2)))* 
    exp(-((xrange - xlocvec[i])^2/(2*sigmap[ix[i]]^2)))
  }
  dem_n = dem_n %*% diag(1/colSums(dem_n))
  
  #demersal fish day
  demmig = param$dvm # ? from matlab
  if ((param$bottom - param$dvm) >= 1200){
  demmig = param$dvm + (param$bottom-param$dvm-1200)
  }
  if ((param$bottom - param$dvm) >= 1500){
  demmig = param$bottom
  }
  dem_d= matrix(nrow=length(xrange), ncol=length(param$ix[[5]]), data=0)
  xlocvec[ixadult:length(xlocvec)] = param$dvm # larvae at surface/ juvenile at bottom/ adult and middle
  for (i in 1: length(ix)) {
  dem_d[, i] = (1/(sqrt(2*pi*sigmap[ix[i]]^2)))* 
    exp(-((xrange - xlocvec[i])^2/(2*sigmap[ix[i]]^2)))
  }
  dem_d = dem_d  %*%  diag(1/colSums(dem_d))
  # ? from matlab
  #if shallower than euphotic depth, adult demersals feed across-habitats
  if (param$bottom <= param$photic) {
  dem_d = (dem_d + dem_n)/2
  dem_n = dem_d
  }
  
  # calculate overlap during day
  depthDay = matrix(nrow=length(xrange), ncol=param$nStates, data=0)
  test = matrix(nrow=length(xrange), ncol=param$nStates, data=0)
  param$dayout = matrix(nrow=param$nStates, ncol=param$nStates, data=0)
  
  depthDay[, 1:2] = zp_d # resources
  depthDay[, 3:4] = bent_dn # resources
  depthDay[, param$ix[[1]]] = spel_dn
  depthDay[, param$ix[[2]]] = mpel_d
  depthDay[, param$ix[[3]]] = lpel_d
  depthDay[, param$ix[[4]]] = bpel_d
  depthDay[, param$ix[[5]]] = dem_d
  
  for (i in 1: param$nStates) {
    for ( j in 1: param$nStates) {
     test[, j] = pmin(depthDay[, i], depthDay[, j])
    }
  param$dayout[, i] = colSums(test)
  }
  
  # calculate overlap during night
  depthNight = matrix(nrow=length(xrange), ncol=param$nStates, data=0)
  # test will be overwritten
  param$nightout = matrix(nrow=param$nStates, ncol=param$nStates, data=0)
  
  depthNight[, 1:2] = zp_n # resources
  depthNight[, 3:4] = bent_dn # resources
  depthNight[, param$ix[[1]]] = spel_dn
  depthNight[, param$ix[[2]]] = mpel_n
  depthNight[, param$ix[[3]]] = lpel_n
  depthNight[, param$ix[[4]]] = bpel_n
  depthNight[, param$ix[[5]]] = dem_n
  
  for (i in 1: param$nStates) {
    for ( j in 1: param$nStates) {
      test[, j] = pmin(depthNight[, i], depthNight[, j])
    }
    param$nightout[, i] = colSums(test)
  }
  
  # visual ability
  # visual predation is good at light, bad in the dark
  visualpred = c(param$ix[[1]], param$ix[[3]]) # small palegic 5 6 always at surface   large pelagic 9 10 11
  param$dayout[visualpred,] = param$dayout[visualpred,]*param$visual # predation enhance during day
  param$nightout[visualpred,] = param$nightout[visualpred,]*(2-param$visual) # predation decrease at night 
  
  # pelagic predators have limited vision in twilight zone during day
  pelpred = param$ix[[3]] # large pelagic   9 10 11
  pelpred = pelpred[ixadult:length(pelpred)] # adult large pelagic  11  at dvm during day
  preytwi = c(param$ix[[2]], param$ix[[4]]) # mesopelagic 7 8   bathypelagic 12 13 14
  param$dayout[pelpred, preytwi] = param$dayout[pelpred, preytwi]/param$visual*(2 - param$visual)  # /1.5 to restore  then *0.5 
  
  # average overlap during the whole day
  param$vertover = (param$dayout + param$nightout)*0.5
  # calculate combined feeding preference matrix
  param$theta = param$sizeprefer*param$vertover
  
  #  specific revision of feeding preference
  idx_be = param$ixFish[1]: (param$ix[[5]][1] + (ixjuv - 2)) # all pelagic and larval demersals
  param$theta[idx_be, 3:4] = 0 # all pelagic and larval demersals do not eat benthos,
                               # only juvenile & adult demersals eat benthos
  
  # small demersals are less preyed on
  idx_smd = (param$ix[[5]][1] + (ixjuv - 1)): (param$ix[[5]][1] + (ixadult - 2)) #
  param$theta[idx_be, idx_smd] = param$theta[idx_be, idx_smd]*0.25
  
  # juvenile & adult demersals do not eat zooplankton
  param$theta[(param$ix[[5]][1] + (ixjuv - 1)) : param$ix[[5]][length(param$ix[[5]])], 1:2] = 0
  
  # provide benefit to forage and mesopelagic fish (predator avoidance)
  pred1 = (param$ix[[3]][1]+ (ixadult - 1)) : param$ix[[3]][length(param$ix[[3]])]
  pred2 = (param$ix[[4]][1]+ (ixadult - 1)) : param$ix[[4]][length(param$ix[[4]])]
  pred3 = (param$ix[[5]][1]+ (ixadult - 1)) : param$ix[[5]][length(param$ix[[5]])]
  prey1 = (param$ix[[1]][1]+ (ixjuv - 1)) : param$ix[[1]][length(param$ix[[1]])]
  prey2 = (param$ix[[2]][1]+ (ixjuv - 1)) : param$ix[[2]][length(param$ix[[2]])]
  idx_predat = c(pred1, pred2, pred3)
  idx_prey= c(prey1, prey2)
  param$theta[idx_predat,idx_prey] = param$theta[idx_predat,idx_prey]*0.5
  
  #...avlocDay  avlocNight    from matlab
  
  # update temperature
  tempdata=read.table("../input/tempdata.dat", sep=',') #
  tempdata[,5]=10 #
  Q10=1.88
  Q10m=1.88
  
  dist=(depthDay+depthNight)/2
  TQ10 =  Q10^((tempdata[1:(param$bottom+1), (region+1)]-10)/10)
  TQ10m =  Q10m^((tempdata[1:(param$bottom+1), (region+1)]-10)/10)
  
  scTemp_step=matrix(0,nrow=size(dist,1),ncol=size(dist,2))
  scTemp_stepm=matrix(0,nrow=size(dist,1),ncol=size(dist,2))
  for (i in 1:size(dist,2)) {
  scTemp_step[,i] = dist[,i] * TQ10
  scTemp_stepm[,i] = dist[,i] * TQ10m
  }
  
  scTemp = colSums(scTemp_step)
  scTempm = colSums(scTemp_stepm)
  
  param$Cmax = scTemp* param$Cmax # maximum consumption rate 
  param$V= scTemp* param$V # clearance rate 
  param$metabolism = scTempm* param$metabolism
  
  #
  # Mortality
  #
  param$mort0 = 0.1
  param$mortF[param$ixFish] = 0 # reset
  param$mortF[length(param$mortF)] = 0.5 # adult large demersal
  
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
  mm = p$Cmax*p$V/(Enc+p$Cmax)*u # temporarily store
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
simulate= function(p = setupBasic(), tEnd = 100, USEdll=TRUE) {
  #
  # Integrate the equations:
  #
  start_time = Sys.time()
  t = seq(0, tEnd,length.out=tEnd+1) 
  #Calc by R or dll
  if (USEdll) {
    loadFEISTYmodel()
    if (p$setup == 1) {
      dummy = .Fortran("f_setupbasic",
        pprod = as.numeric(p$pprod),
        bprod = as.numeric(p$bprod),
        Ts = as.numeric(p$temps),
        Tb = as.numeric(p$tempb)
      )
    } else if (p$setup == 2) {
      dummy = .Fortran("f_setupbasic2",
        pprod = as.numeric(p$pprod),
        bprod = as.numeric(p$bprod),
        nStages = as.integer(p$nSizeGroups),
        Ts = as.numeric(p$temps),
        Tb = as.numeric(p$tempb)
      )
      
    } else if (p$setup == 3) {
      dummy = .Fortran("f_setupvertical",
        pprod = as.numeric(p$pprod),
        nStages = as.integer(p$nSizeGroups),
        region = as.integer(p$region),
        bottom= as.numeric(p$bottom),
        photic= as.numeric(p$photic)
      )
    }
    
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

#!!! PROBABLY INCORRECT
#!!! PROBABLY INCORRECT
#3stages dt=0.1   6stages 0.001        9stages 0.00001(slow)
# setupbasic2 dt=0.01
# setupvertical 6stages 9stages dt=0.1
simulateeuler= function(p = setupVertical(), tEnd = 100, dt=0.1) {
  
  start_time = Sys.time()
  
  loadFEISTYmodel()
  if (p$setup == 1) {
    dummy = .Fortran("f_setupbasic",
                     pprod = as.numeric(p$pprod),
                     bprod = as.numeric(p$bprod),
                     T = as.numeric(p$temp)
    )
  } else if (p$setup == 2) {
    dummy = .Fortran("f_setupbasic2",
                     pprod = as.numeric(p$pprod),
                     bprod = as.numeric(p$bprod),
                     nStages = as.integer(p$nSizeGroups),
                     T = as.numeric(p$temp)
    )
    
  } else if (p$setup == 3) {
    dummy = .Fortran("f_setupvertical",
                     pprod = as.numeric(p$pprod),
                     nStages = as.integer(p$nSizeGroups),
                     region = as.integer(p$region)
    )
  }
  
  res = .Fortran("f_simulateEuler", 
                    u= as.numeric(p$u0),
                    dudt=as.numeric(p$u0),
                    tEnd=as.numeric(tEnd),
                    dt=as.numeric(dt)
                    )
  
  #
  # Assemble output:
  #
  sim = list()
  #p$USEdll=USEdll
  sim$p = p
  sim$R = res$u[p$ixR]
  sim$B = res$u[p$ixFish]
  #sim$t = t
  #sim$nTime = length(t)
  
  #
  # Calculate Spawning Stock Biomass
  #
  # SSB = matrix(nrow=length(t), ncol=p$nGroups)
  # yield = SSB
  # for (i in 1:p$nGroups)
  #   for (j in 1:length(t)) {
  #     SSB[j,i] = sum( u[j, p$ix[[i]]] * p$psiMature[p$ix[[i]]] )
  #     yield[j,i] = sum( u[j, p$ix[[i]]] * p$mortF[p$ix[[i]]] )
  #   }
  # sim$SSB = SSB
  # sim$yield = yield
  
  end_time = Sys.time()
  sim$tictoc=end_time-start_time
  print(sim$tictoc)
  return(sim)

}

#
# for global offline simulation
# design your own loop for all locations
# put loadFEISTYmodel() outside loop in case repeatedly loading library and crash
#
# define 'param$input' shown below before simulation
# p$szprod: small zooplankton production [gww/m2/year]
# p$lzprod: large zooplankton production [gww/m2/year]
# p$bprod:  benthos production [gww/m2/year]
# p$bottom: seafloor depth [m]
# p$photic: euphotic zone depth 
# p$depb: depth grids (boundary data, not grid center) [m] !max 200 layers
# p$Tprof: Temperature grids (boundary data, not grid center) [Celsius] !max 200 layers
# p$nSizeGroups: fish stages !integer
# 
simulateglobal= function(p = setupvertical(), tEnd = 300) {
  #
  # Integrate the equations:
  #
  #start_time = Sys.time()
  t = seq(0, tEnd,length.out=tEnd+1) 
  #Calc by R or dll
  #loadFEISTYmodel()
  dummy = .Fortran("f_setupverticalglobal",
                   szprod=as.numeric(p$szprod),
                   lzprod=as.numeric(p$lzprod),
                   bprod=as.numeric(p$bprod),
                   bottom=as.numeric(p$bottom),
                   photic=as.numeric(p$photic),
                   dgrid=as.numeric(p$depb),
                   tprof=as.numeric(p$Tprof),
                   nStages=as.integer(p$nSizeGroups))
  
  #dudt = assign("dudt", rep(as.double(0),12), envir = .GlobalEnv) 
  u = ode(y=p$u0,
          times = t,
          func = function(t,y,parms) list(calcDerivativesF(t,y,parms)),# Run by dll
          parms = p)
  
  #
  # Assemble output:
  #
  sim = list()
  #p$USEdll=USEdll
  sim$p = p
  sim$R = u[, p$ixR+1]
  sim$B = u[, p$ixFish+1]
  sim$t = t
  sim$nTime = length(t)
  
  #
  # Calculate Spawning Stock Biomass
  #
  # SSB = matrix(nrow=length(t), ncol=p$nGroups)
  # yield = SSB
  # for (i in 1:p$nGroups)
  #   for (j in 1:length(t)) {
  #     SSB[j,i] = sum( u[j, p$ix[[i]]] * p$psiMature[p$ix[[i]]] )
  #     yield[j,i] = sum( u[j, p$ix[[i]]] * p$mortF[p$ix[[i]]] )
  #   }
  # sim$SSB = SSB
  # sim$yield = yield
  
  #end_time = Sys.time()
  #sim$tictoc=end_time-start_time
  #print(sim$tictoc)
  return(sim)
}



