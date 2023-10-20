#===============================================================================
# Functions to creating parameter settings for the Feisty model
#
# Slightly rewritten by Karline Soetaert, based on code from Ken H. Andersen
# 
#===============================================================================

#-------------------------------------------------------------------------------
# A function to estimate feeding preferences based on 
# size differences of prey and predator
#-------------------------------------------------------------------------------

sizePrefFeeding <- function(
    p,           # parameter settings 
    beta = 400,  # preferred predator/prey mass ratio
    sigma = 1.3, # width of size preference for feeding
    type = 1) # 1 = normal, 2=errf
  {
  
  theta = matrix(nrow=p$nStages, ncol=p$nStages, data=0)
  rownames(theta) <- colnames(theta) <- p$stagenames
  if (type == 1){
    for (i in p$ixFish) {
      theta[i,] = exp( -(log(p$mc[i]/(beta*p$mc)))^2 / (2*sigma)^2  )
      theta[i, p$mc >p$mc[i]] = 0
    }
  } else {
  # calculate size-preference matrix based on erf
    for (i in p$ixFish){
       for (j in 1: p$nStages){
            theta[i, j] = sqrt(pi/2)*sigma*(
              erf((log(p$mUpper[j]) - log(p$mc[i]/beta))/(sqrt(2)*sigma))
                  - erf((log(p$mLower[j]) - log(p$mc[i]/beta))/(sqrt(2)*sigma)))
            theta[i, j] = theta[i, j]/(log(p$mUpper[j]) - log(p$mLower[j]))
       }
    }
  }  
  theta[is.na(theta)] = 0
  return(theta)  
}

#-------------------------------------------------------------------------------
# Initialize the parameter structure
#
# In:
#  depth : depth (meters)
#  pprod, bprod : Primary and sediment productivity 
#

# Out:
#  A parameter list
#-------------------------------------------------------------------------------

paramInit = function(...) {
  param = list(...)
  
  # Prepare for groups to be added:
  param$nResources = 0
  param$nGroups = 0
  param$nStages = 0
  param$groupnames=NA
  param$stagenames=NA
  param$resources=NA
  param$fishes=NA
  param$groups=NA
  param$theta=NA
  param$ix     = list() # Indexes of each fish group
  param$ixFish = list() # Indices to all fish groups
  param$ixR = NA # Indices to all resource groups
  
  return(param)
}

#-------------------------------------------------------------------------------
# Makes a size grid for predators.
#
# In:
#  mMin, mMax : minimum of first and maximum sizes of last stage (gram).
#  nStages : number of stages
#
# Out:
#  A list with lower and upper sizes, the ratio between the two (z) and
#  the center mass (mc).
#-------------------------------------------------------------------------------

makeGrid = function(mMin,         # min size, gram
                    mMax,         # max size, gram
                    nStages=1) {  # number of stages
#
# Setup mass grid (log scaled)
#

# minimal sizes; the last being the upper size of the last cell
  m = exp(seq(from=log(mMin), to=log(mMax), length.out=(nStages+1)))
  
  mLower = m[1:nStages]
  mUpper = m[2:(nStages+1)]

# The ratio between upper and lower sizes
  z = mUpper / mLower 

# Geometric mean center mass  
  mc = exp( log(mLower) + 0.5*(log(z)) ) 
  
  return(list(mLower=mLower, mUpper=mUpper, z=z, mc=mc))
}

# ------------------------------------------------------------------------------
# Sets up resources
#
# Input: parameter list
#
# Output: parameter list with resource parameters
#
# ------------------------------------------------------------------------------

paramAddResource = function(p,        # parameter to be updated
                         K,           # maximum resource concentration
                         r=1,         # resource nudging rate (/yr)          
                         dynamics=c("chemostat", "logistic"),
                         mc,          # mean weight, gWW/
                         mLower = NA, # weight lower limit
                         mUpper = NA, #upper limit
                         names=NA,
                         u0=NA){
  nR = length(K)
  p$nResources = nR
  if (any(is.na(names))) names <- paste("Resource",1:nR,sep="")
  p$groupnames = names
  p$stagenames = names
  p$resources = data.frame(K=K, r=r, mc=mc, 
                          mLower=mLower, mUpper=mUpper, u0=u0)
  row.names(p$resources) = names
  p$dynamics  <- match.arg(dynamics)
  p$Rtype <- pmatch(match.arg(dynamics), c("chemostat", "logistic"))
  p$K     = K
  p$r     = r
  p$ixR   = 1:nR 
  p$mc    = mc
  p$mLower = p$resource$mLower
  p$mUpper = p$resource$mUpper
  p$u0     = p$resource$u0
  if (any(is.na(p$u0))) p$u0    = K
  names(p$u0)= names(p$mLower) = names(p$mUpper) = names
  p
}

# ------------------------------------------------------------------------------
# Adds a group to the parameter list. 
#
# input:
# The group is defined by its min and
# max sizes, its size at maturation, and the number of stages.
# default parameters may be overruled.
#
# ------------------------------------------------------------------------------

paramAddGroup = function(p ,           # list of parameters to be updated
                         nStages,      # number of stages
                         mMin,         # minimum mass of first stage, g WW
                         mMax,         # minimum mass of last stage, gWW
                         mMature=mMax, # size at which 50% of individuals is mature
                                       # if NA: only last size class is for 50% mature
                         mortF=0,      # mortality imposed due to fishing
                         mort0=0.1,    # natural mortality 
                         u0=1,         # initial conditions
                         name=NA)      # name of the new group
  {
  p$nGroups = p$nGroups + 1                                                     
  n         = p$nGroups
  
  # update names of groups and stages
  if (any(is.na(name))) 
    name <- paste("Group", n, sep="")
  p$groupnames <- c(p$groupnames, name)
  
  p$stagenames <- c(p$stagenames,
               paste(name, 1:nStages, sep="_")) 

  # Setup mass grid:
  grid = makeGrid(mMin, mMax, nStages)
  
  # update indexes
  if (n==1)    # first group: first index = number of resources+1
    ixStart = max(p$ixR)+1
  else         # 
    ixStart = max(p$ix[[n-1]])+1
  
  ix = ixStart:(ixStart+nStages-1)
  
  p$ix[[n]]    = ix
  p$mLower[ix] = grid$mLower 
  p$mUpper[ix] = grid$mUpper
  p$z[ix]      = grid$z    
  p$mc[ix]     = grid$mc
  
  # Setup maturation schedule: part of fishes in each stage that is mature

  p$mMature[n] = mMature  # half-saturation maturation coefficient
  if (is.na(mMature)){    # Only last class mature for 50% 
    p$psiMature[ix] = 0
    p$psiMature[ix[[nStages]]] = 0.5
  } else p$psiMature[ix] = ( 1 + (p$mc[ix]/mMature)^(-5) )^(-1)

# This is the same, but easier to understand (for me - KS)  
#  p$psiMature[ix] = p$mc[ix]^5/(p$mc[ix]^5+mMature^5)
  
  # fishing mortality:

  p$mortF[ix] = mortF
  p$mort0[ix] = mort0
  
  p$u0[ix]        = u0 # Initial condition
  
  names(p$u0[ix]) = paste(name, 1:nStages, sep="_") # Initial conditiona
  
  # indices to all fish states
  p$ixFish = min(p$ix[[1]]) : max(p$ix[[n]])
  
  # total number of states
  p$nStages = max(ix)
  return(p)
}

# ------------------------------------------------------------------------------
# Sets up physiological parameters of each fish stage
#
# Input: parameter list
#
# Output: parameter list with updated rate parameters
#
# ------------------------------------------------------------------------------

paramAddPhysiology = function (p, 
              ac = 20,          # Max. consumption coefficient  [g^(-n)/yr]
              bc = -0.25,       # Max. consumption exponent     [-]
              am = 0.011*365,#4,           # Metabolism coefficient        [g^(-p)/yr]
              bm = -0.175,      # Metabolism exponent           [-]
              ae = 70,      # Coef. for clearance rate      [m2*g^(q-1)/yr] encounter slope
              be = -0.2,        # Clearance rate exponent
              epsRepro = 0.01, # reproduction * recruitment efficiency )
              epsAssim = 0.7   # Assimilation efficiency  
                              )
{
  # index pointing to fish stages and the fish weights:

  ix = p$ixFish
  m  = p$mc[ix]  # size of fish
  
  # maximum consumption rate, [/yr] 
  p$Cmax[ix]       = ac*m^bc          
  
  # standard metabolism (basal respiration), [/yr]
  p$metabolism[ix] = am*m^bm                   
  
  # clearance rate, search rate [m2/g/yr]  
  p$V[ix]          = ae*m^be                       
  
  # reproduction * recruitment efficiency , [-]
  p$epsRepro = rep(epsRepro, times=p$nGroups)  
  
  # Assimilation efficiency, [-]
  p$epsAssim = epsAssim                   
  
  #Oct 2023 add for Temperature effects
  p$Cmaxsave=p$Cmax
  p$Vsave=p$V
  p$metabolismsave=p$metabolism
  
  # remove the NAs
  # remove the NAs
  p$mc       [is.na(p$mc)]        <- 0 
  p$u0       [is.na(p$u0)]        <- 0 
  p$z        [is.na(p$z)]         <- 0 
  p$psiMature[is.na(p$psiMature)] <- 0 
  p$mortF    [is.na(p$mortF)]     <- 0 
  p$mort0    [is.na(p$mort0)]     <- 0 
  p$Cmax     [is.na(p$Cmax)]      <- 0 
  p$metabolism[is.na(p$metabolism)] <- 0 
  p$V        [is.na(p$V)]         <- 0 

  names(p$u0) =  p$stagenames  
  iF = p$ixFish
  p$fishes = data.frame(mc=p$mc[iF], mLower=p$mLower[iF], mUpper=p$mUpper[iF], 
                        z=p$z[iF], psiMature=p$psiMature[iF], mortF=p$mortF[iF], 
                        mort0=p$mort0[iF], Cmax=p$Cmax[iF], 
                        metabolism=p$metabolism[iF], V=p$V[iF] )
  row.names(p$fishes)= p$stagenames[iF]  
  p$groups = data.frame(epsRepro=p$epsRepro, epsAssim=p$epsAssim)  
  row.names(p$groups) = p$groupnames[-p$ixR]
  return(p)
}

#--------------------------
#
#
#--------------------------
paramTeffect = function (p, # only for setupbasic & 2
                         Tref,
                         Q10,
                         Q10m,
                         u=NA){ # B container: R+fish
  
  p$Q10  = Q10
  p$Q10m = Q10m
  p$Tref = Tref
  # 
  p$fT = Q10^((p$Tp - Tref)/10) # factor of T effects on V Cmax
  p$fT_met = Q10m^((p$Tp - Tref)/10) # factor of T effects on metabolism
  p$fT_dem = Q10^((p$Tb - Tref)/10) # demersal
  p$fT_met_dem = Q10m^((p$Tb - Tref)/10)
  
  if (p$depth<200){
  
  lambda = (sum(u[p$ix[[1]][p$mc[p$ix[[1]]]>10]]) + sum(u[p$ix[[2]][p$mc[p$ix[[2]]]>10 & p$mc[p$ix[[2]]]<5000]]))/
         (sum(u[p$ix[[1]][p$mc[p$ix[[1]]]>10]]) + sum(u[p$ix[[2]][p$mc[p$ix[[2]]]>10 & p$mc[p$ix[[2]]]<5000]]) + 
          sum(u[p$ix[[3]][p$mc[p$ix[[3]]]>10 & p$mc[p$ix[[3]]]<5000]]) + sum(u[3:4])) #Eq. 15
  
  lambda=0.5 #
  p$eT = p$Tp * lambda + p$Tb * (1-lambda) # effect temperature for adult demersals in shallow water (<200m) Eq. 16
  
  p$fT_dem_shallow=Q10^((p$eT - Tref)/10) # demersal
  p$fT_met_dem_shallow=Q10m^((p$eT - Tref)/10)
  }
  
  ix = c(p$ix[[1]],p$ix[[2]])

  #pelagics
  p$V[ix]= p$fT * p$V[ix]
  
  p$Cmax[ix]= p$fT * p$Cmax[ix]
  
  p$metabolism[ix] = p$fT_met * p$metabolism[ix]
  
  #demersal
  #small
  p$V[p$ix[[3]]][p$mc[p$ix[[3]]]<=10] = p$fT * p$Vsave[p$ix[[3]]][p$mc[p$ix[[3]]]<=10] # small demersal are pelagic
  p$Cmax[p$ix[[3]]][p$mc[p$ix[[3]]]<=10] = p$fT * p$Cmaxsave[p$ix[[3]]][p$mc[p$ix[[3]]]<=10]
  p$metabolism[p$ix[[3]]][p$mc[p$ix[[3]]]<=10] = p$fT_met * p$metabolismsave[p$ix[[3]]][p$mc[p$ix[[3]]]<=10]
  #medium
  p$V[p$ix[[3]]][p$mc[p$ix[[3]]]>10 & p$mc[p$ix[[3]]]<5000] = p$fT_dem * p$Vsave[p$ix[[3]]][p$mc[p$ix[[3]]]>10 & p$mc[p$ix[[3]]]<5000] # small demersal are pelagic
  p$Cmax[p$ix[[3]]][p$mc[p$ix[[3]]]>10 & p$mc[p$ix[[3]]]<5000] = p$fT_dem * p$Cmaxsave[p$ix[[3]]][p$mc[p$ix[[3]]]>10 & p$mc[p$ix[[3]]]<5000]
  p$metabolism[p$ix[[3]]][p$mc[p$ix[[3]]]>10 & p$mc[p$ix[[3]]]<5000] = p$fT_met_dem * p$metabolismsave[p$ix[[3]]][p$mc[p$ix[[3]]]>10 & p$mc[p$ix[[3]]]<5000]
    
  if (p$depth < 200){
    p$V[p$ix[[3]]][p$mc[p$ix[[3]]]>5000] = p$fT_dem_shallow * p$Vsave[p$ix[[3]]][p$mc[p$ix[[3]]]>5000] #
    p$Cmax[p$ix[[3]]][p$mc[p$ix[[3]]]>5000] = p$fT_dem_shallow * p$Cmaxsave[p$ix[[3]]][p$mc[p$ix[[3]]]>5000]
    p$metabolism[p$ix[[3]]][p$mc[p$ix[[3]]]>5000] = p$fT_met_dem_shallow * p$metabolismsave[p$ix[[3]]][p$mc[p$ix[[3]]]>5000]
  }else{
    p$V[p$ix[[3]]][p$mc[p$ix[[3]]]>5000] = p$fT_dem * p$Vsave[p$ix[[3]]][p$mc[p$ix[[3]]]>5000] #
    p$Cmax[p$ix[[3]]][p$mc[p$ix[[3]]]>5000] = p$fT_dem * p$Cmaxsave[p$ix[[3]]][p$mc[p$ix[[3]]]>5000]
    p$metabolism[p$ix[[3]]][p$mc[p$ix[[3]]]>5000] = p$fT_met_dem * p$metabolismsave[p$ix[[3]]][p$mc[p$ix[[3]]]>5000]
    
  }
  
  
  return(p)
}
  
                         

#===============================================================================
# MODEL APPLICATIONS 
#===============================================================================

# ------------------------------------------------------------------------------
# Make a basic three-species setup as described in Petrik et al (2019): Bottom-up 
# drivers of global patterns of demersal, forage, and pelagic fishes. Progress 
# in Oceanography 176, 102124. doi 10.1016/j.pocean.2019.102124.
#
# Out:
#  An updated parameter list. The list contains:
#   ixResources - the indices to the resources
#   ixGroup - array of nGroups with indices for each group
#   mMature(nGroups) - mass of maturation of each group
# ------------------------------------------------------------------------------

setupBasic = function(pprod = 100, # primary production?
                      bprod = 5) { # benthic production?
  
  # Initialize the parameters:
  param = paramInit(depth=0, pprod=pprod, bprod=bprod)
  
  # Add resource:
  param = paramAddResource(
        param, 
        names= c("smallZoo", "largeZoo", "smallBenthos", "largeBenthos"),
        K    = c(pprod, pprod, bprod, 0),  # g ww/m2  - maximum resource concentration
        r    = c(1, 1, 1, 1),              # [/yr] nudging coefficient
        mc   = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500)))

  # Add fish groups:
  # mMature=NA overrides the generic psiMature-> only adult classes 50% mature
  param = paramAddGroup(param, mMin=0.001, mMax=   250, mMature=NA, 
                        mortF=c(0,0.3),      nStages=2, name="smallPel")
  
  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=NA, 
                        mortF=c(0,0.03,0.3), nStages=3, name="largePel") 
  
  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=NA, 
                        mortF=c(0,0.03,0.3), nStages=3, name="Demersals")
  
  # physiology of all fish stages
  param = paramAddPhysiology(param)

  #
  # Setup size interaction matrix:
  #

  # preference matrix: columns=predator, rows=prey
  param$theta = matrix(nrow=param$nStages, ncol=param$nStages, data=0)
  rownames(param$theta) <- colnames(param$theta) <- param$stagenames
  
  # Small pelagics:
  
  param$theta["smallPel_1", "smallZoo"] = 1 # Small ones eat only small zooplankton

  param$theta["smallPel_2", "smallZoo"] = 0.25
  param$theta["smallPel_2", "largeZoo"] = 1
  param$theta["smallPel_2", "smallPel_1"] = 1
  param$theta["smallPel_2", "largePel_1"] = 1
  param$theta["smallPel_2", "Demersals_1"] = 1

  # Large pelagics:
  param$theta["largePel_1", "smallZoo"] = 1    
  param$theta["largePel_2", "smallZoo"] = 0.25 
  param$theta["largePel_2", "largeZoo"] = 1   
  param$theta["largePel_2", c("smallPel_1", "largePel_1", "Demersals_1")] = 1 
  param$theta["largePel_3", "smallPel_2"] = 0.5 
  param$theta["largePel_3", "largePel_2"] = 1 

  # Demersals:
  param$theta["Demersals_1", "smallZoo"] = 1
  param$theta["Demersals_2", "smallBenthos"] = 1
  param$theta["Demersals_3", "smallPel_2"] = 1
  param$theta["Demersals_3", "smallPel_2"] = 0.75/2
  param$theta["Demersals_3", "largePel_2"] = 0.75 
  param$theta["Demersals_3", "Demersals_2"] = 1 

  return(param)
}

# ------------------------------------------------------------------------------
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
# ------------------------------------------------------------------------------

setupBasic2 = function(pprod = 100, bprod=5, nStages=9) {
  
  # Initialize the parameters:
  param = paramInit(depth=0, pprod=pprod, bprod=bprod)

  # Setup resource groups:
  param = paramAddResource(
        param, 
        names= c("smallZoo", "largeZoo", "smallBenthos", "largeBenthos"),
        K    = c(pprod, pprod, bprod, 0),  # g ww/m2  - maximum resource concentration
        r    = c(1, 1, 1, 1),              # [/yr] nudging coefficient
        mc   = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500)))


# Add fish groups:
  nSmall = round(0.66*nStages)
  param = paramAddGroup(param, mMin=0.001, mMax=   250, mMature=0.5, 
                        mortF=0,      nStages=nSmall, name="smallPel")
  
  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=250, 
                        mortF=0, nStages=nStages, name="largePel") 
  
  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=250, 
                        mortF=0, nStages=nStages, name="Demersals")

  # physiology of all fish stages
  param = paramAddPhysiology(param)

  # Setup size interaction matrix:
  thetaA = 0.5  # Large fish pref for medium forage fish
  thetaD = 0.75 # Pref of large demersal on pelagic prey
  
  param$theta = matrix(nrow=param$nStages, ncol=param$nStages, data=0)
  rownames(param$theta) <- colnames(param$theta) <- param$stagenames

   # Size-based interactions:  
   beta  = 400
   sigma = 1.3

   for (i in param$ixFish) {
     param$theta[i,] = exp( -(log(param$mc[i]/(beta*param$mc)))^2 / (2*sigma)^2  )
     param$theta[i, param$mc > param$mc[i]] = 0
   }
   param$theta[is.na(param$theta)] = 0
  
   #
   # Setup interactions between groups and resources:
   #
   ixSmall = param$ix[[1]]
   ixLarge = param$ix[[2]]
   ixDem   = param$ix[[3]]
  
   mMedium = 10
   mLarge = 5000
   ixSmallSizeDem = ixDem[ (param$mc[ixDem]<=mMedium) ]
   ixMediumSizeDem = ixDem[ (param$mc[ixDem]>mMedium) &
                            (param$mc[ixDem]<mLarge) ]
   
   # Pelagic/demersal indices:
   ixR = param$ixR

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

# ------------------------------------------------------------------------------
# Make a basic four-species setup based up setupBasic(), but generalised to:
# - distinguish between visual and twilight predators
# - with vertical zooplankton
#
# Out:
#  An updated parameter list. The feeding preferences are quite complex
#
# ------------------------------------------------------------------------------

setupVertical = function(pprod = 80, 
                         depth=1500,      # water depth
                         photic=150,      # photic zone depth
                         mesopelagic=250, # mesopelagic depth
                         bent=150,        # ??
                         visual=1.5  ) {# >1 visual predation primarily during the day, = 1 equal day and night
  
  #------------------  
  # Initialize the parameters:
  # habitat and small benthos
  #------------------  
  bprod=0.1*(bent*(depth/photic)^-0.86)

  param = paramInit(bottom=depth, pprod=pprod, bprod=bprod, photic=photic,
                    mesop=mesopelagic, visual=visual, bent=bent)
  
  #------------------  
  # Setup resource groups:
  #------------------  

    param = paramAddResource(
        param, 
        names= c("smallZoo", "largeZoo", "smallBenthos", "largeBenthos"),
        K    = c(pprod, pprod, bprod, 0),  # g ww/m2  - maximum resource concentration
        r    = c(1, 1, 1, 1),              # [/yr] nudging coefficient
        mc   = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500)),
        mLower = c(2e-06,0.001, 0.5e-03, 0.25), # weight lower limit
        mUpper = c(0.001, 0.5, 125, 125),
        u0     = 0.5)

  #------------------  
  # Add fish groups:
  #------------------  
  # mMature=NA overrides the generic psiMature-> only adult classes 50% mature
  u0  = 0.0001
  param = paramAddGroup(param, mMin=0.001, mMax=   250, mMature=NA, u0=u0,
                        mortF=c(0,0.3),      nStages=2, name="smallPel")
  

  u0M = u0  # initial condition = 0 if no mesopelagic zone
  if (param$bottom > param$mesop) u0m <- 0
  
  param = paramAddGroup(param, mMin=0.001, mMax=   250, mMature=NA, u0=u0M,
                        mortF=c(0,0.5),   nStages=2, name="mesoPel")

  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=NA, u0=u0, 
                        mortF=c(0,0,0.5), nStages=3, name="largePel") 
  
  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=NA, u0=u0, 
                        mortF=c(0,0,0.5), nStages=3, name="bathyPel") 

  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=NA, u0=u0,
                        mortF=c(0,0,0.5), nStages=3, name="Demersals")

  #------------------  
  # Setup physiology:
  #------------------  
  param = paramAddPhysiology(param)
    
  #------------------  
  #------------------  
  # theta (preferences):
  #------------------  
  #------------------  
  param$theta  = matrix(nrow=param$nStages, ncol=param$nStages, data=0)
  rownames(param$theta) <- colnames(param$theta) <- param$stagenames

  beta  = 400
  sigma = 1.3
  param$sizeprefer = matrix(nrow=param$nStages, ncol=param$nStages, data=0)
  param$vertover   = matrix(nrow=param$nStages, ncol=param$nStages, data=0)
  
  # calculate size-preference matrix
  for (i in param$ixFish[[1]]: param$nStages){
       for (j in 1: param$nStages){
            param$sizeprefer[i, j] = sqrt(pi/2)*sigma*(
              erf((log(param$mUpper[j]) - log(param$mc[i]/beta))/(sqrt(2)*sigma))
                  - erf((log(param$mLower[j]) - log(param$mc[i]/beta))/(sqrt(2)*sigma)))
  param$sizeprefer[i, j] = param$sizeprefer[i, j]/(log(param$mUpper[j]) - log(param$mLower[j]))
  }
}

  #------------------  
  # overlap from depth distribution
  #------------------  
  ssigma = 10    # width of initial distribution
  tau    = 10    # increase in width
  
  sigmap = ssigma + tau*log10(param$mc/param$mc[1]) # width for each size class
  xrange = 0 : param$bottom
  param$dvm = param$photic + 500 # 650
  
  if (param$bottom < (param$photic + 500)) 
    param$dvm = param$bottom   # migration to bottom in intermediate habitats
  
  if (param$bottom <= param$mesop) 
    param$dvm = 0              # no migration in shallow habitats
  
  ixjuv   = 2   #minloc(abs(sizes-smat)); from matlab
  ixadult = 3   #minloc(abs(sizes-lmat));

  # a function to generate vertical distributions (a normal distribution)
  VertDist <- function(sigma, xloc){
    xloc = rep(xloc, length.out=length(sigma))
    zp_n = matrix(nrow=length(xrange), ncol=length(sigma), data=0) 
    for (i in 1: length(sigma)){      
     zp_n[,i] = (1/(sqrt(2*pi*sigma[i]^2)))* 
       exp(-(((xrange - xloc[i])^2)/(2*sigma[i]^2)))
    }
    zp_n = zp_n %*% diag(1/colSums(zp_n))
    zp_n  
  }
  
  ## zooplankton : small zoo & large zoo
  # at night: zooplankton is close to surface
  zp_n = VertDist(sigmap[1:2], xloc=0)

  # zooplankton day (half at surface, half at dvm depth
  zp_d = VertDist(sigmap[1:2], xloc=param$dvm)
  zp_d = (zp_n + zp_d)/2
  
  ## benthos small and large (at bottom with width ssigma)
   bent_dn = VertDist(c(ssigma, ssigma), xloc=param$bottom)

  ## small pelagic fish (day + night) always at surface
  ix = param$ix[[1]]
  spel_dn = VertDist(sigmap[ix], xloc=0)
 
  ## meso pelagic night   at surface  
  mpel_n = spel_dn
  
  # meso pelagic day (all at dvm)
  ix = param$ix[[2]]
  mpel_d = VertDist(sigmap[ix], xloc=param$dv)

  ## large pelagic fish night (all at surface)
  ix = param$ix[[3]]
  lpel_n = VertDist(sigmap[ix], xloc=0)
  
  # large pelagic fish day (non-adult at surface   adult at dvm)
  xlocvec = rep(0,length(ix)) 
  xlocvec[ixadult:length(xlocvec)] = param$dvm 
  lpel_d = VertDist(sigmap[ix], xloc=xlocvec)
  lpel_d = (lpel_d + lpel_n)/2
  
  ## bathypelagic night (adults in midwater, others at surface)
  ix = param$ix[[4]]
  xlocvec = rep(0,length(ix)) # initialization
  xlocvec[ixadult:length(xlocvec)] = param$dvm # non-adult at surface   adult at dvm
  bpel_n = VertDist(sigmap[ix], xloc=xlocvec)

  # bathypelagic day (all at dvm)
  bpel_d = VertDist(sigmap[ix], xloc=param$dvm)

  ## demersal fish night
  ix = param$ix[[5]]
  xlocvec = rep(0,length(ix)) # initialization
  xlocvec[ixjuv:length(xlocvec)] = param$bottom # larvae at surface   juvenile and adult at bottom
  dem_n = VertDist(sigmap[ix], xlocvec)

  # demersal fish day; larvae at surface/ juvenile at bottom/ adult and middle
  demmig = param$dvm # ? from matlab
  if ((param$bottom - param$dvm) >= 1200) 
    demmig = param$dvm + (param$bottom-param$dvm-1200)
  if ((param$bottom - param$dvm) >= 1500)
   demmig = param$bottom
   
  dem_d= matrix(nrow=length(xrange), ncol=length(param$ix[[5]]), data=0)

  xlocvec[ixadult:length(xlocvec)] = param$dvm ### or demmig???
  dem_d =  VertDist(sigmap[ix], xlocvec)

  #if shallower than euphotic depth, adult demersals feed across-habitats
  if (param$bottom <= param$photic) {
  dem_d = (dem_d + dem_n)/2
  dem_n = dem_d
  }
  
  # calculate overlap during day
  depthDay = matrix(nrow=length(xrange), ncol=param$nStages, data=0)
  test     = matrix(nrow=length(xrange), ncol=param$nStages, data=0)
  param$dayout = matrix(nrow=param$nStages, ncol=param$nStages, data=0)
  
  depthDay[, 1:2] = zp_d # resources
  depthDay[, 3:4] = bent_dn # resources
  depthDay[, param$ix[[1]]] = spel_dn
  depthDay[, param$ix[[2]]] = mpel_d
  depthDay[, param$ix[[3]]] = lpel_d
  depthDay[, param$ix[[4]]] = bpel_d
  depthDay[, param$ix[[5]]] = dem_d
  
  for (i in 1: param$nStages) {
    for ( j in 1: param$nStages) {
     test[, j] = pmin(depthDay[, i], depthDay[, j])
    }
  param$dayout[, i] = colSums(test)
  }
  
  # calculate overlap during night
  depthNight = matrix(nrow=length(xrange), ncol=param$nStages, data=0)
  # test will be overwritten
  param$nightout = matrix(nrow=param$nStages, ncol=param$nStages, data=0)
  
  depthNight[, 1:2] = zp_n # resources
  depthNight[, 3:4] = bent_dn # resources
  depthNight[, param$ix[[1]]] = spel_dn
  depthNight[, param$ix[[2]]] = mpel_n
  depthNight[, param$ix[[3]]] = lpel_n
  depthNight[, param$ix[[4]]] = bpel_n
  depthNight[, param$ix[[5]]] = dem_n
  
  for (i in 1: param$nStages) {
    for ( j in 1: param$nStages) {
      test[, j] = pmin(depthNight[, i], depthNight[, j])
    }
    param$nightout[, i] = colSums(test)
  }
  
  #------------------  
  # visual ability
  #------------------  

  # visual predatars: good at light, bad in the dark
  visualpred = c(param$ix[[1]], # small palegic 5 6 always at surface
                 param$ix[[3]]) # large pelagic 9 10 11
  param$dayout[visualpred,]   = param$dayout[visualpred,]*param$visual       # predation enhanced during day
  param$nightout[visualpred,] = param$nightout[visualpred,]*(2-param$visual) # predation decreased at night 
  
  # pelagic predators: limited vision in twilight zone during day
  pelpred = param$ix[[3]]                    # large pelagic   9 10 11
  pelpred = pelpred[ixadult:length(pelpred)] # adult large pelagic  11  at dvm during day
  preytwi = c(param$ix[[2]], param$ix[[4]])  # mesopelagic 7 8   bathypelagic 12 13 14
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
  prey1 = (param$ix[[1]][1]+ (ixjuv   - 1)) : param$ix[[1]][length(param$ix[[1]])]
  prey2 = (param$ix[[2]][1]+ (ixjuv   - 1)) : param$ix[[2]][length(param$ix[[2]])]
  idx_predat = c(pred1, pred2, pred3)
  idx_prey   = c(prey1, prey2)
  param$theta[idx_predat,idx_prey] = param$theta[idx_predat,idx_prey]*0.5
  
return(param)  
}


# ------------------------------------------------------------------------------
# Make a basic setup with just pelagic fish, 5 stages
#
# Out:
#  An updated parameter list. The feeding preferences are quite complex
#
# ------------------------------------------------------------------------------

setupPelagicSpecies = function(depth=500, pprod=100, bprod=5, 
                               nStages=6, mInf=125000, names=NA, demersal=TRUE, mort0=0.5) {

  # Initialize the parameters:
  param = paramInit(depth=depth, pprod=pprod, bprod=bprod)

  # Setup resource groups:
  param = paramAddResource(
        param, 
        names= c("smallZoo", "largeZoo", "smallBenthos", "largeBenthos"),
        K    = c(pprod, pprod, bprod, 0),  # g ww/m2  - maximum resource concentration
        r    = c(1, 1, 1, 1),              # [/yr] nudging coefficient
        mc   = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500)),
        mLower = c(2e-06,0.001, 0.5e-03, 0.25), # weight lower limit)  
        mUpper = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500))
  )


  # Add fish groups:
  names = rep(names, length.out=length(mInf))
  for (iGroup in 1:length(mInf))
    param = paramAddGroup(param, 
                          mMin=0.001, mMax=mInf[[iGroup]], 
                          mMature=0.25*mInf[[iGroup]], 
                          nStages=nStages, 
                          name=names[iGroup])

  # physiology of all fish stages
  param = paramAddPhysiology(param)

  # Setup size interaction matrix:
  param$theta = matrix(nrow=param$nStages, ncol=param$nStages, data=0)
  rownames(param$theta) <- colnames(param$theta) <- param$stagenames

  beta = 400   # preferred predator/prey mass ratio
  sigma = 1.3  # width of size preference for feeding
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

  ixR = param$ixR

  if (!demersal)
    param$theta[,ixR[3:4]] = 0  # No demersal feeding

  param$mort0 = mort0 # NOTE: set pretty high to give a stable population
  
  return(param)
}


