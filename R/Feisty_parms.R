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
  } else if (type == 2) {
  # calculate size-preference matrix based on erf used in van Denderen et al., 2020
    for (i in p$ixFish){
       for (j in 1: p$nStages){
            theta[i, j] = sqrt(pi/2)*sigma*(
              erf((log(p$mUpper[j]) - log(p$mc[i]/beta))/(sqrt(2)*sigma))
                  - erf((log(p$mLower[j]) - log(p$mc[i]/beta))/(sqrt(2)*sigma)))
            theta[i, j] = theta[i, j]/(log(p$mUpper[j]) - log(p$mLower[j]))
       }
    }
  } else if (type == 3) {
  #  
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
  
  param$my_palette <- c("smallZoo" = "#FFEE58",
                    "largeZoo" = "#F9A825",
                    "smallBenthos" = "#795548",
                    "largeBenthos" = "#F57C0D",
                    "smallPel" = "#BBDEFB",
                    "mesoPel" = "#9E9E9E",
                    "largePel" = "#2196F3",
                    "bathyPel" =  "#0D47A1",
                    "demersals" =  "#800080")
  
  param$my_names <- c("smallZoo" = "Small zooplankton",
                  "largeZoo" = "Large zooplankton",
                  "smallBenthos" = "Small Benthos",
                  "largeBenthos" = "Large Benthos",
                  "smallPel" = "Small pelagics",
                  "mesoPel" = "Mesopelagics",
                  "largePel" = "Large pelagics",
                  "bathyPel" = "Bathypelagics",
                  "demersals" =  "Demersals")
  
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
  p$mLower = p$resources$mLower
  p$mUpper = p$resources$mUpper
  p$u0     = p$resources$u0
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
  
  lambda = (sum(u[p$ix[[1]][p$mc[p$ix[[1]]]>p$mMedium]]) + sum(u[p$ix[[2]][p$mc[p$ix[[2]]]>p$mMedium & p$mc[p$ix[[2]]]<p$mLarge]]))/
         (sum(u[p$ix[[1]][p$mc[p$ix[[1]]]>p$mMedium]]) + sum(u[p$ix[[2]][p$mc[p$ix[[2]]]>p$mMedium & p$mc[p$ix[[2]]]<p$mLarge]]) + 
          sum(u[p$ix[[3]][p$mc[p$ix[[3]]]>p$mMedium & p$mc[p$ix[[3]]]<p$mLarge]]) + sum(u[3:4])) #Eq. 15
  
  lambda=0.5 #
  p$eT = p$Tp * lambda + p$Tb * (1-lambda) # effect temperature for adult demersals in shallow water (<200m) Eq. 16
  
  p$fT_dem_shallow=Q10^((p$eT - Tref)/10) # demersal
  p$fT_met_dem_shallow=Q10m^((p$eT - Tref)/10)
  }
  
  ix = c(p$ix[[1]],p$ix[[2]])

  #pelagics
  p$V[ix]= p$fT * p$Vsave[ix]
  
  p$Cmax[ix]= p$fT * p$Cmaxsave[ix]
  
  p$metabolism[ix] = p$fT_met * p$metabolismsave[ix]
  
  #demersal
  #small
  p$V[p$ix[[3]]][p$mc[p$ix[[3]]]<=p$mMedium] = p$fT * p$Vsave[p$ix[[3]]][p$mc[p$ix[[3]]]<=p$mMedium] # small demersal are pelagic
  p$Cmax[p$ix[[3]]][p$mc[p$ix[[3]]]<=p$mMedium] = p$fT * p$Cmaxsave[p$ix[[3]]][p$mc[p$ix[[3]]]<=p$mMedium]
  p$metabolism[p$ix[[3]]][p$mc[p$ix[[3]]]<=p$mMedium] = p$fT_met * p$metabolismsave[p$ix[[3]]][p$mc[p$ix[[3]]]<=p$mMedium]
  #medium
  p$V[p$ix[[3]]][p$mc[p$ix[[3]]]>p$mMedium & p$mc[p$ix[[3]]]<p$mLarge] = p$fT_dem * p$Vsave[p$ix[[3]]][p$mc[p$ix[[3]]]>p$mMedium & p$mc[p$ix[[3]]]<p$mLarge] # small demersal are pelagic
  p$Cmax[p$ix[[3]]][p$mc[p$ix[[3]]]>p$mMedium & p$mc[p$ix[[3]]]<p$mLarge] = p$fT_dem * p$Cmaxsave[p$ix[[3]]][p$mc[p$ix[[3]]]>p$mMedium & p$mc[p$ix[[3]]]<p$mLarge]
  p$metabolism[p$ix[[3]]][p$mc[p$ix[[3]]]>p$mMedium & p$mc[p$ix[[3]]]<p$mLarge] = p$fT_met_dem * p$metabolismsave[p$ix[[3]]][p$mc[p$ix[[3]]]>p$mMedium & p$mc[p$ix[[3]]]<p$mLarge]
  #Large  
  if (p$depth < 200){
    p$V[p$ix[[3]]][p$mc[p$ix[[3]]]>=p$mLarge] = p$fT_dem_shallow * p$Vsave[p$ix[[3]]][p$mc[p$ix[[3]]]>=p$mLarge] #
    p$Cmax[p$ix[[3]]][p$mc[p$ix[[3]]]>=p$mLarge] = p$fT_dem_shallow * p$Cmaxsave[p$ix[[3]]][p$mc[p$ix[[3]]]>=p$mLarge]
    p$metabolism[p$ix[[3]]][p$mc[p$ix[[3]]]>=p$mLarge] = p$fT_met_dem_shallow * p$metabolismsave[p$ix[[3]]][p$mc[p$ix[[3]]]>=p$mLarge]
  }else{
    p$V[p$ix[[3]]][p$mc[p$ix[[3]]]>=p$mLarge] = p$fT_dem * p$Vsave[p$ix[[3]]][p$mc[p$ix[[3]]]>=p$mLarge] #
    p$Cmax[p$ix[[3]]][p$mc[p$ix[[3]]]>=p$mLarge] = p$fT_dem * p$Cmaxsave[p$ix[[3]]][p$mc[p$ix[[3]]]>=p$mLarge]
    p$metabolism[p$ix[[3]]][p$mc[p$ix[[3]]]>=p$mLarge] = p$fT_met_dem * p$metabolismsave[p$ix[[3]]][p$mc[p$ix[[3]]]>=p$mLarge]
    
  }
  
  
  return(p)
}

