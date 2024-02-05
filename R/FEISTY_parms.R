#===============================================================================
# Functions to creating parameter settings for the FEISTY model
#
# Slightly rewritten by Karline Soetaert, based on code from Ken H. Andersen
# 
#===============================================================================

#' Size Preference Matrix Calculation
#'
#' This function calculates a size preference matrix of a predator for a prey based on the size.
#' 
#' @usage paramSizepref(p, beta = 400, sigma = 1.3, type = 1)
#' 
#' @param p The parameter list to be updated.
#' @param beta The preferred predator/prey mass ratio. Default 400.
#' @param sigma The width of size preference for feeding. Default 1.3.
#' @param type The type of size preference function. 1 for a normal distribution 
#' based function, 2 for a function based on error function, and 3 for a 
#' more complex function integrating over both predator and prey size groups.
#' Default 1.
#'
#' @details
#' This function computes a size preference matrix \code{theta}.
#' Each value represents the preference of a predator (row) for a prey (column) based on their sizes. \cr
#' Three types of functions can be selected:
#' \itemize{
#' \item Type 1 (normal): Based on a normal distribution around the preferred size ratio.
#' \item Type 2 (error function): Based on the error function as described in van Denderen et al., 2020.
#' \item Type 3 (double integration): Based on Andersen and Visser 2023,
#' this method integrates over both predator and prey size groups for a detailed preference calculation.
#' }
#' 
#' The preference is set to 0 for prey larger than the predator.
#' 
#' This function only needs to be called once in each simulation. It needs to be called after all functional types are added, i.e., after the last call of the function \code{\link{paramAddGroup}}.
#' The returned parameter list can be used for further updates.
#' 
#' @return theta: feeding preference matrix for each predator x to each prey y.
#' 
#' @examples
#' # Just an example, data may not make sense.
#' p <- paramInit()
#' # add three resources
#' p <- paramAddResource(p,
#'      names= c("smallZoo", "largeZoo", "smallBenthos"),
#'      K    = c(100, 120, 80),
#'      r    = c(1, 1, 1),
#'      mLower = c(2e-06,0.001, 0.5e-03),
#'      mUpper = c(0.001, 0.5, 125),
#'      mc   = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000)))
#' # add the first functional type
#' p <- paramAddGroup(p, nStages = 3, mMin = 0.1, mMax = 100, mMature = 100*0.25, mortF=0, mort0 = 0.1, name = "smallfish")
#' # add the second functional type
#' p <- paramAddGroup(p, nStages = 6, mMin = 0.1, mMax = 100000, mMature = 100000*0.25, mortF=0, mort0 = 0.1, name = "largefish")
#' # add the third functional type
#' p <- paramAddGroup(p, nStages = 9, mMin = 0.1, mMax = 500000, mMature = 500000*0.25, mortF=0, mort0 = 0.1, name = "giantfish")
#' # add physiological parameters for three functional types (smallfish, largefish and giantfish)
#' p <- paramAddPhysiology(p, 
#'      ac = 20, bc = -0.25,       
#'      am = 0.011*365, bm = -0.175,      
#'      ae = 70, be = -0.2,        
#'      epsRepro = 0.01, 
#'      epsAssim = 0.7)
#' # Calculate size preference matrix
#' p$theta <- paramSizepref(p=p,
#'            beta = 400,
#'            sigma = 1.3,
#'            type = 3)
#' # Run the simulation for this customized FEISTY setup.           
#'  sim=simulateFEISTY(bCust=TRUE,p=p, tEnd=100)
#' # Plot dynamics of resources and small fish.
#'  plot(sim, sim, which=1:6, lty=1, mfrow=c(2,3), subset=time >= 0)       
#'
#' @author Yixin Zhao
#' 
#' @references 
#' van Denderen, P. D., Petrik, C. M., Stock, C. A., & Andersen, K. H. (2021). Emergent global biogeography of marine fish food webs. Global Ecology and Biogeography, 30(9), 1822-1834.
#' 
#' Andersen, K. H., & Visser, A. W. (2023). From cell size and first principles to structure and function of unicellular plankton communities. Progress in Oceanography, 213, 102995.
#' 
#' @aliases paramSizepref
#'
#' @export
#' 

#-------------------------------------------------------------------------------
# A function to estimate feeding preferences based on 
# size differences of prey and predator
#-------------------------------------------------------------------------------

paramSizepref <- function(
    p,           # parameter settings 
    beta = 400,  # preferred predator/prey mass ratio
    sigma = 1.3, # width of size preference for feeding
    type = 1) # 1 = normal, 2=errf
  {
  
  theta = matrix(nrow=p$nStages, ncol=p$nStages, data=0)
  rownames(theta) <- colnames(theta) <- p$stagenames
  if (type == 1) {
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
  # Calculate size preference function by integrating over both predator
  # and prey size groups. See Andersen and Visser 2023, appendix A.  
    s = 2*sigma*sigma
    
    for (i in p$ixFish[1]:p$nStages) {
      for (j in 1:p$nStages){
        z=p$mc[i]/p$mc[j]
        Delta=p$mUpper[i]/p$mLower[i]
        
        theta[i,j] = max (0, (sqrt(Delta)*(((exp(-log((beta*Delta)/z)**2/s) - 2/exp(log(z/beta)**2/s) + 
                                  exp(-log((Delta*z)/beta)**2/s))*s)/2. - 
                                (sqrt(pi)*sqrt(s)*(erf((-log(beta*Delta) + log(z))/sqrt(s))*log((beta*Delta)/z) + 
                                   2*erf(log(z/beta)/sqrt(s))*log(z/beta) + 
                                   erf((log(beta) - log(Delta*z))/sqrt(s))*log((Delta*z)/beta)))/2.))/ ((-1 + Delta)*log(Delta)) 
        )
        
        theta[i,j] = ifelse(p$mc[j]>p$mc[i],0,theta[i,j])
      }
    }
  }
  
  theta[is.na(theta)] = 0
  return(theta)  
}

#' Initialize Parameters for FEISTY
#'
#' This function initializes a list of parameters for the FEISTY model.
#'
#' @usage paramInit(...)
#'
#' @param ... Additional parameters needed to be added to the list, which will be included in the returned parameter list.
#'
#' @return A parameter list is initialized with some basic blank values, vectors, matrices, and sublists.
#' \itemize{
#' \item ..., as input parameters in \code{paramInit(...)}
#' \item nResources, the total number of resources, will be updated later
#' \item nGroups, the total number of fish functional groups, will be updated later
#' \item nStages, the total number of resources + fish stages, will be updated later
#' \item groupnames, a character vector of names of resources and each fish functional group, will be updated later
#' \item stagenames, a character vector of names of resources and each fish stage (functional group), will be updated later
#' \item theta, a matrix for the feeding preference of each predator x to each prey y, will be updated later
#' \item ix, index lists for size classes of each fish functional type (sublists), will be updated later
#' \item ixFish, a vector containing indices for all fish, will be updated later
#' \item ixR, a vector containing indices for all resources, will be updated later
#' \item my_palette, default color palettes for all functional types
#' \item my_names, default full name and acronyms of all functional types
#' }
#'
#' @details This function prepares some basic values, vectors, matrices, and sublists for the following parameter setting up.
#' Default color palettes and names of all resources and functional types (all stages) are defined, which will be used in visualization.
#' The returned parameter list will be used for further updates.
#'
#' @examples
#' # Initialize parameters
#' p <- paramInit(size_1 = 500, size_2 = 1000, depth_1 = 200, depth_2 = 250) # add four extra parameters to the list
#' 
#' @author Ken H. Andersen, Karline Soetaert <karline.soetaert@nioz.nl>, Yixin Zhao
#' 
#' @aliases paramInit
#'
#' @seealso 
#' \code{\link{paramAddResource}} Add Resource Parameters
#'
#' @export

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

#' Add Resource Parameters
#'
#' This function updates the parameter list by adding resource-related parameters.
#' 
#' @usage paramAddResource (p, K, r=1, dynamics=c("chemostat", "logistic"), mc, mLower = NA, mUpper = NA, names=NA, u0=NA)
#'
#' @param p Parameter list to be updated.
#' @param K A vector of carrying capacities of all resources [gww/m2].
#' @param r A vector containing each resource nudging rate (growth rate) [1/year], default 1, does not recommend other values.
#' @param dynamics The type of resource dynamics, either "chemostat" or "logistic".
#' @param mc A vector containing each resource geometric mean weight [gww].
#' @param mLower A vector containing the lower limit of each resource weight [gww]. Optional, depending on the size-based preference calculation function.
#' @param mUpper A vector containing the upper limit of each resource weight [gww]. Optional, depending on the size-based preference calculation function.
#' @param names A character vector of each resource name (acronym). Optional, if not provided, default names are assigned, e.g., Resource_1 and Resource_2.
#' @param u0 A vector of the initial concentration of each resource. If not provided, defaults to the value of \code{K}.
#'
#' @return The updated parameter list \code{p}:
#' \itemize{
#' \item nResources, the total number of resources, which is the length of the input \code{K}.
#' \item groupnames, a character vector of names of resources and each fish functional group, will be updated later.
#' \item stagenames, a character vector of names of resources and each fish stage (functional group), will be updated later.
#' \item dynamics, from parameter input.
#' \item Rtype, resource growth strategy number. Default 1: chemostat, 2: logistic growth. If 'dynamics' is not specified, Rtype=1.
#' \item K, from parameter input.
#' \item r, from parameter input.
#' \item ixR, indices for all resources, start from 1.
#' \item mc, resource geometric mean weight, and fish data will be added later.
#' \item mLower, the lower limit of each resource weight, fish data will be added later.
#' \item mUpper, the upper limit of each resource weight, fish data will be added later.
#' \item u0, from parameter input or same as \code{K}, fish data will be added later.
#' }
#' 
#' @details
#' This function is designed to add parameters of resources to the parameter list.
#' This function only can be called once in each simulation. All resources are added to the parameter list at one time.
#' Generally, this function needs to be called after the function \code{\link{paramInit}}.
#' The returned parameter list can be used for further updates.
#'
#' @examples
#' # Initialize a paramter list
#' p = paramInit()
#' # Add five resources. Just an example, data may not make sense.
#' p = paramAddResource (p=p, K=c(50,100,150,200,250), r=c(1,1,1,1,1), dynamics="chemostat",
#'    mc= c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500), 1e-07*sqrt(500)),
#'    mLower = c(2e-06,0.001, 0.5e-03, 0.25, 1e-07),
#'    mUpper = c(0.001, 0.5, 125, 125, 5e-05),
#'    names=c("smallZoo", "largeZoo", "smallBenthos", "largeBenthos", "Phytoplankton"), u0=c(100,100,200,200,200))
#'
#' @author Ken H. Andersen, Karline Soetaert <karline.soetaert@nioz.nl>, Yixin Zhao
#'
#' @aliases paramAddResource
#' 
#' @seealso 
#' \code{\link{paramAddGroup}} Add Parameters of One Functional Type
#'
#' @export
#' 
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
                         dynamics = "chemostat", # either "chemostat" or "logistic",
                         mc,          # geometric mean weight, gWW
                         mLower = NA, # weight lower limit
                         mUpper = NA, # upper limit
                         names=NA,    # resource name vector, character
                         u0=NA,
                         ixpelR=NA,  # pelagic resource indices
                         ixbenR=NA){ # benthic resource indices
  nR = length(K)
  p$nResources = nR
  if (any(is.na(names))) names <- paste("Resource",1:nR,sep="")
  p$groupnames = names
  p$stagenames = names
  p$resources = data.frame(K=K, r=r, mc=mc, 
                          mLower=mLower, mUpper=mUpper, u0=u0)
  row.names(p$resources) = names
  p$dynamics  <- match.arg(dynamics, c("chemostat", "logistic"))
  p$Rtype <- pmatch(p$dynamics, c("chemostat", "logistic"))
  p$K     = K
  p$r     = r
  p$ixR   = 1:nR 
  p$mc    = mc
  p$mLower = p$resources$mLower
  p$mUpper = p$resources$mUpper
  p$u0     = p$resources$u0
  if (any(is.na(p$u0))) p$u0    = K
  names(p$u0)= names(p$mLower) = names(p$mUpper) = names
  
  if (is.na(ixpelR)) ixpelR <- grep("zoo", names, ignore.case = TRUE)
  p$ixpelR=ixpelR
  if (any(is.na(ixpelR)) || isempty(ixpelR)) stop("Cannot automaticlly get pelagic resource indices. Please check or assign the argument 'ixpelR'.")
  
  if (is.na(ixbenR)) ixbenR <- grep("ben", names, ignore.case = TRUE)
  p$ixbenR=ixbenR
  if (any(is.na(ixbenR)) || isempty(ixbenR)) stop("Cannot automaticlly get benthic resource indices. Please check or assign the argument 'ixbenR'.")
  
  return(p)
}

#' Add Parameters of One Functional Type
#'
#' This function updates the parameter list by adding parameters of one functional type.
#'
#' @usage paramAddGroup(p, nStages, mMin, mMax, mMature=NA, mortF=0, mort0=0.1, u0=1, name=NA)
#'
#' @param p The parameter list to be updated.
#' @param nStages Number of stages for this functional type.
#' @param mMin Minimum size (mass) of this functional type (boundary value), in gram wet weight [gWW].
#' @param mMax Maximum size (mass) of this functional type (boundary value), in gram wet weight [gWW].
#' @param mMature Size [gWW] at which fish has a 50\% maturity level, which will be used for maturity level calculation based on an S-shape function.\cr
#' If NA, only the last size class is 50\% mature; others are 0.
#' @param mortF Fishing mortality for all size classes [1/year]. Default value 0, indicating no fishing.
#' @param mort0 Natural mortality (background mortality) for all size groups [1/year]. Default value 0.1.
#' @param u0 Initial biomass value of all size classes [gWW/m2]. Default value 1.
#' @param name The name (acronym) of the functional type. If not provided, a default name is assigned.
#'
#' @details
#' This function is designed to add parameters of one functional type to the parameter list.
#' Generally, this function needs to be called after the function \code{\link{paramAddResource}}.
#' Every call of this function can add \bold{one} functional group, which means this function needs to be called multiple times to add more functional types.
#'
#' @return The updated parameter list \code{p}:
#' \itemize{
#' \item nGroups, number of all functional types
#' \item groupnames, a character vector of names of resources and each fish functional groups
#' \item stagenames, a character vector of names of resources and each fish stage (functional groups)
#' \item ix, a list including sublists which contains indices of the size classes of each functional type.
#' \item mLower, a vector containing the lower limit of size of each size class of resource and functional type
#' \item mUpper, a vector containing the upper limit of size of each size class of resource and functional type
#' \item z, a vector containing the ratio between \code{mUpper} and \code{mLower}, all resources (NA) + all size classes
#' \item mc, a vector containing the geometric mean size of resources and all stages of functional types
#' \item mMature, a vector containing the size with a 50\% maturity level of each functional type, from parameter input
#' \item psiMature, a vector containing the maturity level, all resources (NA) + all size classes
#' \item mortF, a vector containing fishing mortality, all resources (NA) + all size classes
#' \item mort0, a vector containing all background mortality, all resources (NA) + all size classes
#' \item u0, a vector containing the initial value of biomass, all resources (from \code{\link{paramAddResource}}) + all size classes
#' \item ixFish, a vector of indices of all functional type size classes, starting from \code{nResources+1}
#' \item nStages, an integer, total number of resources and size classes.
#' }
#' All the parameters above will be updated if a new functional type is added by another call of \code{paramAddGroup}. \cr
#' NAs in \code{z}, \code{psiMature}, \code{mortF}, and \code{mort0} will be revised to 0 in \code{\link{paramAddPhysiology}}.
#' 
#' @examples
#' # Just an example, data may not make sense.
#' p <- paramInit()
#' # add four resources
#' p <- paramAddResource(p,
#'      names= c("smallZoo", "largeZoo", "smallBenthos", "largeBenthos"),
#'      K    = c(100, 120, 80, 0),
#'      r    = c(1, 1, 1, 1),
#'      mLower = c(2e-06,0.001, 0.5e-03, 0.25),
#'      mUpper = c(0.001, 0.5, 125, 125),
#'      mc   = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500)))
#' # add the first functional type
#' p <- paramAddGroup(p, nStages = 3, mMin = 0.1, mMax = 100, mMature = 100*0.25, mortF=0, mort0 = 0.1, name = "smallfish")
#' # add the second functional type
#' p <- paramAddGroup(p, nStages = 6, mMin = 0.1, mMax = 100000, mMature = 100000*0.25, mortF=0, mort0 = 0.1, name = "largefish")
#'
#' @author Ken H. Andersen, Karline Soetaert <karline.soetaert@nioz.nl>, Yixin Zhao
#'
#' @aliases paramAddGroup
#' 
#' @seealso 
#' \code{\link{paramAddPhysiology}} Add Physiological Parameters
#' 
#' \code{\link{paramSizepref}} Size Preference Matrix Calculation
#' 
#' @export
#'

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
                         mMature=NA, # size at which 50% of individuals is mature
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
  
  names(p$u0[ix]) = paste(name, 1:nStages, sep="_") # Initial condition
  
  # indices to all fish states
  p$ixFish = min(p$ix[[1]]) : max(p$ix[[n]])
  
  # total number of states
  p$nStages = max(ix)
  return(p)
}

#' Add Physiological Parameters
#' 
#' This function updates the parameter list by adding physiological parameters of all size classes of all functional types.
#'
#' @usage paramAddPhysiology(p, ac = 20, bc = -0.25, am = 0.011*365, bm = -0.175, ae = 70, be = -0.2, epsRepro = 0.01, epsAssim = 0.7)
#'
#' @param p  The parameter to be updated.
#' @param ac Maximum consumption coefficient [g^bc year-1)].
#' @param bc Maximum consumption exponent [-].
#' @param am Metabolism coefficient [g^bm year-1].
#' @param bm Metabolism exponent [-].
#' @param ae Clearance rate (encounter rate) coefficient [m2 g^be year-1].
#' @param be Clearance rate (encounter rate) exponent [-].
#' @param epsRepro Reproduction efficiency [-].
#' @param epsAssim Assimilation efficiency [-].
#'
#' @return The updated parameter list \code{p}:
#' \itemize{
#' \item Cmax, a vector containing maximum consumption rate values, resources (0) + all size classes of all functional types
#' \item metabolism, a vector containing standard metabolism values, resources (0) + all size classes of all functional types
#' \item V, a vector containing clearance rate values, resources (0) + all size classes of all functional types
#' \item epsRepro, a vector of reproduction efficiency values of each functional type (length is the number of functional types). This value applies to all size classes.
#' \item epsAssim, the assimilation efficiency value. This value applies to all size classes.
#' \item Cmaxsave, same as initial values of \code{Cmax}, saved for temperature effect calculation.
#' \item Vsave, same as initial values of \code{V}, saved for temperature effect calculation.
#' \item metabolismsave, same as initial values of \code{metabolism}, saved for temperature effect calculation.
#' }
#'
#' @details 
#' This function is designed to add physiological parameters of all size classes of all functional types to the parameter list.
#' It sets the maximum consumption rate, standard metabolism, and clearance rate based on the size.
#' It ensures there are no NA values in critical parameters. \cr
#' This function only needs to be called once in each simulation. The physiological parameters of all size classes are added to the parameter list at one time.
#' Generally, this function needs to be called after all functional types are added, i.e., after the last call of the function \code{\link{paramAddGroup}}.
#' The returned parameter list can be used for further updates.
#'
#' @examples
#' # Just an example, data may not make sense.
#' # Initialize the parameter list
#' p <- paramInit()
#' # add three resources
#' p <- paramAddResource(p,
#'      names= c("smallZoo", "largeZoo", "smallBenthos"),
#'      K    = c(100, 120, 80),
#'      r    = c(1, 1, 1),
#'      mLower = c(2e-06,0.001, 0.5e-03),
#'      mUpper = c(0.001, 0.5, 125),
#'      mc   = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000)))
#' # add the first functional type
#' p <- paramAddGroup(p, nStages = 3, mMin = 0.1, mMax = 100, mMature = 100*0.25, mortF=0, mort0 = 0.1, name = "smallfish")
#' # add the second functional type
#' p <- paramAddGroup(p, nStages = 6, mMin = 0.1, mMax = 100000, mMature = 100000*0.25, mortF=0, mort0 = 0.1, name = "largefish")
#' # add physiological parameters for two functional types (smallfish and largefish)
#' p <- paramAddPhysiology(p, 
#'      ac = 20, bc = -0.25,       
#'      am = 0.011*365, bm = -0.175,      
#'      ae = 70, be = -0.2,        
#'      epsRepro = 0.01, 
#'      epsAssim = 0.7)
#'
#' @author Ken H. Andersen, Karline Soetaert <karline.soetaert@nioz.nl>, Yixin Zhao
#'
#' @aliases paramAddPhysiology
#' 
#' @seealso 
#' \code{\link{paramTeffect}} Add temperature effects
#'
#' @export

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


#' Add temperature effects
#'
#' This function adjusts temperature-dependent physiological rates based on environmental temperatures. 
#' It applies Q10 temperature coefficients to modify metabolic rates, clearance rates, and maximum consumption rates.
#' Currently, only works on \code{setupBasic}, \code{setupBasic2}, \code{setupVertical}, and \code{setupVertical2}.
#' 
#' @usage paramTeffect(p, Tref, Q10, Q10m, u=NA)
#' 
#' @param p A parameter list. Must be used after \code{\link{paramAddPhysiology}}. \cr
#' `Tp` and `Tb` must be included. See what are `Tp` and `Tb` in \code{\link{setupBasic}} and \code{\link{setupBasic2}}.
#' @param Tref Reference temperature. Default 10 Celsius, generally cannot be other values, unless users define physiological rates based on another reference temperature.
#' @param Q10 Q10 factor for the maximum consumption rate \code{Cmax} and clearance rate \code{V} [-]
#' @param Q10m Q10 factor for metabolism rates \code{metabolism} [-]
#' @param u Temporarily not used. To be developed...
#'
#' @return An updated parameter list:
#' \itemize{
#' \item Tref, from parameter input
#' \item Q10, from parameter input
#' \item Q10m, from parameter input
#' \item fT: factor of T effects on \code{V} and \code{Cmax} of pelagic fish.
#' \item fT_met: factor of T effects on \code{metabolism}.
#' \item fT_dem: factor of T effects on \code{V} and \code{Cmax} of all demersal fish in deep water (depth >= 200m), and small and medium demersal fish in shallow water (depth < 200m).
#' \item fT_met_dem: factor of T effects on \code{metabolism} of all demersal fish in deep water (depth >= 200m), and small and medium demersal fish in shallow water (depth < 200m).
#' \item eT: effective temperature in shallow water (depth < 200m) for large demersal fish. \code{fT_dem_shallow} and \code{fT_met_dem_shallow} are calculated based on \code{eT}.
#' \item fT_dem_shallow: factor of T effects on \code{V} and \code{Cmax} of large demersal fish in shallow water (depth < 200m).
#' \item fT_met_dem_shallow: factor of T effects on the \code{metabolism} of large demersal fish in shallow water (depth < 200m).
#' }
#'
#' @details
#' The function calculates temperature effects on physiological rates (metabolism, clearance rate, and maximum consumption rate) using 
#' the Q10 theory, which states that a physiological process rate increases by a factor of Q10 for every 10 Celsius rise in temperature.
#' The effects are applied differently in terms of living habitats (pelagic and demersal) and water depth (shallow and deep).
#'
#' The function requires a specific structure in the input parameter list `p`, such as
#' indices of all fish and vectors for saved baseline values of physiological rates (\code{Vsave}, \code{Cmaxsave}, \code{metabolismsave}).
#' Therefore it only can be called after \code{\link{paramAddPhysiology}}.
#' `Tp` and `Tb` are required in the input parameter list, or this function will not work.
#' It might not work in customized setups other than \code{\link{setupBasic}} and \code{\link{setupBasic2}}.
#'
#' @examples
#' # Just an example, data may not make sense.
#' # Create a FEISTY setup by setupBasic2.
#' p1 = setupBasic2(szprod = 100, lzprod = 100, bprod = 5, depth = 100, Tp = 10, Tb = 8, 
#' nStages=9, etaMature=0.25, F=0, etaF=0.05)
#' 
#' # change environmental temperatures and update T effects on physiological rates
#' p1$Tp=15
#' p1$Tb=12
#' p2=paramTeffect(p1, p1$Tref, p1$Q10, p1$Q10m, u=NA)
#' # Compare the original and updated physiological rates 
#' p1$V
#' p2$V
#' p1$metabolism
#' p2$metabolism
#' p1$Cmax
#' p2$Cmax
#' 
#' @author Yixin Zhao
#' 
#' @references
#' Petrik, C. M., Stock, C. A., Andersen, K. H., van Denderen, P. D., & Watson, J. R. (2019). Bottom-up drivers of global patterns of demersal, forage, and pelagic fishes. Progress in oceanography, 176, 102124.
#' 
#' van Denderen, P. D., Petrik, C. M., Stock, C. A., & Andersen, K. H. (2021). Emergent global biogeography of marine fish food webs. Global Ecology and Biogeography, 30(9), 1822-1834.
#' 
#' @aliases paramTeffect
#' 
#' @seealso 
#' \code{\link{setupBasic}} The setup following Petrik et al. (2019) \cr
#' \code{\link{setupBasic2}} A revised setup based on `setupBasic`
#' 
#' @export
#' 

paramTeffect = function (p, # only for setupbasic & 2
                         Tref,
                         Q10,
                         Q10m,
                         pelgroupidx, # pelagic group idx, NA means no pelagic groups
                         demgroupidx){# demersal group idx, NA means no demersal groups
  
  p$Q10  = Q10
  p$Q10m = Q10m
  p$Tref = Tref
  p$pelgroupidx=pelgroupidx
  p$demgroupidx=demgroupidx
  # 
  p$fT = Q10^((p$Tp - Tref)/10) # factor of T effects on V Cmax
  p$fT_met = Q10m^((p$Tp - Tref)/10) # factor of T effects on metabolism
  p$fT_dem = Q10^((p$Tb - Tref)/10) # demersal
  p$fT_met_dem = Q10m^((p$Tb - Tref)/10)
  
  #pelagics
  if (!any(is.na(pelgroupidx))) {
    ix = unlist(lapply(pelgroupidx, function(i) p$ix[[i]]))
    
    p$V[ix] = p$fT * p$Vsave[ix]
    p$Cmax[ix] = p$fT * p$Cmaxsave[ix]
    p$metabolism[ix] = p$fT_met * p$metabolismsave[ix]
  }
  
  #demersal
  if (!any(is.na(demgroupidx))) {
    ix = unlist(lapply(demgroupidx, function(i) p$ix[[i]]))
    
    ix_small  = ix[which(p$mc[ix] <= p$mMedium)]
    ix_medium = ix[which(p$mc[ix] > p$mMedium & p$mc[ix] < p$mLarge)]
    ix_large  = ix[which(p$mc[ix] >= p$mLarge)]
    
    p$smdemidx = ix_small
    p$lgdemidx = ix_large
  
    #small
    p$V[ix_small] = p$fT * p$Vsave[ix_small] # small demersal are pelagic
    p$Cmax[ix_small] = p$fT * p$Cmaxsave[ix_small]
    p$metabolism[ix_small] = p$fT_met * p$metabolismsave[ix_small]
    #medium
    p$V[ix_medium] = p$fT_dem * p$Vsave[ix_medium] # medium demersal stay at bottom
    p$Cmax[ix_medium] = p$fT_dem * p$Cmaxsave[ix_medium]
    p$metabolism[ix_medium] = p$fT_met_dem * p$metabolismsave[ix_medium]
    #Large  
    if (p$depth < 200) {
      # effective temperature for large demersals in shallow water (<200m) Eq. 16
      # default no effective T eT=(Tp+Tb)*0.5
      lambda = 0.5 #
      eT = p$Tp * lambda + p$Tb * (1 - lambda)
      
      p$fT_dem_shallow = Q10 ^ ((eT - Tref) / 10) # demersal
      p$fT_met_dem_shallow = Q10m ^ ((eT - Tref) / 10)
      
      p$V[ix_large] = p$fT_dem_shallow * p$Vsave[ix_large] #
      p$Cmax[ix_large] = p$fT_dem_shallow * p$Cmaxsave[ix_large]
      p$metabolism[ix_large] = p$fT_met_dem_shallow * p$metabolismsave[ix_large]
    } else{
      p$V[ix_large] = p$fT_dem * p$Vsave[ix_large] #
      p$Cmax[ix_large] = p$fT_dem * p$Cmaxsave[ix_large]
      p$metabolism[ix_large] = p$fT_met_dem * p$metabolismsave[ix_large]
    }
    
  }
  
  p$pelgrididx = c(p$ixpelR, unlist(lapply(p$pelgroupidx, function(i) p$ix[[i]])), p$smdemidx) # zoop + pelagic fish + small dem
  p$allgrididx = c(p$ixR, p$ixFish) #(1:ncol(p$theta))
  
  return(p)
}

updateET = function (p, # 
                     u=NA){ # B container: R+fish

  # depth < 200m
  
  pelgrididx = c(p$ixpelR, unlist(lapply(p$pelgroupidx, function(i) p$ix[[i]])), p$smdemidx) # zoop + pelagic fish + small dem
  allgrididx = c(p$ixR, p$ixFish) #(1:ncol(p$theta))
  
  for (ilgdem in p$lgdemidx) {
    pelpreyidx = (pelgrididx)[which(p$theta[ilgdem, pelgrididx] != 0)]
    allpreyidx = (allgrididx)[which(p$theta[ilgdem, ] != 0)]
    
    # lambda = total pelagic pery B / total prey B  # 
    # lambda = (sum(u[p$ix[[1]][p$mc[p$ix[[1]]]>p$mMedium]]) + sum(u[p$ix[[2]][p$mc[p$ix[[2]]]>p$mMedium & p$mc[p$ix[[2]]]<p$mLarge]]))/
    #   (sum(u[p$ix[[1]][p$mc[p$ix[[1]]]>p$mMedium]]) + sum(u[p$ix[[2]][p$mc[p$ix[[2]]]>p$mMedium & p$mc[p$ix[[2]]]<p$mLarge]]) +
    #      sum(u[p$ix[[3]][p$mc[p$ix[[3]]]>p$mMedium & p$mc[p$ix[[3]]]<p$mLarge]]) + sum(u[3:4]))
    
    lambda = sum(u[pelpreyidx]) / sum(u[allpreyidx]) # Eq. 15
    
    eT = p$Tp * lambda + p$Tb * (1 - lambda) # effect temperature for adult demersals in shallow water (<200m) Eq. 16
    
    p$fT_dem_shallow = p$Q10 ^ ((eT - p$Tref) / 10) # demersal
    p$fT_met_dem_shallow = p$Q10m ^ ((eT - p$Tref) / 10)
    
    #demersal
    
    p$V[ilgdem] = p$fT_dem_shallow * p$Vsave[ilgdem] #
    p$Cmax[ilgdem] = p$fT_dem_shallow * p$Cmaxsave[ilgdem]
    p$metabolism[ilgdem] = p$fT_met_dem_shallow * p$metabolismsave[ilgdem]
  
}

  return(p)
}
