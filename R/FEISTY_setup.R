
#===============================================================================
# Functions to setup model applications
# (generate specific MODEL parameter lists)
#===============================================================================


#' setupBasic (Petrik et al., 2019)
#' 
#' \code{setupBasic} creates a basic three-species setup as described by Petrik et al. (2019).
#' 
#' @details The setupBasic makes a basic three-species setup (small pelagic fish, large pelagic fish, and demersal fish) as described in Petrik et al (2019). 
#' There are four resources: small zooplankton, large zooplankton, small benthos, and large benthos. Large benthos actually do \bold{not exist} (always 0).
#' 
#' @author Ken H Andersen, Karline Soetaert <karline.soetaert@nioz.nl>, Yixin Zhao
#'
#' @usage p=setupBasic(szprod = 100, lzprod = 100, bprod = 5, depth = 100, Tp = 10, Tb = 8)
#' 
#' @param szprod Small zooplankton productivity. \cr
#' Actually, this represents small zooplankton carrying capacity [gww/m2] but it will multiply the growth rate \bold{r} which is always 1 [1/yr]. 
#' Therefore, it is described as small zooplankton productivity [gww/m2/year]. \code{lzprod} and \code{bprod} are same.
#' @param lzprod Large zooplankton productivity.
#' @param bprod Benthic organism productivity.
#' @param depth Water column depth [meter]. depth>=200 is characterized as deep water, and depth<200 is characterized as shallow water. 
#' depth=300 and depth=1000 do not have any different effects on simulations.
#' @param Tp Pelagic water temperature, representing the top 100m average temperature [Celsius].
#' @param Tb Bottom water temperature [Celsius].
#' 
#' @return
#' Additional parameters added by the function \code{\link{paramInit}}:
#' \itemize{
#' \item szprod, Small zooplankton productivity, from parameter input.
#' \item lzprod, Large zooplankton productivity, from parameter input.
#' \item bprod, Benthic organism productivity, from parameter input.
#' \item depth, Water column depth, from parameter input.
#' \item Tp, Pelagic water temperature, from parameter input.
#' \item Tb, Bottom water temperature, from parameter input.
#' \item mMedium, The boundary weight (mass) between small fish (mc <= mMedium) and medium fish (mMedium < mc < mLarge).
#' \item mLarge, The boundary weight (mass) between medium fish (mMedium < mc < mLarge) and large fish (mc >= mLarge).
#'}
#' 
#' Added by the  function \code{\link{setupBasic}}:
#' \itemize{
#' \item setup, name (character) of this setup
#' \item theta, size preference matrix added manually.
#' }
#' 
#' Other parameters returned can be found in \code{\link{paramInit}}, \code{\link{paramAddResource}},
#' \code{\link{paramAddGroup}}, \code{\link{paramAddPhysiology}}, and \code{\link{paramTeffect}}.
#' 
#' @examples 
#' p=setupBasic(szprod = 200, lzprod = 150, bprod = 15, depth = 300, Tp = 10, Tb = 9)
#' sim=simulateFEISTY(p=p)
#' plotSimulation(sim)
#' 
#' @references
#' Petrik, C. M., Stock, C. A., Andersen, K. H., van Denderen, P. D., & Watson, J. R. (2019). Bottom-up drivers of global patterns of demersal, forage, and pelagic fishes. Progress in oceanography, 176, 102124.
#' 
#' @seealso
#' \code{\link{paramInit}} 	Initialize parameters for FEISTY \cr
#' \code{\link{paramAddResource}} 	Add resource parameters \cr
#' \code{\link{paramAddGroup}} 	Add parameters of one functional type \cr
#' \code{\link{paramAddPhysiology}} 	Add physiological parameters \cr
#' \code{\link{paramTeffect}} 	Add temperature effects \cr
#' \code{\link{simulateFEISTY}} The main function to run FEISTY simulations
#' 
#' @aliases setupBasic
#' 
#' @export
#' 

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

setupBasic = function(szprod = 100, # small zoo production?
                      lzprod = 100, # large zoo production?
                      bprodin  = NA,  # benthos production?
                      dfbot  = NA,  # detrital flux reaching the bottom
                      depth  = 100, # water column depth [m]
                      Tp     = 10,  # pelagic layer averaged temperature [Celsius]
                      Tb     = 8,   # bottom layer depth [Celsius]
                      eT     = TRUE)# effective T on large demersals, logical  
  {
  # benthic production calc
  if (is.na(bprodin) & is.na(dfbot)){ # if all benthic arguments are NA, assign bprod to 5
    bprod = 5; bprodin =  -1; dfbot = -1
  } else {
    if (sum(!is.na(c(bprodin, dfbot)))>1) stop('Please check "bprod" and "dfbot" input. Only one of them should be assigned values, others should be kept as "NA".')
    if (!is.na(bprodin)) {bprod = bprodin} else {bprodin = -1}
    if (!is.na(dfbot)) {bprod = dfbot*0.1} else {dfbot = -1}
  }
  
  
  # Initialize the parameters:
  param = paramInit(depth=depth, szprod=szprod, lzprod=lzprod, bprodin=bprodin, dfbot=dfbot, bprod=bprod,Tp=Tp,Tb=Tb,
                    mMedium = 0.5, mLarge = 250, eT=eT)
  
  # Add resource:
  param = paramAddResource(
        param, 
        names= c("smallZoo", "largeZoo", "smallBenthos", "largeBenthos"),
        K    = c(szprod, lzprod, bprod, 0),  # g ww/m2  - maximum resource concentration
        r    = c(1, 1, 1, 1),              # [/yr] nudging coefficient
        mc   = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500)))

  # Add fish groups:
  # mMature=NA overrides the generic psiMature-> only adult classes 50% mature
  param = paramAddGroup(param, mMin=0.001, mMax=   250, mMature=NA, 
                        mortF=c(0,0.3),      nStages=2, name="smallPel")
  
  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=NA, 
                        mortF=c(0,0.03,0.3), nStages=3, name="largePel") 
  
  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=NA, 
                        mortF=c(0,0.03,0.3), nStages=3, name="demersals")
  
  # physiology of all fish stages
  param = paramAddPhysiology(param)
  
  param=paramTeffect(param, # only for setupbasic & 2
                      Tref=10,
                      Q10=1.88,
                      Q10m=2.35)#,
                      #u=NA)

  # Add fishing mortality
  # if F=0 No further process, return the input param set
  #setparam=setFishing(param, F=0, etaF=0.05)
  
  
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
  param$theta["smallPel_2", "demersals_1"] = 1

  # Large pelagics:
  param$theta["largePel_1", "smallZoo"] = 1    
  param$theta["largePel_2", "smallZoo"] = 0.25 
  param$theta["largePel_2", "largeZoo"] = 1   
  param$theta["largePel_2", c("smallPel_1", "largePel_1", "demersals_1")] = 1 
  param$theta["largePel_3", "smallPel_2"] = 0.5 
  param$theta["largePel_3", "largePel_2"] = 1 

  # Demersals:
  param$theta["demersals_1", "smallZoo"] = 1
  param$theta["demersals_2", "smallBenthos"] = 1
  # Large demersal fish have reduced feeding preference on large-size small pelagic fish and medium-size large pelagic fish in shallow water,
  # but do not eat them in deep water (depth>=200m).
  if (param$depth < 200){ 
  param$theta["demersals_3", "smallPel_2"] = 0.75/2
  param$theta["demersals_3", "largePel_2"] = 0.75 
  }
  param$theta["demersals_3", "smallBenthos"] = 1
  param$theta["demersals_3", "demersals_2"] = 1 
param$setup="setupBasic"
  return(param)
}

#' setupBasic2 
#' 
#' \code{setupBasic2} creates a revised setup based on \code{setupBasic}.
#' 
#' @details The setupBasic2 makes a revised three-species setup (small pelagic fish, large pelagic fish, and demersal fish) based on Petrik et al. (2019). 
#' There are four resources: small zooplankton, large zooplankton, small benthos, and large benthos. Large benthos actually do \bold{not exist} (always 0).\cr
#' Main revision:
#' \itemize{
#' \item Allowing more size numbers in each functional type. See \code{\link{paramAddGroup}}.
#' \item Generalized size-based maturity
#' \deqn{maturity level = (1 + ({mc}/mMature)^{-5})^{-1}}{maturity level = (1 + (mc/mMature)^(-5))^(-1)} 
#' where `mc` is the vector containing the geometric mean size of each class of a functional type, and `mMature` is the body size with a 50\% maturity level. \cr
#' See \code{\link{paramAddGroup}}.
#' \item Generalized size-based feeding preference
#' \deqn{\theta_{i,j} = \exp\left( -\frac{(\log(\frac{mc{i}}{\beta \cdot mc{_j}}))^2}{(2 \cdot \sigma)^2} \right)}{\theta[i,j] = exp(-(log(mc[i]/(beta*mc[j])))^2/(2*sigma)^2)}
#' \deqn{\theta_{i,j} = 0, if {mc_j} > {mc_i}}{\theta[i,j] = 0, if mc[j] > mc[i]}
#' where \eqn{\theta} is the size preference for a predator i to a prey j. The  'mc[i]' and 'mc[j]' are the geometric mean size of a size class of predator `i` and prey `j`.
#' `beta` and `sigma` are the preferred predator/prey mass ratio and the width of size preference for feeding, respectively. Default 400 and 1.3. \cr
#' The second equation indicates that the predator cannot prey on an organism larger than itself.
#' See \code{\link{paramSizepref}}.
#' \item Allowing the size-based fishing mortality. See \code{\link{setFishing}}.
#' \item Further process of size preference \eqn{\theta} for feeding preference. \cr
#' Fish of same size can conduct cannibalism (except medium-size demersal fish). \cr
#' Small pelagics and large pelagics do \bold{not} prey on benthic resources and medium-size demersal fish. \cr
#' The size preference for the large pelagic fish on small pelagic fish (functional type) is reduced by multiplying a coefficient 0.5 (\eqn{\theta}*0.5). \cr
#' Medium-size demersal fish do \bold{not} prey on zooplankton and all fish (no cannibalism as well); only eat benthic resources. \cr
#' 
#' Shallow water and deep water:
#' \itemize{
#' \item In shallow water (depth<200): \cr
#' Large demersal fish eat everything. They have reduced feeding preference on all small pelagics (\eqn{\theta}*0.5*0.75) and all large pelagics (\eqn{\theta}*0.75). \cr
#' Large demersal fish have feeding preference on zooplankton although the values are small. \cr
#' 
#' \item In deep water (depth>=200): \cr
#' Large-size demersal fish do not eat pelagic prey, only eat benthic resources, medium- and large-size demersals.
#' 
#' }
#' 
#' Small, medium, and large fish are differentiated according to `mMedium = 0.5gww` and `mLarge = 250gww`. \cr
#' `mMedium`: The boundary weight (mass) between small fish (mc <= mMedium) and medium fish (mMedium < mc < mLarge). \cr
#' `mLarge`: The boundary weight (mass) between medium fish (mMedium < mc < mLarge) and large fish (mc >= mLarge).
#' 
#' }
#' 
#' @author Ken H Andersen, Karline Soetaert <karline.soetaert@nioz.nl>, Yixin Zhao
#'
#' @usage p=setupBasic2(szprod = 100, lzprod = 100, bprod = 5, depth = 100, Tp = 10, Tb = 8, 
#' nStages=9, etaMature=0.25, F=0, etaF=0.05)
#' 
#' @param szprod Small zooplankton productivity. \cr
#' Actually, this represents small zooplankton carrying capacity [gww/m2] but it will multiply the growth rate \bold{r} which is always 1 [1/yr]. 
#' Therefore, it is described as small zooplankton productivity [gww/m2/year]. \code{lzprod} and \code{bprod} are same.
#' @param lzprod Large zooplankton productivity.
#' @param bprod Benthic organism productivity.
#' @param depth Water column depth [meter]. depth>=200 is characterized as deep water, and depth<200 is characterized as shallow water. 
#' depth=300 and depth=1000 do not have any different effects on simulations.
#' @param Tp Pelagic water temperature, representing the top 100m average temperature [Celsius].
#' @param Tb Bottom water temperature [Celsius].
#' @param nStages size number of large fish functional groups (e.g., large pelagic fish, demersal fish, and bathypelagic fish). 
#' The size number of small fish functional groups (e.g., small  pelagic fish and mesopelagic fish) is \code{round(2/3*nStages)}. 
#' Generally, \code{nStages} is multiples of 3 (e.g., \code{nStages = 3, 6, 9, or 12}...).
#' @param etaMature The coefficient determines the fish size \code{mMature} with a 50\% maturity level. 
#' \code{mMature = etaMature * mMax},  where \code{mMax} is the largest fish size (boundary) of a fish functional group. See \code{\link{paramAddGroup}}. 
#' In van Denderen et al. (2021), it was 0.002.
#' @param F Baseline fishing mortality [1/year]. \cr
#' If \code{F} is 0, there is no fishing mortality.\cr 
#' If \code{F} is assigned a value greater than 0, fishing mortality will be set by multiplying the fishing selectivity \code{psi} which is based on a S-shape function.
#' See source code of \code{\link{setFishing}}.
#' @param etaF The coefficient determining the fish size \code{mFishing} with 50\% fishing selectivity. See source code of \code{\link{setFishing}}.
#' 
#' @return 
#' Additional parameters added by the function \code{\link{paramInit}}:
#' \itemize{
#' \item szprod, Small zooplankton productivity, from parameter input.
#' \item lzprod, Large zooplankton productivity, from parameter input.
#' \item bprod, Benthic organism productivity, from parameter input.
#' \item depth, Water column depth, from parameter input.
#' \item Tp, Pelagic water temperature, from parameter input.
#' \item Tb, Bottom water temperature, from parameter input.
#' \item etaMature, The coefficient determines the fish size \code{mMature} with a 50\% maturity level, from parameter input.
#' \item mMedium, The boundary weight (mass) between small fish (mc <= mMedium) and medium fish (mMedium < mc < mLarge).
#' \item mLarge, The boundary weight (mass) between medium fish (mMedium < mc < mLarge) and large fish (mc >= mLarge).
#'}
#'
#' Added by the function \code{\link{paramSizepref}}:
#' \itemize{
#' \item theta, the size preference matrix
#' }
#'
#' Added by the function \code{\link{setupBasic2}}:
#' \itemize{
#' \item setup, name (character) of this setup
#' }
#' 
#' Other parameters returned can be found in \code{\link{paramInit}}, \code{\link{paramAddResource}}, \code{\link{paramAddGroup}},
#' \code{\link{paramAddPhysiology}}, \code{\link{paramTeffect}}, \code{\link{setFishing}}, and \code{\link{paramSizepref}}.
#' 
#' @examples 
#' p=setupBasic2(szprod = 200, lzprod = 150, bprod = 15, depth = 300, Tp = 10, Tb = 9, 
#' nStages=6, etaMature=0.25, F=0,etaF=0.05)
#' sim=simulateFEISTY(p=p)
#' plotSimulation(sim)
#' 
#' @references
#' Petrik, C. M., Stock, C. A., Andersen, K. H., van Denderen, P. D., & Watson, J. R. (2019). Bottom-up drivers of global patterns of demersal, forage, and pelagic fishes. Progress in oceanography, 176, 102124.
#' 
#' @seealso
#' \code{\link{paramInit}} 	Initialize parameters for FEISTY \cr
#' \code{\link{paramAddResource}} 	Add resource parameters \cr
#' \code{\link{paramAddGroup}} 	Add parameters of one functional type \cr
#' \code{\link{paramAddPhysiology}} 	Add physiological parameters \cr
#' \code{\link{paramSizepref}} 	Size preference matrix calculation \cr
#' \code{\link{paramTeffect}} 	Add temperature effects \cr
#' \code{\link{setFishing}} 	Set fishing mortality \cr
#' \code{\link{simulateFEISTY}} The main function to run FEISTY simulations
#'
#' @aliases setupBasic2
#' @export
#' 

# ------------------------------------------------------------------------------
# Make a basic three-species setup based on setupBasic(), but generalised to:
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

setupBasic2 = function(szprod = 100, # small zoo production?
                       lzprod = 100, # large zoo production?
                       bprodin  = NA,   # benthos production?
                       dfbot  = NA,  # detrital flux reaching the bottom
                       depth  = 100, # water column depth [m]
                       Tp     = 10,  # pelagic layer averaged temperature [Celsius]
                       Tb     = 8,  # bottom layer depth [Celsius]
                       nStages=9,
                       etaMature=0.25,
                       F=0,
                       etaF=0.05,
                       eT     = TRUE) { # effective T on large demersals, logical)
  # benthic production calc
  if (is.na(bprodin) & is.na(dfbot)){ # if all benthic arguments are NA, assign bprod to 5
    bprod = 5; bprodin = -1; dfbot = -1
   } else {
    if (sum(!is.na(c(bprodin, dfbot)))>1) stop('Please check "bprod" and "dfbot" input. Only one of them should be assigned values, others should be kept as "NA".')
    if (!is.na(bprodin)) {bprod = bprodin} else {bprodin = -1}
    if (!is.na(dfbot)) {bprod = dfbot*0.1} else {dfbot = -1}
   }
  
    
  # Initialize the parameters:
  param = paramInit(depth=depth, szprod=szprod, lzprod=lzprod, bprodin=bprodin, dfbot=dfbot, bprod=bprod, Tp=Tp,Tb=Tb,etaMature=etaMature,
                    mMedium = 0.5, mLarge = 250, eT=eT)

  # Setup resource groups:
  param = paramAddResource(
        param, 
        names= c("smallZoo", "largeZoo", "smallBenthos", "largeBenthos"),
        K    = c(szprod, lzprod, bprod, 0),  # g ww/m2  - maximum resource concentration
        r    = c(1, 1, 1, 1),              # [/yr] nudging coefficient
        mLower = c(2e-06,0.001, 0.5e-03, 0.25), # weight lower limit
        mUpper = c(0.001, 0.5, 125, 125),
        mc   = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500)))


# Add fish groups:
  nSmall = round(0.66*nStages)
  param = paramAddGroup(param, mMin=0.001, mMax=   250, mMature=etaMature*250, 
                        mortF=0,      nStages=nSmall, name="smallPel")
  
  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=etaMature*125000, 
                        mortF=0, nStages=nStages, name="largePel") 
  
  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=etaMature*125000, 
                        mortF=0, nStages=nStages, name="demersals")

  # physiology of all fish stages
  param = paramAddPhysiology(param)
  
  param=paramTeffect(param, # only for setupbasic & 2
                     Tref=10,
                     Q10=1.88,
                     Q10m=2.35)#,
                     #u=NA)
  
  # Add fishing mortality
  param=setFishing(param, F=F, etaF=etaF)

  # Setup size interaction matrix:
  thetaA = 0.5  # Large fish pref for medium forage fish
  thetaD = 0.75 # Pref of large demersal on pelagic prey
  
  #param$theta = matrix(nrow=param$nStages, ncol=param$nStages, data=0)
  #rownames(param$theta) <- colnames(param$theta) <- param$stagenames

   # Size-based interactions:  
   beta  = 400
   sigma = 1.3
   
   param$theta=paramSizepref(p=param,           # parameter settings 
                               beta = 400,  # preferred predator/prey mass ratio
                               sigma = 1.3, # width of size preference for feeding
                               type = 1)
   
  # update preference from fish to resources with the simple function
   # for (i in param$ixFish) {
   #   for(j in param$ixR){
   #     param$theta[i,j] = exp( -(log(param$mc[i]/(beta*param$mc[j])))^2 / (2*sigma)^2  )
   #   if (param$mc[j] >param$mc[i]) param$theta[i, j] = 0
   #   }
   # }
   
   
   #
   # Setup interactions between groups and resources:
   #
   ixSmall = param$ix[[1]]
   ixLarge = param$ix[[2]]
   ixDem   = param$ix[[3]]
  
   ixSmallSizeDem = ixDem[ (param$mc[ixDem]<=param$mMedium) ]
   ixMediumSizeDem = ixDem[ (param$mc[ixDem]>param$mMedium) &
                            (param$mc[ixDem]<param$mLarge) ]
   ixLargeSizeDem = ixDem[ (param$mc[ixDem]>=param$mLarge) ]
   
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
   
   # Large demersal fish have reduced feeding preference on all small pelagic fish and all large pelagic fish in shallow water
   if(param$depth<200){
          param$theta[ixLargeSizeDem, ixSmall] = thetaA * thetaD * param$theta[ixLargeSizeDem,ixSmall] 
          param$theta[ixLargeSizeDem, ixLarge] = thetaD * param$theta[ixLargeSizeDem, ixLarge] 
        
          #param$theta[ixLargeSizeDem, 1:2] = 
          #param$theta[ixLargeSizeDem, ixSmallSizeDem] =
   }else{ # Large-size demersal fish do not eat pelagic prey in deep water (depth>=200m).
          param$theta[ixLargeSizeDem, ixSmall] = 0
          param$theta[ixLargeSizeDem, ixLarge] = 0 
     
          param$theta[ixLargeSizeDem, 1:2] = 0
          param$theta[ixLargeSizeDem, ixSmallSizeDem] = 0
   }

  param$metabolism[is.na(param$metabolism)]=0
  param$mortF[is.na(param$mortF)]=0
  param$psiMature[is.na(param$psiMature)]=0
  param$z[is.na(param$z)]=0
  param$setup="setupBasic2"
  return(param)
}

#' setupVertical (van Denderen et al., 2021)
#' 
#' \code{setupVertical} creates a basic five-species setup with vertical distribution as described in van Denderen et al. (2021).
#' 
#' @details The setupVertical makes a basic five-species setup (small pelagic fish, mesopelagic fish, large pelagic fish, bathypelagic fish, and demersal fish) as described in van Denderen et al. (2021). 
#' There are four resources: small zooplankton, large zooplankton, small benthos, and large benthos. Large benthos actually do \bold{not exist} (always 0).
#' 
#' @author Ken H Andersen, Karline Soetaert <karline.soetaert@nioz.nl>, Yixin Zhao
#'
#' @usage p=setupVertical(szprod = 80, lzprod = 80, bent = 150, region = 4, depth = 1500, photic = 150)
#' 
#' @param szprod Small zooplankton productivity. \cr
#' Actually, this represents small zooplankton carrying capacity [gww/m2] but it will multiply the growth rate \bold{r} which is always 1 [1/yr]. 
#' Therefore, it is described as small zooplankton productivity [gww/m2/year]. \code{lzprod} is the same.
#' @param lzprod Large zooplankton productivity.
#' @param bent Detrital flux out of the photic zone [gww/m2/year].
#' \code{bent} will be further calculated based on the Martin curve to get detrital flux reaching the bottom for driving benthic communities. See source code of \code{setupVertical}.
#' @param region Different regions: 1 Tropical, 2 Temperate, 3 Boreal, 4 Default 10 Celsius.
#' It represents the water column temperature profile for three regions. 
#' The default is 10 Celcius for the whole water column (\code{region = 4}). The source file is in .../data/tempdata.dat. It is the same dataset used in van Denderen et al. (2021).
#' @param depth  water column depth [meter]. \cr 
#' Different \code{depth} values will influence fish vertical overlap and temperature-dependent physiological rates. See source code of \code{setupVertical}
#' @param photic Photic zone depth [meter]. The value affects the diel vertical migration depth. See source code of \code{setupVertical}.
#' 
#' @return
#' Additional parameters added by function \code{\link{paramInit}}:
#' \itemize{
#' \item szprod, Small zooplankton productivity, from parameter input.
#' \item lzprod, Large zooplankton productivity, from parameter input.
#' \item bent,  Detrital flux out of the photic zone, from parameter input.
#' \item bottom, Water column depth, from parameter input (depth).
#' \item photic, Photic zone depth, from parameter input.
#' \item mesop, mesopelagic zone depth.
#' \item visual, visual=1.5: visual predator. visual=1: non-visual predator.
#' \item etaMature, The coefficient determines the fish size \code{mMature} with a 50\% maturity level. 
#' It is a constant 0.002, following van Denderen et al. (2021).
#' \item region, water region index, from parameter input.
#'}
#'
#' Added by function \code{\link{paramSizepref}}:
#' \itemize{
#' \item sizepref, the size preference matrix for each predator x to each prey y.
#' }
#'
#' Added by function \code{\link{setupVertical}}:
#' \itemize{
#' \item setup, name (character) of this setup
#' \item dvm, diel vertical migration depth [m]
#' \item depthDay, a matrix containing vertical distribution data during daytime for each resource and size class (column) in water (row)
#' \item dayout, a matrix containing overlap data during daytime for each predator x to each prey y
#' \item depthNight, a matrix containing vertical distribution data during the night for each resource and size class (column) in water (row)
#' \item nightout, a matrix containing overlap data during the night for each predator x to each prey y
#' \item vertover, the average vertical overlap matrix for each predator x to each prey y. `(dayout+nightout)/2`
#' \item theta, the feeding preference matrix for each predator x to each prey y. It is the product of `sizeprefer` and `vertover`.
#' }
#' 
#' Other parameters returned can be found in \code{\link{paramInit}}, \code{\link{paramAddResource}}, \code{\link{paramAddGroup}}, and
#' \code{\link{paramAddPhysiology}}.
#' 
#' 
#' @references
#' van Denderen, P. D., Petrik, C. M., Stock, C. A., & Andersen, K. H. (2021). Emergent global biogeography of marine fish food webs. Global Ecology and Biogeography, 30(9), 1822-1834.
#' 
#' @examples 
#' p=setupVertical(szprod = 200, lzprod = 150, bent = 100, region = 1, depth = 1000, photic = 120)
#' sim=simulateFEISTY(p=p)
#' plotSimulation(sim)
#' 
#' @seealso
#' \code{\link{paramInit}} 	Initialize parameters for FEISTY \cr
#' \code{\link{paramAddResource}} 	Add resource parameters \cr
#' \code{\link{paramAddGroup}} 	Add parameters of one functional type \cr
#' \code{\link{paramAddPhysiology}} 	Add physiological parameters \cr
#' \code{\link{paramSizepref}} 	Size preference matrix calculation \cr
#' \code{\link{simulateFEISTY}} The main function to run FEISTY simulations
#' 
#' @aliases setupVertical
#' @export
#' 

# ------------------------------------------------------------------------------
# Make a basic four-species setup based up setupBasic(), but generalised to:
# - distinguish between visual and twilight predators
# - with vertical zooplankton distribution
#
# Out:
#  An updated parameter list. The feeding preferences are quite complex
#
# ------------------------------------------------------------------------------

setupVertical = function(szprod= 80,lzprod = 80, # Pelagic productivities
                         bprodin  = NA, # benthos production
                         dfbot  = NA, # detrital flux reaching the bottom
                         dfpho  = NA, # detrital flux out of photic zone
                         #nStages=6, # No. of size groups    it is 6 in van Denderen et al., 2020
                         region = 4, # Temperature profile regions: 1 Tropical, 2 Temperate, 3 Boreal, 4 Default 10 Celsius 
                         depth=1500, # Bottom depth
                         photic=150 # Photic zone depth
                         ){
  # benthic production calc
  if (is.na(bprodin) & is.na(dfbot) & is.na(dfpho)){ # if all benthic arguments are NA, assign bprod to 5
    bprodin = -1; dfbot = -1; dfpho = 150
    bprod=0.1*(dfpho*(depth/photic)^-0.86)
    if(bprod>=0.1*dfpho) bprod=0.1*dfpho
  } else {
    if (sum(!is.na(c(bprodin, dfbot, dfpho)))>1) stop('Please check "bprod" and "dfbot" input. Only one of them should be assigned values, others should be kept as "NA".')
    if (!is.na(bprodin)) {bprod = bprodin} else {bprodin = -1}
    if (!is.na(dfbot)) {bprod = dfbot*0.1} else {dfbot = -1}
    if (!is.na(dfpho)) {bprod=0.1*(dfpho*(depth/photic)^-0.86); if(bprod>=0.1*dfpho) bprod=0.1*dfpho} else {dfpho = -1}
  }
  
  #------------------  
  # Initialize the parameters:
  # habitat and small benthos
  #------------------  
  etaMature=0.002
  # bprod=0.1*(bent*(depth/photic)^-0.86)
  # if(bprod>=0.1*bent)bprod=0.1*bent

  param = paramInit(bottom=depth, szprod=szprod, lzprod=lzprod, photic=photic,
                    mesop=250, visual=1.5, bprodin=bprodin, dfbot=dfbot, dfpho=dfpho, bprod=bprod, etaMature=etaMature,region=region)
  
  #------------------  
  # Setup resource groups:
  #------------------  

    param = paramAddResource(
        param, 
        names= c("smallZoo", "largeZoo", "smallBenthos", "largeBenthos"),
        K    = c(szprod, lzprod, bprod, 0),  # g ww/m2  - maximum resource concentration
        r    = c(1, 1, 1, 1),              # [/yr] nudging coefficient
        mc   = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500)),
        mLower = c(2e-06,0.001, 0.5e-03, 0.25), # weight lower limit
        mUpper = c(0.001, 0.5, 125, 125),
        u0     = c(0.5,0.5,0.5,0))

  #------------------  
  # Add fish groups:
  #------------------  
  nStages=6
  nSmall = round(0.66*nStages)
  # mMature=NA overrides the generic psiMature-> only adult classes 50% mature
  u0  = 0.0001
  param = paramAddGroup(param, mMin=0.001, mMax=   250, mMature=etaMature*250, u0=u0,
                        mortF=0,      nStages=nSmall, name="smallPel")
  

  u0M = u0  # initial condition = 0 if no mesopelagic zone
  if (param$bottom <= param$mesop) u0M <- 0
  
  param = paramAddGroup(param, mMin=0.001, mMax=   250, mMature=etaMature*250, u0=u0M,
                        mortF=0,   nStages=nSmall, name="mesoPel")

  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=etaMature*125000, u0=u0, 
                        mortF=0, nStages=nStages, name="largePel") 
  
  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=etaMature*125000, u0=u0M, 
                        mortF=0, nStages=nStages, name="bathyPel") 

  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=etaMature*125000, u0=u0,
                        mortF=0, nStages=nStages, name="demersals")
  #param$mortF[length(param$mortF)]=0.5
  
  #------------------  
  # Setup physiology:
  #------------------  
  param = paramAddPhysiology(param,am = 0.2*20) # 20% * Max. consumption coefficient)
  
  
  #overwrite psiMature only for setupVertical
  nsize= nStages+1
  sizes = 10^(linspace(log10(0.001), log10(125000), nsize))
  matstageS =   which.min(abs(sizes - etaMature*250))
  matstageL =   which.min(abs(sizes - etaMature*125000))
  param$psiMature=param$psiMature*0
  param$psiMature[param$ix[[1]][matstageS:max(param$ix[[1]])]]=0.5
  param$psiMature[param$ix[[2]][matstageS:max(param$ix[[2]])]]=0.5
  param$psiMature[param$ix[[3]][matstageL:max(param$ix[[3]])]]=0.5
  param$psiMature[param$ix[[4]][matstageL:max(param$ix[[4]])]]=0.5
  param$psiMature[param$ix[[5]][matstageL:max(param$ix[[5]])]]=0.5
    
  #------------------  
  #------------------  
  # theta (preferences):
  #------------------  
  #------------------  
  #param$theta  = matrix(nrow=param$nStages, ncol=param$nStages, data=0)
  #rownames(param$theta) <- colnames(param$theta) <- param$stagenames

  beta  = 400
  sigma = 1.3

  param$vertover   = matrix(nrow=param$nStages, ncol=param$nStages, data=0)
  
  # calculate size-preference matrix
  param$sizeprefer=paramSizepref(p=param,           # parameter settings 
                              beta = 400,  # preferred predator/prey mass ratio
                              sigma = 1.3, # width of size preference for feeding
                              type = 2)

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

  ixjuv = which.min(abs(sizes-0.5))# - etaMature*250)) # -0.5))
  ixadult = which.min(abs(sizes-250))# - etaMature*125000)) # -250))
    
  # ixjuv = which.min(abs(param$mLower[param$ix[[5]]] - etaMature*250))
  # ixadult = which.min(abs(param$mLower[param$ix[[5]]] - etaMature*125000))

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

  xlocvec[ixadult:length(xlocvec)] = demmig #param$dvm ### or demmig???
  dem_d =  VertDist(sigmap[ix], xlocvec)

  #if shallower than euphotic depth, adult demersals feed across-habitats
  if (param$bottom <= param$photic) {
  dem_d = (dem_d + dem_n)/2
  dem_n = dem_d
  }
  
  # calculate overlap during day
  param$depthDay = matrix(nrow=length(xrange), ncol=param$nStages, data=0)
  test     = matrix(nrow=length(xrange), ncol=param$nStages, data=0)
  param$dayout = matrix(nrow=param$nStages, ncol=param$nStages, data=0)
  
  param$depthDay[, 1:2] = zp_d # resources
  param$depthDay[, 3:4] = bent_dn # resources
  param$depthDay[, param$ix[[1]]] = spel_dn
  param$depthDay[, param$ix[[2]]] = mpel_d
  param$depthDay[, param$ix[[3]]] = lpel_d
  param$depthDay[, param$ix[[4]]] = bpel_d
  param$depthDay[, param$ix[[5]]] = dem_d
  
  for (i in 1: param$nStages) {
    for ( j in 1: param$nStages) {
     test[, j] = pmin(param$depthDay[, i], param$depthDay[, j])
    }
  param$dayout[, i] = colSums(test)
  }
  
  # calculate overlap during night
  param$depthNight = matrix(nrow=length(xrange), ncol=param$nStages, data=0)
  # test will be overwritten
  param$nightout = matrix(nrow=param$nStages, ncol=param$nStages, data=0)
  
  param$depthNight[, 1:2] = zp_n # resources
  param$depthNight[, 3:4] = bent_dn # resources
  param$depthNight[, param$ix[[1]]] = spel_dn
  param$depthNight[, param$ix[[2]]] = mpel_n
  param$depthNight[, param$ix[[3]]] = lpel_n
  param$depthNight[, param$ix[[4]]] = bpel_n
  param$depthNight[, param$ix[[5]]] = dem_n
  
  for (i in 1: param$nStages) {
    for ( j in 1: param$nStages) {
      test[, j] = pmin(param$depthNight[, i], param$depthNight[, j])
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
  idx_be = param$ixFish[1]: (param$ix[[5]][1] + (ixjuv - 2)) # all pelagic and small demersals
  param$theta[idx_be, 3:4] = 0 # all pelagic and small demersals do not eat benthos,
                               # only small & large demersals eat benthos
  
  # small demersals are less preyed on
  idx_smd = (param$ix[[5]][1] + (ixjuv - 1)): (param$ix[[5]][1] + (ixadult - 2)) #
  param$theta[idx_be, idx_smd] = param$theta[idx_be, idx_smd]*0.25
  
  # medium & large demersals do not eat zooplankton
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
  
  # update temperature
  tempdata=read.table(system.file("data", "tempdata.dat", package = "FEISTY"), sep=',') #
  tempdata[,5]=10 #
  Q10=1.88
  Q10m=1.88
  
  dist=(param$depthDay+param$depthNight)/2
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
  param$setup="setupVertical"
  
return(param)  
}

#' setupVertical2
#' 
#' \code{setupVertical2} creates a revised setup based on \code{setupVertical}.
#' 
#' @details The setupVertical2 makes a revised five-species setup (small pelagic fish, mesopelagic fish, large pelagic fish, bathypelagic fish, and demersal fish) with vertical distribution based on van Denderen et al. (2021). 
#' There are four resources: small zooplankton, large zooplankton, small benthos, and large benthos. Large benthos actually do \bold{not exist} (always 0).\cr
#' Main revision:
#' \itemize{
#' \item Allowing more size numbers in each functional type. See \code{\link{paramAddGroup}}.
#' \item Generalized size-based maturity
#' \deqn{maturity level = (1 + ({mc}/mMature)^{-5})^{-1}}{maturity level = (1 + (mc/mMature)^(-5))^(-1)} 
#' where `mc` is the vector containing the geometric mean size of each class of a functional type, and `mMature` is the body size with a 50\% maturity level. \cr
#' See \code{\link{paramAddGroup}}.
# \item Generalized size-based feeding preference.
#' \item Allowing the size-based fishing mortality. See \code{\link{setFishing}}.
#' 
#' }
#' 
#' @author Ken H Andersen, Karline Soetaert <karline.soetaert@nioz.nl>, Yixin Zhao
#'
#' @usage p=setupVertical2(szprod = 80, lzprod = 80, bent = 150, nStages = 6, region = 4, depth = 1500, photic = 150, 
#' mesopelagic = 250, visual = 1.5, etaMature = 0.25, F = 0, etaF=0.05)
#' 
#' @param szprod Small zooplankton productivity. \cr
#' Actually, this represents small zooplankton carrying capacity [gww/m2] but it will multiply the growth rate \bold{r} which is always 1 [1/yr]. 
#' Therefore, it is described as small zooplankton productivity [gww/m2/year]. \code{lzprod} is the same.
#' @param lzprod Large zooplankton productivity.
#' @param bent Detrital flux out of the photic zone [gww/m2/year].
#' \code{bent} will be further calculated based on the Martin curve to get detrital flux reaching the bottom for driving benthic communities. See source code of \code{setupVertical}.
#' @param region Different regions: 1 Tropical, 2 Temperate, 3 Boreal, 4 Default 10 Celsius.
#' It represents the water column temperature profile for three regions. 
#' The default is 10 Celcius for the whole water column (\code{region = 4}). The source file is in .../data/tempdata.dat. It is the same dataset used in van Denderen et al. (2021).
#' @param depth  water column depth [meter]. \cr 
#' Different \code{depth} values will influence fish vertical overlap and temperature-dependent physiological rates. See source code of \code{setupVertical}
#' @param photic Photic zone depth [meter]. The value affects the diel vertical migration depth. See source code of \code{setupVertical}.
#' @param mesopelagic Mesopelagic depth [meter]. The value affects the diel vertical migration depth. See source code of \code{setupVertical} or \code{setupVertical2}.
#' @param visual \code{visual=1.5}: visual predator, predation ability is enhanced during the day and decreased in the twilight zone during the day. \cr 
#' \code{visual=1}: non-visual predator, predation abilities are equal in day and night. \cr
#' It must \bold{be careful} to assign other values to \code{visual}, or the setup could crash. See source code of \code{setupVertical} or \code{setupVertical2}.
#' @param etaMature The coefficient determines the fish size \code{mMature} with a 50\% maturity level. \code{mMature = etaMature * mMax},  where \code{mMax} is the largest fish size (boundary) of a fish functional group. See \code{\link{paramAddGroup}}. 
#' In van Denderen et al. (2021), it was 0.002.
#' @param F Baseline fishing mortality [1/year]. \cr
#' If \code{F} is 0, there is no fishing mortality. \cr
#' If \code{F} is assigned a value greater than 0, fishing mortality will be set by multiplying the fishing selectivity \code{psi} which is based on a S-shape function. See source code of \code{\link{setFishing}}.
#' @param etaF the coefficient determining the fish size \code{mFishing} with 50\% fishing selectivity. See source code of \code{\link{setFishing}}.
#' 
#' @return
#' Additional parameters added by function \code{\link{paramInit}}:
#' \itemize{
#' \item szprod, Small zooplankton productivity, from parameter input.
#' \item lzprod, Large zooplankton productivity, from parameter input.
#' \item bent,  Detrital flux out of the photic zone, from parameter input.
#' \item bottom, Water column depth, from parameter input (depth).
#' \item photic, Photic zone depth, from parameter input.
#' \item mesop, mesopelagic zone depth. from parameter input.
#' \item visual, visual=1.5: visual predator. visual=1: non-visual predator. from parameter input.
#' \item etaMature, The coefficient determines the fish size \code{mMature} with a 50\% maturity level, from parameter input.
#' \item region, water region index, from parameter input.
#'}
#'
#' Added by function \code{\link{paramSizepref}}:
#' \itemize{
#' \item sizepref, the size preference matrix for each predator x to each prey y.
#' }
#'
#' Added by function \code{\link{setupVertical}}:
#' \itemize{
#' \item setup, name (character) of this setup
#' \item dvm, diel vertical migration depth [m]
#' \item depthDay, a matrix containing vertical distribution data during daytime for each resource and size class (column) in water (row)
#' \item dayout, a matrix containing overlap data during daytime for each predator x to each prey y
#' \item depthNight, a matrix containing vertical distribution data during the night for each resource and size class (column) in water (row)
#' \item nightout, a matrix containing overlap data during the night for each predator x to each prey y
#' \item vertover, the average vertical overlap matrix for each predator x to each prey y. `(dayout+nightout)/2`
#' \item theta, the feeding preference matrix for each predator x to each prey y. It is the product of `sizeprefer` and `vertover`.
#' }
#' 
#' Other parameters returned can be found in \code{\link{paramInit}}, \code{\link{paramAddResource}}, \code{\link{paramAddGroup}},
#' \code{\link{paramAddPhysiology}}, and \code{\link{setFishing}}.
#' 
#' @examples 
#' p=setupVertical2(szprod = 200, lzprod = 150, bent = 100, nStages = 6, region = 1, depth = 1000, photic = 120,
#' mesopelagic = 250, visual = 1.5, etaMature = 0.25, F = 0, etaF = 0.05)
#' sim=simulateFEISTY(p=p)
#' plotSimulation(sim)
#' 
#' @references
#' van Denderen, P. D., Petrik, C. M., Stock, C. A., & Andersen, K. H. (2021). Emergent global biogeography of marine fish food webs. Global Ecology and Biogeography, 30(9), 1822-1834.
#' 
#' @seealso
#' \code{\link{paramInit}} 	Initialize parameters for FEISTY \cr
#' \code{\link{paramAddResource}} 	Add resource parameters \cr
#' \code{\link{paramAddGroup}} 	Add parameters of one functional type \cr
#' \code{\link{paramAddPhysiology}} 	Add physiological parameters \cr
#' \code{\link{paramSizepref}} 	Size preference matrix calculation \cr
#' \code{\link{setFishing}} 	Set fishing mortality \cr
#' \code{\link{simulateFEISTY}} The main function to run FEISTY simulations
#' 
#' @aliases setupVertical2
#' @export
#' 

# ------------------------------------------------------------------------------
# Make a basic four-species setup based up setupBasic(), but generalised to:
# - distinguish between visual and twilight predators
# - with vertical zooplankton distribution
#
# Out:
#  An updated parameter list. The feeding preferences are quite complex
#
# ------------------------------------------------------------------------------

setupVertical2 = function(szprod= 80,lzprod = 80, # Pelagic productivities
                         bprodin  = NA, # benthos production
                         dfbot  = NA, # detrital flux reaching the bottom
                         dfpho  = NA, # detrital flux out of photic zone
                         nStages=6, # No. of size groups
                         region = 4, # Temperature profile regions: 1 Tropical, 2 Temperate, 3 Boreal, 4 Default 10 Celsius 
                         depth=1500, # Bottom depth
                         photic=150, # Photic zone depth
                         mesopelagic=250, # mesopelagic depth
                         visual=1.5,# >1 visual predation primarily during the day, = 1 equal day and night
                         etaMature = 0.25, # Size of matureation relative to
                         # asymptotic size. Different from
                         # van Denderen (2021), where it is 0.002
                         F=0,
                         etaF=0.05) {
  # benthic production calc
  if (is.na(bprodin) & is.na(dfbot) & is.na(dfpho)){ # if all benthic arguments are NA, assign bprod to 5
    bprodin = -1; dfbot = -1; dfpho = 150
    bprod=0.1*(dfpho*(depth/photic)^-0.86)
    if(bprod>=0.1*dfpho) bprod=0.1*dfpho
  } else {
    if (sum(!is.na(c(bprodin, dfbot, dfpho)))>1) stop('Please check "bprod" and "dfbot" input. Only one of them should be assigned values, others should be kept as "NA".')
    if (!is.na(bprodin)) {bprod = bprodin} else {bprodin = -1}
    if (!is.na(dfbot)) {bprod = dfbot*0.1} else {dfbot = -1}
    if (!is.na(dfpho)) {bprod=0.1*(dfpho*(depth/photic)^-0.86); if(bprod>=0.1*dfpho) bprod=0.1*dfpho} else {dfpho = -1}
  }
    
  #------------------  
  # Initialize the parameters:
  # habitat and small benthos
  #------------------  
  # bprod=0.1*(bent*(depth/photic)^-0.86)
  # if(bprod>=0.1*bent)bprod=0.1*bent
  
  param = paramInit(bottom=depth, szprod=szprod, lzprod=lzprod, photic=photic,
                    mesop=mesopelagic, visual=visual, bprodin=bprodin, dfbot=dfbot, dfpho=dfpho, bprod=bprod, etaMature=etaMature,region=region)
  
  #------------------  
  # Setup resource groups:
  #------------------  
  
  param = paramAddResource(
    param, 
    names= c("smallZoo", "largeZoo", "smallBenthos", "largeBenthos"),
    K    = c(szprod, lzprod, bprod, 0),  # g ww/m2  - maximum resource concentration
    r    = c(1, 1, 1, 1),              # [/yr] nudging coefficient
    mc   = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500)),
    mLower = c(2e-06,0.001, 0.5e-03, 0.25), # weight lower limit
    mUpper = c(0.001, 0.5, 125, 125),
    u0     = c(0.5,0.5,0.5,0))
  
  #------------------  
  # Add fish groups:
  #------------------  
  nSmall = round(0.66*nStages)
  # mMature=NA overrides the generic psiMature-> only adult classes 50% mature
  u0  = 0.0001
  param = paramAddGroup(param, mMin=0.001, mMax=   250, mMature=etaMature*250, u0=u0,
                        mortF=0,      nStages=nSmall, name="smallPel")
  
  
  u0M = u0  # initial condition = 0 if no mesopelagic zone
  if (param$bottom <= param$mesop) u0M <- 0
  
  param = paramAddGroup(param, mMin=0.001, mMax=   250, mMature=etaMature*250, u0=u0M,
                        mortF=0,   nStages=nSmall, name="mesoPel")
  
  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=etaMature*125000, u0=u0, 
                        mortF=0, nStages=nStages, name="largePel") 
  
  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=etaMature*125000, u0=u0M, 
                        mortF=0, nStages=nStages, name="bathyPel") 
  
  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=etaMature*125000, u0=u0,
                        mortF=0, nStages=nStages, name="demersals")
  #param$mortF[length(param$mortF)]=0.5
  
  #------------------  
  # Setup physiology:
  #------------------  
  param = paramAddPhysiology(param,am = 0.2*20) # 20% * Max. consumption coefficient)
  
  # Add fishing mortality
  param=setFishing(param, F=F, etaF=etaF)
  
  #------------------  
  #------------------  
  # theta (preferences):
  #------------------  
  #------------------  
  #param$theta  = matrix(nrow=param$nStages, ncol=param$nStages, data=0)
  #rownames(param$theta) <- colnames(param$theta) <- param$stagenames
  
  beta  = 400
  sigma = 1.3

  param$vertover   = matrix(nrow=param$nStages, ncol=param$nStages, data=0)
  
  # calculate size-preference matrix
  param$sizeprefer=paramSizepref(p=param,           # parameter settings 
                                   beta = 400,  # preferred predator/prey mass ratio
                                   sigma = 1.3, # width of size preference for feeding
                                   type = 1)
  
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
  
  ixjuv = which.min(abs(param$mLower[param$ix[[5]]] - 0.5))# which.min(abs(param$mLower[param$ix[[5]]] - etaMature*250))
  ixadult = which.min(abs(param$mLower[param$ix[[5]]] - 250))# which.min(abs(param$mLower[param$ix[[5]]] - etaMature*125000))
  
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
  
  xlocvec[ixadult:length(xlocvec)] = demmig #param$dvm ### or demmig???
  dem_d =  VertDist(sigmap[ix], xlocvec)
  
  #if shallower than euphotic depth, adult demersals feed across-habitats
  if (param$bottom <= param$photic) {
    dem_d = (dem_d + dem_n)/2
    dem_n = dem_d
  }
  
  # calculate overlap during day
  param$depthDay = matrix(nrow=length(xrange), ncol=param$nStages, data=0)
  test     = matrix(nrow=length(xrange), ncol=param$nStages, data=0)
  param$dayout = matrix(nrow=param$nStages, ncol=param$nStages, data=0)
  
  param$depthDay[, 1:2] = zp_d # resources
  param$depthDay[, 3:4] = bent_dn # resources
  param$depthDay[, param$ix[[1]]] = spel_dn
  param$depthDay[, param$ix[[2]]] = mpel_d
  param$depthDay[, param$ix[[3]]] = lpel_d
  param$depthDay[, param$ix[[4]]] = bpel_d
  param$depthDay[, param$ix[[5]]] = dem_d
  
  for (i in 1: param$nStages) {
    for ( j in 1: param$nStages) {
      test[, j] = pmin(param$depthDay[, i], param$depthDay[, j])
    }
    param$dayout[, i] = colSums(test)
  }
  
  # calculate overlap during night
  param$depthNight = matrix(nrow=length(xrange), ncol=param$nStages, data=0)
  # test will be overwritten
  param$nightout = matrix(nrow=param$nStages, ncol=param$nStages, data=0)
  
  param$depthNight[, 1:2] = zp_n # resources
  param$depthNight[, 3:4] = bent_dn # resources
  param$depthNight[, param$ix[[1]]] = spel_dn
  param$depthNight[, param$ix[[2]]] = mpel_n
  param$depthNight[, param$ix[[3]]] = lpel_n
  param$depthNight[, param$ix[[4]]] = bpel_n
  param$depthNight[, param$ix[[5]]] = dem_n
  
  for (i in 1: param$nStages) {
    for ( j in 1: param$nStages) {
      test[, j] = pmin(param$depthNight[, i], param$depthNight[, j])
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
  idx_be = param$ixFish[1]: (param$ix[[5]][1] + (ixjuv - 2)) # all pelagic and small demersals
  param$theta[idx_be, 3:4] = 0 # all pelagic and small demersals do not eat benthos,
  # only medium & adult demersals eat benthos
  
  # medium demersals are less preyed on
  idx_smd = (param$ix[[5]][1] + (ixjuv - 1)): (param$ix[[5]][1] + (ixadult - 2)) #
  param$theta[idx_be, idx_smd] = param$theta[idx_be, idx_smd]*0.25
  
  # small & large demersals do not eat zooplankton
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
  
  # update temperature
  tempdata=read.table(system.file("data", "tempdata.dat", package = "FEISTY"), sep=',') #
  tempdata[,5]=10 #
  Q10=1.88
  Q10m=1.88
  
  dist=(param$depthDay+param$depthNight)/2
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
  
  param$setup="setupVertical2"
  
  return(param)  
}

# ------------------------------------------------------------------------------
# Make a basic setup with just pelagic fish, 
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


