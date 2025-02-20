
#===============================================================================
# Functions to setup model applications
# (generate specific MODEL parameter lists)
#===============================================================================


#' setupBasic
#' 
#' \code{setupBasic} creates a basic three-species setup as described by Petrik et al. (2019).
#' 
#' @details The setupBasic makes a basic three-species setup (small pelagic fish, large pelagic fish, and demersal fish) as described in Petrik et al. (2019). 
#' There are three resources: small mesozooplankton, large mesozooplankton, benthos, and a spare position in the state variable vector.
#' 
#' @author Ken H. Andersen, Karline Soetaert, Yixin Zhao
#'
#' @usage setupBasic(szprod = 100, 
#'                   lzprod = 100, 
#'                   bprodin = NA, 
#'                   dfbot = NA, 
#'                   depth = 100, 
#'                   Tp = 10, 
#'                   Tb = 8)
#' 
#' @param szprod Small mesozooplankton productivity. \cr
#' The parameter represents small mesozooplankton carrying capacity [g/m2] and multiplied with the growth rate \bold{r}, which is always 1 [1/yr], it gives the maximum resource productivity [g/m2/year]. 
#' Therefore, it is termed small mesozooplankton productivity [g/m2/year]. \code{lzprod} and \code{bprod} are same.
#' @param lzprod Large mesozooplankton productivity [g/m2/year]. 
#' @param bprodin Benthic productivity input [g/m2/year]. Default NA. Input either of `bprodin` or `dfbot`.
#' The benthic productivity `bprod` equals `bprodin`.
#' @param dfbot Detrital flux reaching the bottom [g/m2/year]. Default NA. Input either of `bprodin` or `dfbot`.
#' It will multiply the trophic transfer efficiency (10\%) to get the benthic productivity `bprod`. If both are NAs then `bprod = 5`.
#' @param depth Water column depth [meter]. depth>=200 is characterized as deep water, and depth<200 is characterized as shallow water. 
#' depth=300 and depth=1000 do not have any different effects on simulations.
#' @param Tp Pelagic water temperature, representing the top 100m average temperature [Celsius].
#' @param Tb Bottom water temperature [Celsius].
#' 
#' @return
#' Additional parameters added by the function \code{\link{paramInit}}:
#' \itemize{
#' \item szprod, Small mesozooplankton productivity, from parameter input.
#' \item lzprod, Large mesozooplankton productivity, from parameter input.
#' \item bprodin, Benthic productivity input, from parameter input. If input is NA, the returned value is -1 for passing to FORTRAN.
#' \item dfbot, Detrital flux reaching the bottom, from parameter input. If input is NA, the returned value is -1 for passing to FORTRAN.
#' \item bprod, Benthic productivity input, from calculation based on `bprodin` or `dfbot`.
#' \item depth, Water column depth, from parameter input.
#' \item Tp, Pelagic water temperature, from parameter input.
#' \item Tb, Bottom water temperature, from parameter input.
#' \item mMedium, The boundary weight (mass) between small fish (mc <= mMedium) and medium fish (mMedium < mc < mLarge).
#' \item mLarge, The boundary weight (mass) between medium fish (mMedium < mc < mLarge) and large fish (mc >= mLarge).
#' \item bET, Logical flag, which is always TRUE as described in Petrik et al., 2019. See `updateET` source code in FEISTY_parms.R.
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
#' p=setupBasic(szprod = 200, 
#'              lzprod = 150,
#'              bprodin = 15, 
#'              dfbot = NA, 
#'              depth = 300, 
#'              Tp = 10, 
#'              Tb = 9)
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
setupBasic = function(szprod = 100, # small zoo production
                      lzprod = 100, # large zoo production
                      bprodin = NA, # benthos production
                      dfbot  = NA,  # detrital flux reaching the bottom
                      depth  = 100, # water column depth [m]
                      Tp     = 10,  # pelagic layer averaged temperature [Celsius]
                      Tb     = 8)   # bottom layer depth [Celsius]
{
  # benthic production calc
  if (is.na(bprodin) & is.na(dfbot)){ # if all benthic arguments are NA, assign bprod to 5
    bprod = 5; bprodin =  -1; dfbot = -1
  } else {
    if (sum(!is.na(c(bprodin, dfbot)))>1) stop('Please check "bprodin" and "dfbot" input. Only one of them should be assigned values, others should be kept as "NA".')
    if (!is.na(bprodin)) {bprod = bprodin} else {bprodin = -1}
    if (!is.na(dfbot)) {bprod = dfbot*0.1} else {dfbot = -1}
  }
  
  # Initialize the parameters:
  param = paramInit(depth=depth, szprod=szprod, lzprod=lzprod, bprodin=bprodin, dfbot=dfbot, bprod=bprod,Tp=Tp,Tb=Tb,
                    mMedium = 0.5, mLarge = 250, bET=TRUE) # bET: boolean, effective T on large demersals in Petrik et al., 2019
  
  # Add resource:
  param = paramAddResource(
    param, 
    names= c("smallZoo", "largeZoo", "benthos", "Spare_position"),
    K    = c(szprod, lzprod, bprod, 0),  # g ww/m2  - maximum resource concentration
    r    = c(1, 1, 1, 1),              # [/yr] nudging coefficient
    mc   = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500)))
  
  # Add fish groups:
  # mMature=NA overrides the generic psiMature-> only adult classes 50% mature
  u0  = 1E-5
  param = paramAddGroup(param, mMin=0.001, mMax=   250, mMature=NA, u0=u0,
                        mortF=0,      nStages=2, name="smallPel") #mortF=c(0,0.03,0.3)
  
  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=NA, u0=u0,
                        mortF=0, nStages=3, name="largePel") 
  
  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=NA, u0=u0,
                        mortF=0, nStages=3, name="demersals")
  
  # physiology of all fish stages
  param = paramAddPhysiology(param)
  
  param=paramTeffect(param, # only for setupbasic & 2
                     Tref=10,
                     Q10=1.88,
                     Q10m=2.35,
                     pelgroupidx=c(1:(param$nGroups-1)),
                     demgroupidx=param$nGroups)
  
  # Add fishing mortality
  # Has been assigned by param = paramAddGroup(..., mortF=c()), hard coded
  
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
  param$theta["demersals_2", "benthos"] = 1
  # Large demersal fish have reduced feeding preference on large-size small pelagic fish and medium-size large pelagic fish in shallow water,
  # but do not eat them in deep water (depth>=200m).
  if (param$depth < 200){ 
    param$theta["demersals_3", "smallPel_2"] = 0.75/2
    param$theta["demersals_3", "largePel_2"] = 0.75 
  }
  param$theta["demersals_3", "benthos"] = 1
  param$theta["demersals_3", "demersals_2"] = 1 
  param$setup="setupBasic"
  return(param)
}

#' setupBasic2 
#' 
#' \code{setupBasic2} creates a revised setup based on \code{setupBasic}.
#' 
#' @details The setupBasic2 makes a revised three-species setup (small pelagic fish, large pelagic fish, and demersal fish) based on Petrik et al. (2019). 
#' There are three resources: small mesozooplankton, large mesozooplankton, benthos, and a spare position in the state variable vector. \cr
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
#' Small, medium, and large fish are differentiated according to `mMedium = 0.5g` and `mLarge = 250g`. \cr
#' `mMedium`: The boundary weight (mass) between small fish (mc <= mMedium) and medium fish (mMedium < mc < mLarge). \cr
#' `mLarge`: The boundary weight (mass) between medium fish (mMedium < mc < mLarge) and large fish (mc >= mLarge).
#' 
#' }
#' 
#' @author Ken H. Andersen, Karline Soetaert, Yixin Zhao
#'
#' @usage setupBasic2(szprod = 100, 
#'                    lzprod = 100, 
#'                    bprodin = NA, 
#'                    dfbot = NA, 
#'                    depth = 100, 
#'                    Tp = 10, 
#'                    Tb = 8, 
#'                    nStages=9, 
#'                    etaMature=0.25, 
#'                    Fmax=0, 
#'                    etaF=0.05, 
#'                    bET=TRUE)
#' 
#' @param szprod Small mesozooplankton productivity. \cr
#' The parameter represents small mesozooplankton carrying capacity [g/m2] and multiplied with the growth rate \bold{r}, which is always 1 [1/yr], it gives the maximum resource productivity [g/m2/year].
#' Therefore, it is described as small mesozooplankton productivity [g/m2/year]. \code{lzprod} and \code{bprod} are same.
#' @param lzprod Large mesozooplankton productivity.
#' @param bprodin Benthic productivity input [g/m2/year]. Default NA. Input either of `bprodin` or `dfbot`.
#' The benthic productivity `bprod` equals `bprodin`.
#' @param dfbot Detrital flux reaching the bottom [g/m2/year]. Default NA. Input either of `bprodin` or `dfbot`.
#' It will multiply the trophic transfer efficiency (10\%) to get the benthic productivity `bprod`. If both are NAs then `bprod = 5`.
#' @param depth Water column depth [meter]. depth>=200 is characterized as deep water, and depth<200 is characterized as shallow water. 
#' depth=300 and depth=1000 do not have any different effects on simulations.
#' @param Tp Pelagic water temperature, representing the top 100m average temperature [Celsius].
#' @param Tb Bottom water temperature [Celsius].
#' @param nStages Size number of large fish functional groups (e.g., large pelagic fish, demersal fish, and midwater predators). 
#' The size number of small fish functional groups (e.g., small pelagic fish and mesopelagic fish) is \code{round(2/3*nStages)}. 
#' Generally, \code{nStages} is multiples of 3 (e.g., \code{nStages = 3, 6, 9, or 12}...).
#' @param etaMature The coefficient determines the fish size \code{mMature} with a 50\% maturity level. 
#' \code{mMature = etaMature * mMax},  where \code{mMax} is the largest fish size (boundary) of a fish functional group. See \code{\link{paramAddGroup}}. 
#' In van Denderen et al. (2021), it was 0.002.
#' @param Fmax Maximum fishing mortality [1/year]. \cr
#' If \code{Fmax} is 0, there is no fishing mortality.\cr 
#' If \code{Fmax} is assigned a value greater than 0, fishing mortality will be set by multiplying the fishing selectivity \code{psi} which is based on a S-shape function. See source code of \code{\link{setFishing}}.\cr
#' Note here it only allows assigning fishing mortality to all functional types based on the `Fmax` input. 
#' If users want to assign fishing mortality to specific functional types with different `Fmax` and `etaF`, it can be done by calling \code{\link{setFishing}} later on. Details are in \code{\link{setFishing}}. 
#' @param etaF The coefficient determining the fish size \code{mFishing} with 50\% fishing selectivity. See source code of \code{\link{setFishing}}.
#' @param bET Logical flag, controlling whether turns on the effective temperature effects on the large demersal fish. See Petrik et al., 2019. See `updateET` source code in FEISTY_parms.R.
#' 
#' @return 
#' Additional parameters added by the function \code{\link{paramInit}}:
#' \itemize{
#' \item szprod, Small mesozooplankton productivity, from parameter input.
#' \item lzprod, Large mesozooplankton productivity, from parameter input.
#' \item bprodin, Benthic productivity input, from parameter input. If input is NA, the returned value is -1 for passing to FORTRAN.
#' \item dfbot, Detrital flux reaching the bottom, from parameter input. If input is NA, the returned value is -1 for passing to FORTRAN.
#' \item bprod, Benthic productivity, from calculation based on `bprodin` or `dfbot`.
#' \item depth, Water column depth, from parameter input.
#' \item Tp, Pelagic water temperature, from parameter input.
#' \item Tb, Bottom water temperature, from parameter input.
#' \item etaMature, The coefficient determines the fish size \code{mMature} with a 50\% maturity level, from parameter input.
#' \item mMedium, The boundary weight (mass) between small fish (mc <= mMedium) and medium fish (mMedium < mc < mLarge).
#' \item mLarge, The boundary weight (mass) between medium fish (mMedium < mc < mLarge) and large fish (mc >= mLarge).
#' \item bET, Logical flag, TRUE or FALSE. Effective temperature effects on the large demersal fish.
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
#' Other returned parameters can be found in \code{\link{paramInit}}, \code{\link{paramAddResource}}, \code{\link{paramAddGroup}},
#' \code{\link{paramAddPhysiology}}, \code{\link{paramTeffect}}, \code{\link{setFishing}}, and \code{\link{paramSizepref}}.
#' 
#' @examples 
#' p=setupBasic2(szprod = 200, 
#'               lzprod = 150, 
#'               bprodin = 15, 
#'               dfbot = NA, 
#'               depth = 300, 
#'               Tp = 10, 
#'               Tb = 9, 
#'               nStages=6, 
#'               etaMature=0.25, 
#'               Fmax=0, 
#'               etaF=0.05, 
#'               bET=TRUE)
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
setupBasic2 = function(szprod = 100, # small zoo production?
                       lzprod = 100, # large zoo production?
                       bprodin = NA, # benthos production?
                       dfbot  = NA,  # detrital flux reaching the bottom
                       depth  = 100, # water column depth [m]
                       Tp     = 10,  # pelagic layer averaged temperature [Celsius]
                       Tb     = 8,  # bottom layer depth [Celsius]
                       nStages = 9,
                       etaMature = 0.25,
                       Fmax = 0,
                       etaF = 0.05,
                       bET  = TRUE) { # boolean, effective T on large demersals)
  # benthic production calc
  if (is.na(bprodin) & is.na(dfbot)){ # if all benthic arguments are NA, assign bprod to 5
    bprod = 5; bprodin = -1; dfbot = -1
  } else {
    if (sum(!is.na(c(bprodin, dfbot)))>1) stop('Please check "bprodin" and "dfbot" input. Only one of them should be assigned values, others should be kept as "NA".')
    if (!is.na(bprodin)) {bprod = bprodin} else {bprodin = -1}
    if (!is.na(dfbot)) {bprod = dfbot*0.1} else {dfbot = -1}
  }
  
  # Initialize the parameters:
  param = paramInit(depth=depth, szprod=szprod, lzprod=lzprod, bprodin=bprodin, dfbot=dfbot, bprod=bprod, Tp=Tp,Tb=Tb,etaMature=etaMature,
                    mMedium = 0.5, mLarge = 250, bET=bET)
  
  # Setup resource groups:
  param = paramAddResource(
    param, 
    names= c("smallZoo", "largeZoo", "benthos", "Spare_position"),
    K    = c(szprod, lzprod, bprod, 0),  # g ww/m2  - maximum resource concentration
    r    = c(1, 1, 1, 1),              # [/yr] nudging coefficient
    mLower = c(2e-06,0.001, 0.5e-03, 0.25), # weight lower limit
    mUpper = c(0.001, 0.5, 125, 125),
    mc   = c(2e-06*sqrt(500), 0.001*sqrt(500), 1e-04*sqrt(250000), 0.25*sqrt(500)))
  
  # Add fish groups:
  nSmall = round(0.66*nStages)
  u0  = 1E-5
  param = paramAddGroup(param, mMin=0.001, mMax=   250, mMature=etaMature*250, u0=u0,
                        mortF=0,      nStages=nSmall, name="smallPel")
  
  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=etaMature*125000, u0=u0,
                        mortF=0, nStages=nStages, name="largePel") 
  
  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=etaMature*125000, u0=u0,
                        mortF=0, nStages=nStages, name="demersals")
  
  # physiology of all fish stages
  param = paramAddPhysiology(param)
  
  param=paramTeffect(param, # only for setupbasic & 2
                     Tref=10,
                     Q10=1.88,
                     Q10m=2.35,
                     pelgroupidx=c(1:(param$nGroups-1)),
                     demgroupidx=param$nGroups)
  
  # Add fishing mortality
  param=setFishing(param, Fmax=Fmax, etaF=etaF)
  
  # Setup size interaction matrix:
  thetaA = 0.5  # Large fish pref for medium forage fish
  thetaD = 0.75 # Pref of large demersal on pelagic prey
  
  # Size-based interactions:  
  param$theta=paramSizepref(p=param,           # parameter settings 
                            beta = 400,  # preferred predator/prey mass ratio
                            sigma = 1.3, # width of size preference for feeding
                            type = 1)
  
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
  param$theta[ixSmallSizeDem, 3:4] = 0
  
  # ... or on medium-sized demersal fish:
  param$theta[ixSmall, ixMediumSizeDem] = 0
  param$theta[ixLarge, ixMediumSizeDem] = 0
  
  # Large pelagics have reduced feeding efficiency on small pelagics:
  param$theta[ixLarge,ixSmall] = thetaA * param$theta[ixLarge,ixSmall] 
  # ... and do not feed on medium-sized demersal:
  param$theta[ixLarge, ixMediumSizeDem ] = 0
  
  # Medium-size large demersals only feed on benthos and do cannibalism:
  param$theta[ixMediumSizeDem, 1:2] = 0 
  param$theta[ixMediumSizeDem, param$ixFish[!(param$ixFish %in% ixMediumSizeDem)]] = 0 
  
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
  
  param$setup="setupBasic2"
  
  return(param)
}

#' setupVertical
#' 
#' \code{setupVertical} creates a basic five-species setup with vertical distribution as described in van Denderen et al. (2021).
#' 
#' @details The setupVertical makes a basic five-species setup (small pelagic fish, mesopelagic fish, large pelagic fish, midwater predators, and demersal fish) as described in van Denderen et al. (2021). 
#' There are three resources: small mesozooplankton, large mesozooplankton, benthos, and a spare position in the state variable vector. \cr
#' 
#' @author Ken H. Andersen, Karline Soetaert, Yixin Zhao
#'
#' @usage setupVertical(szprod = 80, 
#'                      lzprod = 80, 
#'                      bprodin = NA, 
#'                      dfbot = NA, 
#'                      dfpho = NA, 
#'                      region = 4, 
#'                      depth = 800, 
#'                      photic = 150)
#' 
#' @param szprod Small mesozooplankton productivity. \cr
#' The parameter represents small mesozooplankton carrying capacity [g/m2] and multiplied with the growth rate \bold{r}, which is always 1 [1/yr], it gives the maximum resource productivity [g/m2/year].
#' Therefore, it is described as small mesozooplankton productivity [g/m2/year]. \code{lzprod} is the same.
#' @param lzprod Large mesozooplankton productivity.
#' @param bprodin Benthic productivity input [g/m2/year]. Default NA. Input either of `bprodin`, `dfbot` or `dfpho`.
#' The benthic productivity `bprod` equals `bprodin`.
#' @param dfbot Detrital flux reaching the bottom [g/m2/year]. Default NA. Input either of `bprodin`, `dfbot` or `dfpho`.
#' It will multiply the trophic transfer efficiency (10\%) to get the benthic productivity `bprod`.
#' @param dfpho Detrital flux out of the photic zone [g/m2/year]. Default NA. Input either of `bprodin`, `dfbot` or `dfpho`. If all are NAs then `dfpho = 150`.
#' `dfpho` will be further calculated based on the Martin curve to get detrital flux reaching the bottom and then multiplied the trophic transfer efficiency (10\%) to get benthic productivity `bprod` ultimately .\cr
#' See source code of \code{setupVertical}.
#' @param region Different regions: 1 Tropical, 2 Temperate, 3 Boreal, 4 Default 10 Celsius.
#' It represents the water column temperature profile for three regions. 
#' The default is 10 Celcius for the whole water column (\code{region = 4}). It is the same dataset used in van Denderen et al. (2021).
#' @param depth  water column depth [meter]. \cr 
#' Different \code{depth} values will influence fish vertical overlap and temperature-dependent physiological rates. See source code of \code{setupVertical}
#' @param photic Photic zone depth [meter]. The value affects the diel vertical migration depth and the detrital flux reaching the bottom (if `dfpho` is used). See source code of \code{setupVertical}.
#' 
#' @return
#' Additional parameters added by function \code{\link{paramInit}}:
#' \itemize{
#' \item szprod, Small mesozooplankton productivity, from parameter input.
#' \item lzprod, Large mesozooplankton productivity, from parameter input.
#' \item bprodin, Benthic productivity input, from parameter input. If input is NA, the returned value is -1 for passing to FORTRAN.
#' \item dfbot, Detrital flux reaching the bottom, from parameter input. If input is NA, the returned value is -1 for passing to FORTRAN.
#' \item dfpho, Detrital flux out of the photic zone, from parameter input. If input is NA, the returned value is -1 for passing to FORTRAN.
#' \item bprod, Benthic productivity, from calculation based on `bprodin` or `dfbot`, or `dfpho`.
#' \item bottom, Water column depth, from parameter input (depth).
#' \item photic, Photic zone depth, from parameter input.
#' \item shelfdepth, Continental shelf depth. 250m, cannot be changed in setupVertical. 
#'       If the water is shallow (bottom<=shelfdepth), mesopelagic fish and midwater predators do not exist,
#'       and fish do not conduct diel vertical migration (dvm=0). See source code of \code{setupVertical} or \code{setupVertical2}.
#' \item visual, visual=1.5: visual predator. visual=1: non-visual predator.
#' \item etaMature, The coefficient determines the fish size \code{mMature} with a 50\% maturity level. 
#' It is a constant 0.002, following van Denderen et al. (2021).
#' \item region, Water region index, from parameter input.
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
#' \item ixmedium, an index indicating where medium size fish start.
#' \item ixlarge, an index indicating where large size fish start. E.g., ixmedium = 4, ixlarge = 7: number 1 to 3 represent small fish, number 4 to 6 represent medium fish, number 7 to the last size class represent large fish
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
#' @references
#' van Denderen, P. D., Petrik, C. M., Stock, C. A., & Andersen, K. H. (2021). Emergent global biogeography of marine fish food webs. Global Ecology and Biogeography, 30(9), 1822-1834.
#' 
#' @examples 
#' p=setupVertical(szprod = 200, 
#'                 lzprod = 150, 
#'                 bprodin = NA, 
#'                 dfbot = NA, 
#'                 dfpho = 100, 
#'                 region = 1, 
#'                 depth = 1000, 
#'                 photic = 120)
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
setupVertical = function(szprod = 80, # small zoo production
                         lzprod = 80, # large zoo production
                         bprodin = NA, # benthos production
                         dfbot  = NA, # detrital flux reaching the bottom
                         dfpho  = NA, # detrital flux out of photic zone
                         #nStages=6, # No. of size groups    it is 6 in van Denderen et al., 2020
                         region = 4, # Temperature profile regions: 1 Tropical, 2 Temperate, 3 Boreal, 4 Default 10 Celsius 
                         depth = 800, # Bottom depth
                         photic = 150 # Photic zone depth
){
  # benthic production calc
  if (is.na(bprodin) & is.na(dfbot) & is.na(dfpho)){ # if all benthic arguments are NA, assign bprod to 5
    bprodin = -1; dfbot = -1; dfpho = 350
    bprod=0.1*(dfpho*(depth/photic)^-0.86)
    if(bprod>=0.1*dfpho) bprod=0.1*dfpho
  } else {
    if (sum(!is.na(c(bprodin, dfbot, dfpho)))>1) stop('Please check "bprodin", "dfbot" and "dfpho" input. Only one of them should be assigned values, others should be kept as "NA".')
    if (!is.na(bprodin)) {bprod = bprodin} else {bprodin = -1}
    if (!is.na(dfbot)) {bprod = dfbot*0.1} else {dfbot = -1}
    if (!is.na(dfpho)) {bprod=0.1*(dfpho*(depth/photic)^-0.86); if(bprod>=0.1*dfpho) bprod=0.1*dfpho} else {dfpho = -1}
  }
  
  #------------------  
  # Initialize the parameters:
  # habitat and small benthos
  #------------------  
  etaMature=0.002
  
  param = paramInit(bottom=depth, szprod=szprod, lzprod=lzprod, photic=photic,
                    shelfdepth=250, visual=1.5, bprodin=bprodin, dfbot=dfbot, dfpho=dfpho, bprod=bprod,
                    etaMature=etaMature,region=region)
  
  #------------------  
  # Setup resource groups:
  #------------------  
  
  param = paramAddResource(
    param, 
    names= c("smallZoo", "largeZoo", "benthos", "Spare_position"),
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
  if (param$bottom <= param$shelfdepth) u0M <- 0
  
  param = paramAddGroup(param, mMin=0.001, mMax=   250, mMature=etaMature*250, u0=u0M,
                        mortF=0,   nStages=nSmall, name="mesoPel")
  
  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=etaMature*125000, u0=u0, 
                        mortF=0, nStages=nStages, name="largePel") 
  
  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=etaMature*125000, u0=u0M, 
                        mortF=0, nStages=nStages, name="midwPred") 
  
  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=etaMature*125000, u0=u0,
                        mortF=0, nStages=nStages, name="demersals")
  #param$mortF[length(param$mortF)]=0.5
  
  #------------------  
  # Setup physiology:
  #------------------  
  param = paramAddPhysiology(param)
  
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
  # theta (preferences):
  #------------------  
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
  
  if (param$bottom <= param$shelfdepth) 
    param$dvm = 0              # no migration in shallow habitats
  
  ixmedium = which.min(abs(sizes-0.5))# - etaMature*250)) # -0.5))
  ixlarge = which.min(abs(sizes-250))# - etaMature*125000)) # -250))
  # ixmedium = which.min(abs(param$mLower[param$ix[[5]]] - etaMature*250))
  # ixlarge = which.min(abs(param$mLower[param$ix[[5]]] - etaMature*125000))
  param$ixmedium=ixmedium
  param$ixlarge=ixlarge
  
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
  mpel_d = VertDist(sigmap[ix], xloc=param$dvm)
  
  ## large pelagic fish night (all at surface)
  ix = param$ix[[3]]
  lpel_n = VertDist(sigmap[ix], xloc=0)
  
  # large pelagic fish day (non-large at surface   large half at surface half at dvm)
  xlocvec = rep(0,length(ix)) 
  xlocvec[ixlarge:length(xlocvec)] = param$dvm 
  lpel_d = VertDist(sigmap[ix], xloc=xlocvec)
  lpel_d = (lpel_d + lpel_n)/2
  
  ## bathypelagic night (large in midwater, others at surface)
  ix = param$ix[[4]]
  xlocvec = rep(0,length(ix)) # initialization
  xlocvec[ixlarge:length(xlocvec)] = param$dvm # non-large at surface   large at dvm
  bpel_n = VertDist(sigmap[ix], xloc=xlocvec)
  
  # bathypelagic day (all at dvm)
  bpel_d = VertDist(sigmap[ix], xloc=param$dvm)
  
  ## demersal fish night
  ix = param$ix[[5]]
  xlocvec = rep(0,length(ix)) # initialization
  xlocvec[ixmedium:length(xlocvec)] = param$bottom # small at surface   medium and large at bottom
  dem_n = VertDist(sigmap[ix], xlocvec)
  
  # demersal fish day; small at surface/ medium at bottom/ large at middle
  demmig = param$dvm # ? from matlab
  if ((param$bottom - param$dvm) >= 1200) 
    demmig = param$dvm + (param$bottom-param$dvm-1200)
  if ((param$bottom - param$dvm) >= 1500)
    demmig = param$bottom
  
  dem_d= matrix(nrow=length(xrange), ncol=length(param$ix[[5]]), data=0)
  
  xlocvec[ixlarge:length(xlocvec)] = demmig #param$dvm ### or demmig???
  dem_d =  VertDist(sigmap[ix], xlocvec)
  
  #if shallower than euphotic depth, large demersals feed across-habitats
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
  pelpred = pelpred[ixlarge:length(pelpred)] # large large pelagic  11  at dvm during day
  preytwi = c(param$ix[[2]], param$ix[[4]])  # mesopelagic 7 8   bathypelagic 12 13 14
  param$dayout[pelpred, preytwi] = param$dayout[pelpred, preytwi]/param$visual*(2 - param$visual)  # /1.5 to restore  then *0.5 
  
  # average overlap during the whole day
  param$vertover = (param$dayout + param$nightout)*0.5
  
  # calculate combined feeding preference matrix
  param$theta = param$sizeprefer*param$vertover
  
  #  specific revision of feeding preference
  idx_be = param$ixFish[1]: (param$ix[[5]][1] + (ixmedium - 2)) # all pelagic and small demersals
  param$theta[idx_be, 3:4] = 0 # all pelagic and small demersals do not eat benthos,
  # only small & large demersals eat benthos
  
  # medium demersals are less preyed on
  idx_smd = (param$ix[[5]][1] + (ixmedium - 1)): (param$ix[[5]][1] + (ixlarge - 2)) #
  param$theta[idx_be, idx_smd] = param$theta[idx_be, idx_smd]*0.25
  
  # medium & large demersals do not eat zooplankton
  param$theta[(param$ix[[5]][1] + (ixmedium - 1)) : param$ix[[5]][length(param$ix[[5]])], 1:2] = 0
  
  # provide benefit to forage and mesopelagic fish (predator avoidance)
  pred1 = (param$ix[[3]][1]+ (ixlarge - 1)) : param$ix[[3]][length(param$ix[[3]])]
  pred2 = (param$ix[[4]][1]+ (ixlarge - 1)) : param$ix[[4]][length(param$ix[[4]])]
  pred3 = (param$ix[[5]][1]+ (ixlarge - 1)) : param$ix[[5]][length(param$ix[[5]])]
  prey1 = (param$ix[[1]][1]+ (ixmedium   - 1)) : param$ix[[1]][length(param$ix[[1]])]
  prey2 = (param$ix[[2]][1]+ (ixmedium   - 1)) : param$ix[[2]][length(param$ix[[2]])]
  idx_predat = c(pred1, pred2, pred3)
  idx_prey   = c(prey1, prey2)
  param$theta[idx_predat,idx_prey] = param$theta[idx_predat,idx_prey]*0.5
  
  # update temperature
  tempdata=read.table(system.file("extdata", "tempdata.dat", package = "FEISTY"), sep=',') #
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
#' @details The setupVertical2 makes a revised five-species setup (small pelagic fish, mesopelagic fish, large pelagic fish, midwater predators, and demersal fish) with vertical distribution based on van Denderen et al. (2021). 
#' There are three resources: small mesozooplankton, large mesozooplankton, benthos, and a spare position in the state variable vector. \cr
#' Main revision:
#' \itemize{
#' \item Allowing more size numbers in each functional type. See \code{\link{paramAddGroup}}.
#' \item Generalized size-based maturity
#' \deqn{maturity level = (1 + ({mc}/mMature)^{-5})^{-1}}{maturity level = (1 + (mc/mMature)^(-5))^(-1)} 
#' where `mc` is the vector containing the geometric mean size of each class of a functional type, and `mMature` is the body size with a 50\% maturity level. \cr
#' See \code{\link{paramAddGroup}}.
# \item Generalized size-based feeding preference.
#' \item Allowing the size-based fishing mortality. See \code{\link{setFishing}}.
#' }
#' 
#' @author Ken H. Andersen, Karline Soetaert, Yixin Zhao
#'
#' @usage setupVertical2(szprod = 80, lzprod = 80, 
#'                       bprodin = NA, 
#'                       dfbot = NA, 
#'                       dfpho = NA, 
#'                       nStages = 9, 
#'                       Tp = NA,
#'                       Tm = NA, 
#'                       Tb = NA, 
#'                       depth = 800,
#'                       photic = 150, 
#'                       shelfdepth = 250, 
#'                       visual = 1.5, 
#'                       etaMature = 0.25, 
#'                       Fmax = 0, 
#'                       etaF=0.05)
#' 
#' @param szprod Small mesozooplankton productivity. \cr
#' The parameter represents small mesozooplankton carrying capacity [g/m2] and multiplied with the growth rate \bold{r}, which is always 1 [1/yr], it gives the maximum resource productivity [g/m2/year].
#' Therefore, it is described as small mesozooplankton productivity [g/m2/year]. \code{lzprod} is the same.
#' @param lzprod Large mesozooplankton productivity.
#' @param bprodin Benthic productivity input [g/m2/year]. Default NA. Input either of `bprodin`, `dfbot` or `dfpho`.
#' The benthic productivity `bprod` equals `bprodin`.
#' @param dfbot Detrital flux reaching the bottom [g/m2/year]. Default NA. Input either of `bprodin`, `dfbot` or `dfpho`.
#' It will multiply the trophic transfer efficiency (10\%) to get the benthic productivity `bprod`.
#' @param dfpho Detrital flux out of the photic zone [g/m2/year]. Default NA. Input either of `bprodin`, `dfbot` or `dfpho`. If all are NAs then `dfpho = 150`.
#' `dfpho` will be further calculated based on the Martin curve to get detrital flux reaching the bottom and then multiplied the trophic transfer efficiency (10\%) to get benthic productivity `bprod` ultimately .\cr
#' See source code of \code{setupVertical}.
#' @param nStages size number of large fish functional types (e.g., large pelagic fish, demersal fish, and midwater predators). 
#' The size number of small fish functional types (e.g., small pelagic fish and mesopelagic fish) is \code{round(2/3*nStages)}. 
#' Generally, \code{nStages} is multiples of 3 (e.g., \code{nStages = 3, 6, 9, or 12}...).
#' @param Tp Pelagic water temperature, representing the top 100m average temperature [Celsius]. Default NA. Input NA means Tp = 10.
#' @param Tm Mid-water temperature, representing the average temperature of 500m - up to 1500m. Default NA. Input NA means Tp = Tb.
#' @param Tb Bottom water (the bottom layer) temperature [Celsius]. Default NA. Input NA means Tb = 10.
#' @param depth  Water column depth [meter]. \cr 
#' Different \code{depth} values will influence fish vertical overlap and temperature-dependent physiological rates. See source code of \code{setupVertical}
#' @param photic Photic zone depth [meter]. The value affects the diel vertical migration depth and the detrital flux reaching the bottom (if `dfpho` is used). See source code of \code{setupVertical}.
#' @param shelfdepth Continental shelf depth [meter]. If the water is shallow (depth<=shelfdepth), mesopelagic fish and midwater predators do not exist,
#'                   and fish do not conduct diel vertical migration (dvm=0). See source code of \code{setupVertical} or \code{setupVertical2}.
#' @param visual \code{visual=1.5}: visual predator, predation ability is enhanced during the day and decreased in the twilight zone during the day. \cr 
#' \code{visual=1}: non-visual predator, predation abilities are equal in day and night. \cr
#' It must \bold{be careful} to assign other values to \code{visual}, or the setup could crash. See source code of \code{setupVertical} or \code{setupVertical2}.
#' @param etaMature the coefficient determines the fish size \code{mMature} with a 50\% maturity level. \code{mMature = etaMature * mMax},  where \code{mMax} is the largest fish size (boundary) of a fish functional group. See \code{\link{paramAddGroup}}. 
#' In van Denderen et al. (2021), it was 0.002.
#' @param Fmax Maximum fishing mortality [1/year]. \cr
#' If \code{Fmax} is 0, there is no fishing mortality. \cr
#' If \code{Fmax} is assigned a value greater than 0, fishing mortality will be set by multiplying the fishing selectivity \code{psi} which is based on a S-shape function. See source code of \code{\link{setFishing}}.
#' Note here it only allows assigning fishing mortality to all functional types based on the `Fmax` input. 
#' If users want to assign fishing mortality to specific functional types with different `Fmax` and `etaF`, it can be done by calling \code{\link{setFishing}} later on. Details are in \code{\link{setFishing}}.
#' @param etaF the coefficient determining the fish size \code{mFishing} with 50\% fishing selectivity. See source code of \code{\link{setFishing}}.
#' 
#' @return
#' Additional parameters added by function \code{\link{paramInit}}:
#' \itemize{
#' \item szprod, Small mesozooplankton productivity, from parameter input.
#' \item lzprod, Large mesozooplankton productivity, from parameter input.
#' \item bprodin, Benthic productivity input, from parameter input. If input is NA, the returned value is -1 for passing to FORTRAN.
#' \item dfbot, Detrital flux reaching the bottom, from parameter input. If input is NA, the returned value is -1 for passing to FORTRAN.
#' \item dfpho, Detrital flux out of the photic zone, from parameter input. If input is NA, the returned value is -1 for passing to FORTRAN.
#' \item bprod, Benthic productivity, from calculation based on `bprodin` or `dfbot`, or `dfpho`.
#' \item bottom, Water column depth, from parameter input (depth).
#' \item photic, Photic zone depth, from parameter input.
#' \item shelfdepth, continental shelf depth. from parameter input.
#' \item visual, visual=1.5: visual predator. visual=1: non-visual predator. from parameter input.
#' \item etaMature, The coefficient determines the fish size \code{mMature} with a 50\% maturity level, from parameter input.
#' \item Tp, pelagic water temperature, from parameter input.
#' \item Tm, mid-water temperature, from parameter input.
#' \item Tb, bottom water temperature, from parameter input.
#'}
#'
#' Added by function \code{\link{paramSizepref}}:
#' \itemize{
#' \item sizepref, the size preference matrix for each predator x to each prey y.
#' }
#'
#' Added by function \code{\link{setupVertical2}}:
#' \itemize{
#' \item setup, name (character) of this setup
#' \item dvm, diel vertical migration depth [m]
#' \item ixmedium, an index indicating where medium size fish start.
#' \item ixlarge, an index indicating where large size fish start. E.g., ixmedium = 4, ixlarge = 7: number 1 to 3 represent small fish, number 4 to 6 represent medium fish, number 7 to the last size class represent large fish
#' \item depthDay, a matrix containing vertical distribution data during daytime for each resource and size class (column) in water (row)
#' \item dayout, a matrix containing overlap data during daytime for each predator x to each prey y
#' \item depthNight, a matrix containing vertical distribution data during the night for each resource and size class (column) in water (row)
#' \item nightout, a matrix containing overlap data during the night for each predator x to each prey y
#' \item vertover, the average vertical overlap matrix for each predator x to each prey y. `(dayout+nightout)/2`
#' \item theta, the feeding preference matrix for each predator x to each prey y. It is the product of `sizeprefer` and `vertover`.
#' }
#' 
#' Other returned parameters can be found in \code{\link{paramInit}}, \code{\link{paramAddResource}}, \code{\link{paramAddGroup}},
#' \code{\link{paramAddPhysiology}}, and \code{\link{setFishing}}.
#' 
#' @examples 
#' p=setupVertical2(szprod = 200, lzprod = 150, 
#'                  bprodin = NA, 
#'                  dfbot = NA, 
#'                  dfpho = 100, 
#'                  nStages = 9, 
#'                  Tp=18,
#'                  Tm=15,
#'                  Tb=10,
#'                  depth = 800, 
#'                  photic = 120,
#'                  shelfdepth = 250, 
#'                  visual = 1.5, 
#'                  etaMature = 0.25, 
#'                  Fmax = 0, 
#'                  etaF = 0.05)
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
setupVertical2 = function(szprod = 80, # small zoo production
                          lzprod = 80, # large zoo production
                          bprodin = NA, # benthos production
                          dfbot  = NA, # detrital flux reaching the bottom
                          dfpho  = NA, # detrital flux out of photic zone
                          nStages = 9, # No. of size groups
                          Tp = NA, # Average T of top 100 m (up to 100 m). Default 10 Celsius.
                          Tm = NA, # Average T of 500 - 1500 m (up to 1500 m). Default 10 Celsius. Keep it as NA, if no Tm data. Tm = Tb.
                          Tb = NA, # Bottom T (last layer value). Default 10 Celsius.
                          depth = 800, # Bottom depth
                          photic = 150, # Photic zone depth
                          shelfdepth = 250, # shelf region depth
                          visual = 1.5, # >1 visual predation primarily during the day, = 1 equal day and night
                          etaMature = 0.25, # Size of matureation relative to
                          # asymptotic size. Different from
                          # van Denderen (2021), where it is 0.002
                          Fmax = 0,
                          etaF = 0.05) {
  # benthic production calc
  if (is.na(bprodin) & is.na(dfbot) & is.na(dfpho)){ # if all benthic arguments are NA, assign bprod to 5
    bprodin = -1; dfbot = -1; dfpho = 350
    bprod=0.1*(dfpho*(depth/photic)^-0.86)
    if(bprod>=0.1*dfpho) bprod=0.1*dfpho
  } else {
    if (sum(!is.na(c(bprodin, dfbot, dfpho)))>1) stop('Please check "bprodin", "dfbot" and "dfpho" input. Only one of them should be assigned values, others should be kept as "NA".')
    if (!is.na(bprodin)) {bprod = bprodin} else {bprodin = -1}
    if (!is.na(dfbot)) {bprod = dfbot*0.1} else {dfbot = -1}
    if (!is.na(dfpho)) {bprod=0.1*(dfpho*(depth/photic)^-0.86); if(bprod>=0.1*dfpho) bprod=0.1*dfpho} else {dfpho = -1}
  }
  
  # Temperature initialize
  if (is.na(Tp)) Tp = 10
  if (is.na(Tb)) Tb = 10
  if (is.na(Tm)) Tm = Tb # if Tm is not provided, Tm = Tb
  
  #------------------  
  # Initialize the parameters:
  # habitat and small benthos
  #------------------  
  
  param = paramInit(bottom=depth, szprod=szprod, lzprod=lzprod, photic=photic,
                    shelfdepth=shelfdepth, visual=visual, bprodin=bprodin, dfbot=dfbot, dfpho=dfpho, bprod=bprod,
                    etaMature=etaMature, Tp=Tp, Tm=Tm, Tb=Tb)
  
  #------------------  
  # Setup resource groups:
  #------------------  
  
  param = paramAddResource(
    param, 
    names= c("smallZoo", "largeZoo", "benthos", "Spare_position"),
    K    = c(szprod, lzprod, bprod, 0),  # g ww/m2  - maximum resource concentration
    r    = c(1, 1, 1, 1),              # [/yr] nudging coefficient
    mc   = c(2e-06*sqrt(500), 0.001*sqrt(500), 1e-04*sqrt(250000), 0.25*sqrt(500)),
    mLower = c(2e-06,0.001, 1e-04, 0.25), # weight lower limit
    mUpper = c(0.001, 0.5, 25, 125),
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
  if (param$bottom <= param$shelfdepth) u0M <- 0
  
  param = paramAddGroup(param, mMin=0.001, mMax=   250, mMature=etaMature*250, u0=u0M,
                        mortF=0,   nStages=nSmall, name="mesoPel")
  
  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=etaMature*125000, u0=u0, 
                        mortF=0, nStages=nStages, name="largePel") 
  
  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=etaMature*125000, u0=u0M, 
                        mortF=0, nStages=nStages, name="midwPred") 
  
  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=etaMature*125000, u0=u0,
                        mortF=0, nStages=nStages, name="demersals")
  #param$mortF[length(param$mortF)]=0.5
  
  #------------------  
  # Setup physiology:
  #------------------  
  param = paramAddPhysiology(param)
  
  # Add fishing mortality
  param=setFishing(param, Fmax=Fmax, etaF=etaF)
  
  #------------------  
  # theta (preferences):
  #------------------  
  
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
  
  if (param$bottom <= param$shelfdepth) 
    param$dvm = 0              # no migration in shallow habitats
  
  ixmedium = which.min(abs(param$mLower[param$ix[[5]]] - 0.5))# which.min(abs(param$mLower[param$ix[[5]]] - etaMature*250))
  ixlarge = which.min(abs(param$mLower[param$ix[[5]]] - 250))# which.min(abs(param$mLower[param$ix[[5]]] - etaMature*125000))
  param$ixmedium=ixmedium
  param$ixlarge=ixlarge
  
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
  mpel_d = VertDist(sigmap[ix], xloc=param$dvm)
  
  ## large pelagic fish night (all at surface)
  ix = param$ix[[3]]
  lpel_n = VertDist(sigmap[ix], xloc=0)
  
  # large pelagic fish day (non-large at surface   large half at surface half at dvm)
  xlocvec = rep(0,length(ix)) 
  xlocvec[ixlarge:length(xlocvec)] = param$dvm 
  lpel_d = VertDist(sigmap[ix], xloc=xlocvec)
  lpel_d = (lpel_d + lpel_n)/2
  
  ## bathypelagic night (large in midwater, others at surface)
  ix = param$ix[[4]]
  xlocvec = rep(0,length(ix)) # initialization
  xlocvec[ixlarge:length(xlocvec)] = param$dvm # non-large at surface   large at dvm
  bpel_n = VertDist(sigmap[ix], xloc=xlocvec)
  
  # bathypelagic day (all at dvm)
  bpel_d = VertDist(sigmap[ix], xloc=param$dvm)
  
  ## demersal fish night
  ix = param$ix[[5]]
  xlocvec = rep(0,length(ix)) # initialization
  xlocvec[ixmedium:length(xlocvec)] = param$bottom # small at surface   medium and large at bottom
  dem_n = VertDist(sigmap[ix], xlocvec)
  
  # demersal fish day; small at surface/ medium at bottom/ large at middle
  demmig = param$dvm # ? from matlab
  if ((param$bottom - param$dvm) >= 1200) 
    demmig = param$dvm + (param$bottom-param$dvm-1200)
  if ((param$bottom - param$dvm) >= 1500)
    demmig = param$bottom
  
  dem_d= matrix(nrow=length(xrange), ncol=length(param$ix[[5]]), data=0)
  
  xlocvec[ixlarge:length(xlocvec)] = demmig #param$dvm ### or demmig???
  dem_d =  VertDist(sigmap[ix], xlocvec)
  
  #if shallower than euphotic depth, large demersals feed across-habitats
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
  pelpred = pelpred[ixlarge:length(pelpred)] # large large pelagic  11  at dvm during day
  preytwi = c(param$ix[[2]], param$ix[[4]])  # mesopelagic 7 8   bathypelagic 12 13 14
  param$dayout[pelpred, preytwi] = param$dayout[pelpred, preytwi]/param$visual*(2 - param$visual)  # /1.5 to restore  then *0.5 
  
  # average overlap during the whole day
  param$vertover = (param$dayout + param$nightout)*0.5
  
  # calculate combined feeding preference matrix
  param$theta = param$sizeprefer*param$vertover
  
  #  specific revision of feeding preference
  idx_be = param$ixFish[1]: (param$ix[[5]][1] + (ixmedium - 2)) # all pelagic and small demersals
  param$theta[idx_be, 3:4] = 0 # all pelagic and small demersals do not eat benthos,
  # only medium & large demersals eat benthos
  
  # medium demersals are less preyed on
  idx_smd = (param$ix[[5]][1] + (ixmedium - 1)): (param$ix[[5]][1] + (ixlarge - 2)) #
  param$theta[idx_be, idx_smd] = param$theta[idx_be, idx_smd]*0.25
  
  # small & large demersals do not eat zooplankton
  param$theta[(param$ix[[5]][1] + (ixmedium - 1)) : param$ix[[5]][length(param$ix[[5]])], 1:2] = 0
  
  # provide benefit to forage and mesopelagic fish (predator avoidance)
  pred1 = (param$ix[[3]][1]+ (ixlarge - 1)) : param$ix[[3]][length(param$ix[[3]])]
  pred2 = (param$ix[[4]][1]+ (ixlarge - 1)) : param$ix[[4]][length(param$ix[[4]])]
  pred3 = (param$ix[[5]][1]+ (ixlarge - 1)) : param$ix[[5]][length(param$ix[[5]])]
  prey1 = (param$ix[[1]][1]+ (ixmedium   - 1)) : param$ix[[1]][length(param$ix[[1]])]
  prey2 = (param$ix[[2]][1]+ (ixmedium   - 1)) : param$ix[[2]][length(param$ix[[2]])]
  idx_predat = c(pred1, pred2, pred3)
  idx_prey   = c(prey1, prey2)
  param$theta[idx_predat,idx_prey] = param$theta[idx_predat,idx_prey]*0.5
  
  param=paramTeffect_vet(param)
  
  param$setup="setupVertical2"
  
  return(param)  
}

#' setupTimeseries
#' 
#' \code{setupTimeseries} creates the setup for time-series simulations based on prescribed setups (setupBasic, setupBasic2, or setupVertical2, not for setupVertical).
#' 
#' @author Yixin Zhao
#'
#' @usage setupTimeseries(p = setupVertical2(),
#'                        szbio_ts = NA,
#'                        lzbio_ts = NA,
#'                        szprod_ts = NA,
#'                        lzprod_ts = NA,
#'                        bprodin_ts = NA,
#'                        dfbot_ts  = NA,
#'                        dfpho_ts  = NA,
#'                        Tp_ts = NA,
#'                        Tm_ts = NA,
#'                        Tb_ts = NA,
#'                        Fsmp_ts = NA,
#'                        Fmesop_ts = NA,
#'                        Flgp_ts = NA,
#'                        Fmidwp_ts = NA,
#'                        Fdem_ts = NA,
#'                        benthosK = 80)
#' 
#' @param p Parameter set (setupBasic, setupBasic2, or setupVertical2, not for setupVertical). 
#' Non time-varying data of a grid such as water column depth `depth` and photic zone depth `photic` should be assigned in `p = setupXX()`. Also the non-varying data defined by users should be put here. 
#' For example, there are no time-varying temperature, then the temperature constants should be defined in `p = setupXX()`.
#' @param szbio_ts Small mesozooplankton biomass time-series data [g/m2].
#' @param lzbio_ts Large mesozooplankton biomass time-series data [g/m2].
#' @param szprod_ts Small mesozooplankton productivity time-series data [g/m2/year].
#' @param lzprod_ts Large mesozooplankton productivity time-series data [g/m2/year].
#' @param bprodin_ts Large Benthic productivity time-series data [g/m2/year].
#' @param dfbot_ts Detrital flux reaching the bottom time-series data [g/m2/year]. It will multiply the trophic transfer efficiency (10\%) to get benthic productivity `bprod_ts`.
#' @param dfpho_ts Detrital flux out of the photic zone time-series data [g/m2/year]. Default NA. 
#' It will be further calculated based on the Martin curve to get detrital flux reaching the bottom and then multiplied the trophic transfer efficiency (10\%) to get benthic productivity `bprod_ts` ultimately .
#' See source code of \code{setupTimeseries}.\cr
#' Input either of `bprodin_ts`, `dfbot_ts` or `dfpho_ts`. If all are NAs then `p$bprod` will be used in ts simulation `p$r[3] = p$bprod`. In this case, remember to add benthos arguments in `p = setupXX()`.
#' @param Tp_ts Pelagic water temperature, representing the top 100m average temperature [Celsius].
#' @param Tm_ts Mid-water temperature, representing the average temperature of 500m - up to 1500m [Celsius].
#' @param Tb_ts Bottom water (the bottom layer) temperature [Celsius].
#' @param Fsmp_ts Small pelagic fish maximum fishing mortality time-series data [1/year].
#' @param Fmesop_ts Mesopelagic fish maximum fishing mortality time-series data [1/year].
#' @param Flgp_ts Large pelagic fish maximum fishing mortality time-series data [1/year].
#' @param Fmidwp_ts Mid-water predator maximum fishing mortality time-series data [1/year].
#' @param Fdem_ts Demersal fish maximum fishing mortality time-series data [1/year].
#' @param benthosK Carrying capacity of small benthos used for logistic growth [g/m2]. Default is 80. 
#' If not provided, the value will remain the same as the `p$bprod` set in `p = setupXX()`.
#' 
#' @details
#' The setupTimeseries extends the prescribed setup for time-series simulations. It adds time-series data input and parameters related to the time-series simulations.
#' The main operation on each time-series data input is adding an extra element at the end, which is just a replicate of the last element of the input. It is required for time integration by the ode solver.
#' All time-series data arrays should be the same length. 
#' 
#' Zooplankton biomass of each time step are directly provided for fish consumption. Therefore, there are no zooplankton population dynamics in time-series simulations, which do not follow semi-chemostat or logistic growth.
#' Benthos follow the logistic growth, due to the biomass is hard to get. Productivity and carrying capacity are required.
#' 
#' Temperature data will be used for temperature effects computation in each timestep. See source code \code{derivativesFEISTYR}.
#' In Fortran, the temperature effects have been pre-calculated and stored in a large matrix (See source code \code{buildforcings}), which is transmitted to Fortran. The data of each timestep will be called automatically.
#' 
#' Units of the bio-related time-series data should be in year. For example, if the original data is monthly data [g/m2/month], it must be converted to yearly data [g/m2/year] by multiplying by 12. See examples below.
#' If a time-series data is input, it should not contain any NAs.
#' `szbio_ts` and `lzbio_ts` must be provided for consumption by fish.
#' `szprod_ts` and `lzprod_ts` must be provided for restricting consumption by fish.
#' 
#' Fishing mortality will be assigned to every size class of a functional type by \code{setFishing} based on the maximum fishing mortality time-series input (e.g., \code{Fsmp_ts}).
#' 
#' 
#' @return
#' The last element of the output time-series data array is the replicate element of the last element in input time-series data array.
#' Time-series data added:
#' \itemize{
#' \item szbio_ts, small mesozooplankton biomass time-series data.
#' \item lzbio_ts, large mesozooplankton biomass time-series data.
#' \item szprod_ts, small mesozooplankton productivity time-series data.
#' \item lzprod_ts, large mesozooplankton productivity time-series data.
#' \item bprod_ts, benthic productivity time-series data, based on `bprodin_ts` or `dfbot_ts`, or `dfpho_ts` input.
#' }
#' Time-series temperature data, which will be used for temperature effects every time step (see source code of \code{\link{derivativesFEISTYR}}).
#' \itemize{
#' \item Tp_ts, pelagic water temperature time-series data.
#' \item Tm_ts, mid-water temperature time-series data.
#' \item Tb_ts, bottom water temperature time-series data.
#' }
#' 
#' \item K[3], the third element of K is overwritten by the `benthosK` input, representing the carrying capacity of the small benthos community, which follows the logistic growth.
#' Note the other elements in `K` do not effective since the time-series simulation has no zooplankton dynamics.
#' 
#' Fishery-related data:\cr
#' \itemize{
#' \item Fsmp_ts, Small pelagic fish maximum fishing mortality time-series data.
#' \item Fmesop_ts, Mesopelagic fish maximum fishing mortality time-series data.
#' \item Flgp_ts, Large pelagic fish maximum fishing mortality time-series data.
#' \item Fmidwp_ts, Mid-water predator maximum fishing mortality time-series data.
#' \item Fdem_ts, Demersal fish maximum fishing mortality time-series data.
#' }
#' 
#' \item bTS, boolean flag of time-series simulation.
#' 
#' 
#' @examples
#' # Two example time series data for 1850 (1 year) and 1850-2014 (165 year) of one grid. They are monthly data but converted to yearly by multiplying 12, since FEISTY units are in year.
#' # However in simulations, each data only run for 1/12 year to align with month.
#' # Photic zone depth and water column depth are not time-varying so they are in `setupVertical2()`.
#' # In the example data, zooplankton biomass (`Zbio`) and production (`Zhploss`) are halved to represent small and large zooplankton.
#' 
#' # One year data example of one grid.
#' data(tsinput_example_1850)
#' p = setupTimeseries(p = setupVertical2(photic = photic, depth = depth, nStages = 9),
#'                   Tp_ts = Tp,
#'                   Tm_ts = Tm,
#'                   Tb_ts = Tb,
#'                   szbio = Zbio/2,
#'                   lzbio = Zbio/2,
#'                   szprod_ts = Zhploss/2,
#'                   lzprod_ts = Zhploss/2,
#'                   dfbot_ts  = dfbot)
#' # Run by R
#' # `tEnd = 1` represents simulation time is one year. 
#' # `tStep  = 1/12` represents 1/12 year (month), so the data for simulation will update every 1/12 year from the time-series.
#' # If the user has yearly data for 20 years, the arguments should be `tEnd = 20` and `tStep = 1`.
#' simR = simulateFEISTY(p = p,
#'                      tEnd = 1,
#'                      tStep  = 1/12,
#'                      USEdll = F,
#'                      spinup = T)
#' # Run by Fortran
#' simF = simulateFEISTY(p = p,
#'                      tEnd = 1,
#'                      tStep  = 1/12,
#'                      USEdll = T,
#'                      spinup = T)
#' plotBiomasstime(simR)
#' plotBiomasstime(simF)
#'
#' # Example of 165 years (1850-2014).
#' data(tsinput_example_1850_2014)
#' p = setupTimeseries(p = setupVertical2(photic = photic, depth = depth, nStages = 9),
#'                   tStep_ts = 1/12,
#'                   Tp_ts = Tp,
#'                   Tm_ts = Tm,
#'                   Tb_ts = Tb,
#'                   szbio = Zbio/2,
#'                   lzbio = Zbio/2,
#'                   szprod_ts = Zhploss/2,
#'                   lzprod_ts = Zhploss/2,
#'                   dfbot_ts  = dfbot,
#'                   Fsmp_ts   = Fspel,
#'                   Fmesop_ts = Fmeso,
#'                   Flgp_ts   = Flpel,
#'                   Fmidwp_ts = FmidP,
#'                   Fdem_ts   = Fdem)
#' # When simulation time is long, running by Fortran is much faster than by R.
#' sim = simulateFEISTY(p = p,
#'                      tEnd = 165,
#'                      tStep  = 1/12,
#'                      USEdll = T,
#'                      spinup = T)
#' plotBiomasstime(sim)
#' 
#' @seealso
#' \code{\link{setupBasic}}     \cr
#' \code{\link{setupBasic2}} 	  \cr
#' \code{\link{setupVertical2}} \cr
#' 
#' @aliases setupTimeseries
#' @export
#' 
# p = setupTimeseries(p = setupVertical2(photic = photic, depth = depth, nStages = 9),
# tStep_ts = 1/12,
# tSpin = 1.6, #[year]
# nSpinloop = 4, #[]
# Tp_ts = Tp[1:200],
# Tm_ts = Tm[1:200],
# Tb_ts = Tb[1:200],
# szbio = Zbio[1:200]/2,
# lzbio = Zbio[1:200]/2,
# szprod_ts = Zhploss[1:200]/2,
# lzprod_ts = Zhploss[1:200]/2,
# dfbot_ts  = dfbot[1:200])
#
# sim = simulateFEISTY(p = p, USEdll = F)
# sim2 = simulateFEISTY(p = p, USEdll = T)
# p = setupTimeseries(p = setupVertical2(photic = photic, depth = depth, nStages = 9),
#                     tStep_ts = 1/12,
#                     tSpin = 1.6, #[year]
#                     nSpinloop = 4, #[]
#                     Tp_ts = Tp,
#                     Tm_ts = Tm,
#                     Tb_ts = Tb,
#                     szbio = Zbio/2,
#                     lzbio = Zbio/2,
#                     szprod_ts = Zhploss/2,
#                     lzprod_ts = Zhploss/2,
#                     dfbot_ts  = dfbot)
# sim = simulateFEISTY(p = p, USEdll = F)
# sim2 = simulateFEISTY(p = p, USEdll = T)
setupTimeseries = function (p = setupVertical2(),
                            tStep_ts = 1/12, # [year]
                            tSpin = 10, #[year]
                            nSpinloop = 4, #[]
                            szbio_ts = NA,#Zbio/2, #c(1e3,1e3)
                            lzbio_ts = NA,#Zbio/2,
                            szprod_ts = NA,#Zhploss/2,
                            lzprod_ts = NA,#Zhploss/2,
                            bprodin_ts = NA, # benthos production
                            dfbot_ts  = NA,#dfbot,#NA, # detrital flux reaching the bottom
                            dfpho_ts  = NA, # detrital flux out of photic zone
                            Tp_ts = NA, #Tp,
                            Tm_ts = NA, #Tm,
                            Tb_ts = NA, #Tb,
                            Fsmp_ts = NA, # smallPel
                            Fmesop_ts = NA, # mesoPel
                            Flgp_ts = NA, # largePel
                            Fmidwp_ts = NA, # midwPred
                            Fdem_ts = NA, # demersals
                            benthosK = 80){
  p$bTS = TRUE
  p$tSpin = tSpin
  p$nSpinloop = nSpinloop
  args <- list(
    szbio_ts = szbio_ts, lzbio_ts = lzbio_ts, szprod_ts = szprod_ts,
    lzprod_ts = lzprod_ts, bprodin_ts = bprodin_ts, dfbot_ts = dfbot_ts,
    dfpho_ts = dfpho_ts, Tp_ts = Tp_ts, Tm_ts = Tm_ts, Tb_ts = Tb_ts,
    Fsmp_ts = Fsmp_ts, Fmesop_ts = Fmesop_ts, Flgp_ts = Flgp_ts, Fmidwp_ts = Fmidwp_ts, Fdem_ts = Fdem_ts)
  
  # Check which inputs are NOT NA
  not_na_args <- names(args)[!unlist(lapply(args, function(x) identical(x, NA)))]
  # Check if all time-series inputs have same length
  if (length(unique(sapply(args[not_na_args], length))) != 1) stop("All time-series inputs must have same length.")
  
  if (length(not_na_args) > 0) {
    # Check if the given time-series inputs contain NA values
    if ( any(sapply(args[not_na_args], function(x) any(is.na(x)))) ) {
      input_with_na=not_na_args[sapply(args[not_na_args], function(x) any(is.na(x)))]
      stop(paste("Please check time-series inputs. The following inputs contain NA values:", 
                 paste(input_with_na, collapse = ", ")))
    }
    # Check if szbio_ts, lzbio_ts, szprod_ts, and lzprod_ts are provided
    args2 <- list(szbio_ts = szbio_ts, lzbio_ts = lzbio_ts, szprod_ts = szprod_ts,lzprod_ts = lzprod_ts)
    if ( any(sapply(args2, function(x) identical(x, NA))) ) {
      na_args <- names(args2)[unlist(lapply(args2, function(x) identical(x, NA)))]
      stop(paste("The following time-series data must be provided:", paste(na_args, collapse = ", ")))
      }
    # print all valid time-series inputs names
    cat(sprintf("Time-series input: %s. \n", paste(not_na_args, collapse = ", "))) 
  } else {
    stop("No time-series inputs are provided.")
  }
  
  tslength = length(unlist(args[not_na_args[1]])) # length of time-series input
  p$tStep_ts = tStep_ts
  p$tEnd_ts  = tslength/(1/tStep_ts)
  # print time parameters
  cat(sprintf("Time-series length: %s %s. \n", p$tEnd_ts, "years"))
  cat(sprintf("Time-series time step: %s %s. \n", MASS::fractions(p$tStep_ts), "years"))
  p$szbio_ts = szbio_ts #seq(from=100, to=800, length.out=12)
  p$szbio_ts[tslength+1] = szbio_ts[tslength] #p$zbio_ts[13] = 800
  p$lzbio_ts = lzbio_ts
  p$lzbio_ts[tslength+1] = lzbio_ts[tslength] 
  p$szprod_ts = szprod_ts
  p$szprod_ts[tslength+1] = szprod_ts[tslength] 
  p$lzprod_ts = lzprod_ts
  p$lzprod_ts[tslength+1] = lzprod_ts[tslength] 
  p$Tp_ts=Tp_ts
  if (all(!is.na(Tp_ts))) p$Tp_ts[tslength+1] = Tp_ts[tslength] 
  p$Tm_ts=Tm_ts
  if (all(!is.na(Tm_ts))) p$Tm_ts[tslength+1] = Tm_ts[tslength] 
  p$Tb_ts=Tb_ts
  if (all(!is.na(Tb_ts))) p$Tb_ts[tslength+1] = Tb_ts[tslength] 
  
  # benthic production calc
  if (all(is.na(bprodin_ts)) & all(is.na(dfbot_ts)) & all(is.na(dfpho_ts))){
  p$bprod_ts = NA
  } else {
  if (sum(all(!is.na(bprodin_ts)), all(!is.na(dfbot_ts)), all(!is.na(dfpho_ts)))>1) stop('Please check "bprodin_ts", "dfbot_ts" and "dfpho_ts" input. 
                                                                                         Only one of them should be assigned values, others should be kept as "NA".')
  if (all(!is.na(bprodin_ts))) {bprod_ts = bprodin_ts} #else {bprodin_ts = -1}
  if (all(!is.na(dfbot_ts))) {bprod_ts = dfbot_ts*0.1} #else {dfbot_ts = -1}
  if (all(!is.na(dfpho_ts))) {bprod_ts = 0.1*(dfpho_ts*(depth/photic)^-0.86)
                              bprod_ts[bprod_ts >= 0.1*dfpho_ts ] = 0.1*dfpho_ts[bprod_ts >= 0.1 * dfpho_ts]} #else {dfpho_ts = -1}
  #  
  p$bprod_ts=bprod_ts
  p$bprod_ts[tslength+1] = bprod_ts[tslength]
  }
  
  # fishing ts
  if (p$setup %in% c("setupBasic", "setupBasic2")) {
    if (any(!is.na(c(Fmesop_ts, Fmidwp_ts)))) {
      stop("In setupBasic and setupBasic2, mesopelagic fish (Fmesop_ts) and midwater predators (Fmidwp_ts) should not exist.")
    }
  }
  p$Fsmp_ts = Fsmp_ts
  if (all(!is.na(Fsmp_ts))) p$Fsmp_ts[tslength+1] = Fsmp_ts[tslength]
  p$Fmesop_ts = Fmesop_ts
  if (all(!is.na(Fmesop_ts))) p$Fmesop_ts[tslength+1] = Fmesop_ts[tslength]
  p$Flgp_ts = Flgp_ts
  if (all(!is.na(Flgp_ts))) p$Flgp_ts[tslength+1] = Flgp_ts[tslength]
  p$Fmidwp_ts = Fmidwp_ts
  if (all(!is.na(Fmidwp_ts))) p$Fmidwp_ts[tslength+1] = Fmidwp_ts[tslength]
  p$Fdem_ts = Fdem_ts
  if (all(!is.na(Fdem_ts))) p$Fdem_ts[tslength+1] = Fdem_ts[tslength]
  
  #update benthos carrying capacity, benthos biomass cannot beyond this value.
  p$K[3] = if (is.na(benthosK)) 80 else benthosK 
  
  return(p)
}


# ------------------------------------------------------------------------------
# Make a basic setup with just pelagic fish. Currently not functional
#
# Out:
#  An updated parameter list. The feeding preferences are quite complex
#
# ------------------------------------------------------------------------------

# setupPelagicSpecies = function(depth=500, pprod=100, bprod=5, 
#                                nStages=6, mInf=125000, names=NA, demersal=TRUE, mort0=0.5) {
#   
#   # Initialize the parameters:
#   param = paramInit(depth=depth, pprod=pprod, bprod=bprod)
#   
#   # Setup resource groups:
#   param = paramAddResource(
#     param, 
#     names= c("smallZoo", "largeZoo", "benthos", "largeBenthos"),
#     K    = c(pprod, pprod, bprod, 0),  # g ww/m2  - maximum resource concentration
#     r    = c(1, 1, 1, 1),              # [/yr] nudging coefficient
#     mc   = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500)),
#     mLower = c(2e-06,0.001, 0.5e-03, 0.25), # weight lower limit)  
#     mUpper = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500))
#   )
#   
#   
#   # Add fish groups:
#   names = rep(names, length.out=length(mInf))
#   for (iGroup in 1:length(mInf))
#     param = paramAddGroup(param, 
#                           mMin=0.001, mMax=mInf[[iGroup]], 
#                           mMature=0.25*mInf[[iGroup]], 
#                           nStages=nStages, 
#                           name=names[iGroup])
#   
#   # physiology of all fish stages
#   param = paramAddPhysiology(param)
#   
#   # Setup size interaction matrix:
#   param$theta = matrix(nrow=param$nStages, ncol=param$nStages, data=0)
#   rownames(param$theta) <- colnames(param$theta) <- param$stagenames
#   
#   beta = 400   # preferred predator/prey mass ratio
#   sigma = 1.3  # width of size preference for feeding
#   for (i in param$ixFish) {
#     param$theta[i,] = exp( -(log(param$mc[i]/(beta*param$mc)))^2 / (2*sigma)^2  )
#     param$theta[i,param$mc>param$mc[i]] = 0
#   }
#   param$theta[is.na(param$theta)] = 0
#   
#   #
#   # Setup interactions between groups and resources:
#   #
#   
#   mMedium = 10
#   mLarge = 5000
#   
#   ixR = param$ixR
#   
#   if (!demersal)
#     param$theta[,ixR[3:4]] = 0  # No demersal feeding
#   
#   param$mort0 = mort0 # NOTE: set pretty high to give a stable population
#   
#   return(param)
# }
# 
# 
