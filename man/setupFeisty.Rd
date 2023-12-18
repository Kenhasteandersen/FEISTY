\name{setupFeisty}
\alias{setupFeisty}
\alias{setupPelagicSpecies}
\alias{setupBasic}
\alias{setupBasic2}
\alias{setupVertical} 
\alias{setupVertical2} 


\title{
  Feisty model parameter sets.
}
\description{

\code{setupPelagicSpecies} sets up model parameters for pelagic species only 

\code{setupVertical} sets up model parameters for pelagic+benthic species, where
pelagic species feed from zooplankton that has a different vertical distribution day and night

\code{setupVertical2} sets up model parameters for pelagic+benthic species, where
pelagic species feed from zooplankton that has a different vertical distribution day and night

\code{setupBasic} creates a basic three-species setup as described in Petrik et al (2019) 
 
\code{setupBasic2} creates a basic three-species setup generalised to 
more realistic sizes, and size-based feeding
}

\usage{
p=setupBasic(szprod = 100,
             lzprod = 100,
             bprod  = 5,
             depth  = 100,
             Tp     = 10,
             Tb     = 8)
                      
p=setupBasic2(szprod = 100,
              lzprod = 100,
              bprod  = 5,
              depth  = 100,
              Tp     = 10,
              Tb     = 8,
              nStages=9,
              etaMature=0.25,
              F=0,
              etaF=0.05)

p=setupVertical(szprod= 80,
                lzprod = 80,
                bent = 150,
                region = 4,
                depth=1500,
                photic=150)

p=setupVertical2(szprod= 80,
                 lzprod = 80,
                 bent = 150,
                 nStages=6,
                 region = 4,
                 depth=1500,
                 photic=150,
                 mesopelagic=250,
                 visual=1.5,
                 etaMature = 0.25, # Size of matureation relative to
                         # asymptotic size. Different from
                         # van Denderen (2021), where it is 0.002
                 F=0,
                 etaF=0.05)


%setupPelagicSpecies(depth = 500, pprod = 100, bprod = 5, nStages = 6,     
%     mInf, names = NA, demersal = TRUE, mort0 = 0.5)
}
\arguments{
  \item{szprod}{small zooplankton productivity\cr Actually this represents small zooplankton carrying capacity [gww/m2] but it will mutiply growth rate \bold{r} which is always 1 [1/yr]. Therefore, it is described as small zooplankton productivity [gww/m2/year]. \code{lzprod} and \code{bprod} are same.}
  \item{lzprod}{large zooplankton productivity}
  \item{bprod}{benthic community productivity}  
  \item{depth}{water column depth [meter]. \cr 
  For \code{setupBasic} and \code{setupBasic2}, \code{depth}>=200 is characterized as deep water, and \code{depth}<200 is characterized as shallow water. \code{depth}=300 and \code{depth}=1000 do not have any different effects on simulations. \cr
  For \code{setupVertical} and \code{setupVertical2}, different \code{depth} values will influence fish vertical overlap and temperature-dependent physiological rates.}
  \item{Tp}{pelagic water temperature representing the top 100m average temperature [Celsius].}
  \item{Tb}{bottom water temperature [Celsius].}
  \item{nStages}{size number of large fish functional groups (e.g., large pelagic fish, demersal fish, and bathypelagic fish). The size number of small fish functional groups (e.g., small  pelagic fish and mesopelagic fish) is \code{round(2/3*nStages)}. Generally, \code{nStages} is multiples of 3 (e.g., \code{nStages = 3, 6, 9, or 12}...).}
  \item{etaMature}{the coefficient determines the fish size \code{mMature} with a 50\% maturity level. \code{mMature = etaMature * mMax},  where \code{mMax} is the largest fish size (boundary) of a fish functional group. See \code{\link{paramAddGroup}}. In van Denderen et al. (2021), it was 0.002.}
    
  \item{region}{Different regions: 1 Tropical, 2 Temperate, 3 Boreal, 4 Default 10 Celsius. It represents the water column temperature profile for three regions. The default is 10 Celcius for the whole water column (\code{region = 4}). The source file is in .../data/tempdata.dat. It is the same dataset used in van Ddenderen et al. (2021).}
  \item{F}{fishing mortality [1/year]. If \code{F} is 0, there is no fishing mortality. If \code{F} is assigned a value greater than 0, fishing mortality will be set by multiplying the fishing selectivity \code{psi} which is based on a S-shape function. See source code of \code{\link{setFishing}}.}
  \item{etaF}{the coefficient determining the fish size \code{mFishing} with 50\% fishing selectivity. See source code of \code{\link{setFishing}}.}
  \item{photic}{photic zone depth [meter]. The value affects the diel vertical migration depth. See source code of \code{setupVertical} or \code{setupVertical2}.}
  \item{mesopelagic}{mesopelagic depth [meter]. The value affects the diel vertical migration depth. See source code of \code{setupVertical} or \code{setupVertical2}.}
  \item{bent}{Detrital flux out of photic zone [gww/m2/year]. \code{bent} will be further calculated based on the Martin curve to get detrital flux reaching to the bottom for driving benthic communities. See source code \code{setupVertical}.}
  \item{visual}{\code{visual=1.5}: visual predator, predation ability is enhanced during the day and decreased in the twilight zone during the day, \code{visual=1}: non-visual predator, predation abilities are equal in day and night. It must \bold{be careful} to assign other values to \code{visual}, or the setup could crash. See source code of \code{setupVertical} or \code{setupVertical2}.}
  
%  \item{demersal }{bb  }
%  \item{mInf }{bb}
%  \item{names }{the names of the functional group - will be used to label the model output
%  }
%  \item{mort0 }{natural mortality [/year]. 
%  }
}

\value{

a list with parameter values and model settings containing:

\itemize{
 \item  ..,         any argument passed with paramInit
 \item  nResources,    the total number of resources
 \item  nGroups,    the total number of fish groups
 \item  nStages, the total number of resource+fish stages
 \item  ix,         for each fish group, the index to the stage class
 \item  ixFish      the index to the stage class for all fish stages
 \item  ixR      the index to the stage class for all resources

}

parameter settings for the resource dynamics are:

\itemize{
 \item  dynamics, Rtype: the type of dynamics of the resources, as text and integer respectively
 \item  K, r,     the resource dynamic parameters
 \item mc, mLower, mUpper: the mean, lower and upper size of the resource
 \item  resources, a data.frame with all the above resource information (redundant)
}
  
parameter settings for the fish dynamics are:

\itemize{
 \item mc, mLower, mUpper: the mean, lower and upper size of the stages
 \item z: the size ratio of the stages
 \item psiMature: the fraction of each stage that is mature
 \item mortF, mort0, Cmax, metabolism, V: the rate parameters of the stages
 \item  fishes, a data.frame with all the above fish information (redundant)
}
}
\seealso{ 
\link{simulateFeisty} for running the model.
}

\author{
Ken H Andersen,
Karline Soetaert <karline.soetaert@nioz.nl>,
Yixin Zhao
}

\examples{ 
#-----------------------------------------------
# Create default settings
#-----------------------------------------------

p1 <- setupBasic()
names(p1)

p1$resources
p1$fishes

RR <- simulateFeisty()

#-----------------------------------------------
# run model with altered parameter settings
#-----------------------------------------------
p <- setupBasic()

# increase fraction that is mature 
p$psiMature[p$psiMature == 0.5] <- 0.75

R2 <- simulateFeisty(p=p)

# show last 50 years of simulation
plot(RR, R2, which=1:12, lty=1, mfrow=c(3,4), subset=time >50)
}

\references{

Van Denderen et al. 2020. Emergent global biogeography of marine fish food webs. 
Global Ecology and Biogeography, DOI: 10.1111/geb.13348 

Petrik, CM, Stock, CA, Andersen, KH, van Denderen, PD, Watson, JR 2019. 
Bottom-up drivers of global patterns of demersal, forage, and pelagic fishes. 
Progress in Oceanography, 176, 102124. DOI: 10.1016/j.pocean.2019.102124

De Roos, A.M., Schellekens, T., Van Kooten, T., Van De Wolfshaar, K., Claessen, D.,
Persson, L., 2008. Simplifying a physiologically structured population model to a
stage-structured biomass model. Theor. Popul Biol. 73, 47--62.
}
\details{
Functions setupBasic, setupBasic2, setupVertical, setupPelagicSpecies create parameter input files.

\itemize{
 \item setupBasic makes a basic three-species setup (small pelagics, large pelagics and demersals) 
 as described in Petrik et al (2019). There are 5 resources (small and large zooplankton and small and large benthos)
 \item setupBasic2 makes a basic three-species setup as described in Petrik et al  but..
 }

}

\keyword{misc}

