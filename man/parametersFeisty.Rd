\name{parametersFeisty}
\alias{parametersFeisty}
\alias{paramInit}
\alias{paramAddResource}
\alias{paramAddGroup}
\alias{paramAddPhysiology}
\alias{paramTeffect}
\alias{paramSizepref}

\title{
Creating parameter sets of the FEISTY model
}
\description{

Functions paramInit, paramAddResource, paramAddGroup are helper functions to 
create the input parameter lists. Instructions telling users how to customize their own setups.
}

\usage{
paramInit(...)

paramAddResource(p, K, r=1, dynamics=c("chemostat", "logistic"),
                 mc, mLower = NA, mUpper = NA, names=NA, u0=NA)
                 
paramAddGroup(p, nStages, mMin, mMax, mMature=mMax, 
              mortF=0, mort0=0.1, u0=1, name=NA)
              
paramAddPhysiology(p, ac = 20, bc = -0.25, am = 4, bm = -0.175,
              ae = 70, be = -0.2, epsRepro = 0.01, epsAssim = 0.7)  

paramTeffect(p, Tref, Q10, Q10m, u=NA)

paramSizepref(p, beta = 400, sigma = 1.3, type = 1)
              
}
\arguments{
  \item{p}{the model parameters, as created with e.g. \link{setupBasic}, and that need to be updated
  }
  \item{K, r }{model parameters for resource dynamics, i.e. maximum resource concentration (K) and resource nudging rate (in case type="chemostat") or resource growth rate ("logistic").
  }
  \item{u0 }{the initial value of the resources (paramAddResource) or fish density biomass (paramAddGroup). [gWW/m2]
  }
  \item{dynamics }{how the resource should be modeled: one of "chemostat" (dR/dt=r*(K-R)) or "logistic" (dR/dt=r*R*(1-R/K)). 
  }
  \item{mc }{geometric mean weight of resources - can be used for estimating food preferences [gWW].
  }
  \item{mLower, mUpper }{lower and upper limits of resource weights [gWW].
  }
  \item{names, name }{the names of the resources or functional groups - will be used to label the model output
  }
  \item{mMin, mMax }{minimal fish weight of the first (mMin) and max weight of last stage (mMax), [gWW].
  }
  \item{mMature }{size at which half of the individuals is mature; if NA then only the last size class is assumed 50 percent mature, [gWW]. 
  }
  \item{mortF }{fishing mortality for all size groups [/year]. Can be revised manually or by function \code{\link{setFishing}} later. 
  }
  \item{mort0 }{natural mortality for all size groups (background mortality) [/year]. defult value 0.1. 
  }
  \item{nStages }{Number of stages of the functional group.
  }
  \item{ac, bc }{max. consumption coefficient  [g^(-n)/yr] and exponent     [-]
  }
  \item{am, bm }{metabolism coefficient        [g^(-p)/yr] and exponent     [-]
  }
  \item{ae, be }{clearance rate coefficient (or encounter slope)  [m2*g^(q-1)/yr] and exponent [-]
  }
  \item{epsRepro }{reproduction efficiency [-]
  }
  \item{epsAssim }{assimilation efficiency [-]
  }
  \item{Tref }{Reference temperature. 10 Celsius, generally cannot be other values, unless users define physiological rates based on other reference temperature.
  }
  \item{Q10 }{Q10 factor for the maximum consumption rate \code{Cmax} and clearance rate \code{V} [-]
  }
  \item{Q10m }{Q10 factor for metabolism rates \code{metabolism} [-]
  }
  \item{u=NA }{To be developed...
  }
  \item{beta }{preferred predator/prey mass ratio
  }
  \item{sigma }{width of size preference for feeding
  }
  \item{type }{1 = normal (?), 2 = errf (vanDenderen et al., 2021), 3 = double integration (Andersen and Visser 2023, appendix A.) 
  }
  \item{... }{any input argument that needs to be stored with the parameter list
  }  
   
}

\value{

a list with parameter values and model settings containing:

\itemize{
 \item  ..,         any argument passed with paramInit
 \item  nResources,    the total number of resources
 \item  nGroups,    the total number of fish functional groups
 \item  nStages, the total number of resource + fish stages
 \item  ix,         for each fish functional type, the index to the stage classes
 \item  ixFish      the index to the stage class for all fish stages
 \item  ixR      the index for all resources

}

parameters for the resource dynamics are:

\itemize{
 \item  dynamics, Rtype: the type of dynamics of the resources, as text and integer respectively
 \item  K, r,     the resource dynamic parameters
 \item mc, mLower, mUpper: the mean, lower and upper size of the resource
 \item  resources, a data.frame with all the above resource information (redundant)
}
  
parameters for the fish dynamics are:

\itemize{
 \item mc, mLower, mUpper: the mean, lower and upper size of the stages
 \item z: the size ratio of the stages
 \item psiMature: the fraction of each stage that is mature
 \item mortF, mort0, Cmax, metabolism, V: the rate parameters of the stages
 \item  fishes, a data.frame with all the above fish information (redundant)
}

parameters for the temperature effects are:

\itemize{
 \item Cmaxsave, Vsave, metabolismsave: same as \code{Cmax}, \code{V}, \code{metabolism} saved for temperature effect calculation.
 \item fT: factor of T effects on \code{V} and \code{Cmax} of pelagic fish.
 \item fT_met: factor of T effects on \code{metabolism}.
 \item fT_dem: factor of T effects on \code{V} and \code{Cmax} of all demersal fish in deep water (depth>=200), and small and medium demersal fish in shallow water (depth<200).
 \item fT_met_dem: factor of T effects on \code{metabolism} of all demersal fish in deep water (depth>=200), and small and medium demersal fish in shallow water (depth<200).
 \item eT: effecitve temperature in shallow water (depth<200) for large demersal fish.
 \item fT_dem_shallow: factor of T effects on \code{V} and \code{Cmax} of large demersal fish in shallow water (depth<200).
 \item fT_met_dem_shallow: factor of T effects on \code{metabolism} of large demersal fish in shallow water (depth<200).
}


}
\seealso{ 
\link{simulateFeisty} for running the model.
}

\author{
Ken H Andersen
Karline Soetaert <karline.soetaert@nioz.nl>
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
Functions setupBasic, setupBasic2, setupVertical, setupPelagicSpecies create parameter sets for FEISTY simulations.

\itemize{
 \item \bold{setupBasic} makes a basic 3-species setup (small pelagics, large pelagics and demersals) 
 as described in Petrik et al (2019). There are 4 resources (small and large zooplankton and small and large benthos)
 \item \bold{setupBasic2} is a revised version based on \bold{setupBasic}
 \item \bold{setupVertical} makes a 5-species setup as described in van Denderen et al., 2021.
 \item \bold{setupVertical2} is a revised version based on \bold{setupVertical}.
 }

Functions paramInit, paramAddResource, paramAddGroup and paramAddPhysiology are helper functions to 
create new input parameter lists.

Function \code{paramInit} should be called first, followed by \code{paramAddResource}, \code{paramAddGroup} and finally \code{paramAddPhysiology}

\itemize{
 \item Function \code{paramInit} initialises some elements in the list and prepares for resources and fish groups to be added. Default palette and functional group names are in it.
 \item Function \code{paramAddResource} adds all resources with resource parameters to the list. 
 \item Function \code{paramAddGroup} adds one fish group to the list. 
 \item Function \code{paramAddPhysiology} adds the physiological parameters to the list. maximum consumption rate [/year], standard metabolism (basal respiration) [/year], clearance rate [m2/gWW/year] 
}

\bold{How to customize your own FEISTY setup} \cr



}

\note{
There are two versions of the same model:

For efficiency reasons, the model was written in Fortran.

For model development, an R-version of the model is also available.
}
\keyword{misc}

