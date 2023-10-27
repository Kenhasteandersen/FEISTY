\name{setupFeisty}
\alias{setupFeisty}
\alias{setupPelagicSpecies} 
\alias{setupVertical} 
\alias{setupBasic} 
\alias{setupBasic2}

\title{
  Feisty model parameter sets.
}
\description{

\code{setupPelagicSpecies} sets up model parameters for pelagic species only 

\code{setupVertical} sets up model parameters for pelagic+benthic species, where
pelagic species feed from zooplankton that has a different vertical distribution day and night

\code{setupBasic} creates a basic three-species setup as described in Petrik et al (2019) 
 
\code{setupBasic2} creates a basic three-species setup generalised to 
more realistic sizes, and size-based feeding
}

\usage{
setupPelagicSpecies(depth = 500, pprod = 100, bprod = 5, nStages = 6,      
     mInf, names = NA, demersal = TRUE, mort0 = 0.5)
setupVertical(pprod = 80, depth = 1500, photic = 150, mesopelagic = 250, 
     bent = 150, visual = 1.5)
setupBasic(pprod = 100, bprod = 5)
setupBasic2(pprod = 100, bprod = 5, nStages = 9)
}
\arguments{
  \item{depth }{bb  }
  \item{pprod }{bb  }
  \item{bprod }{bb  }
  \item{mInf }{bb  }
  \item{demersal }{bb  }
  \item{photic }{bb  }
  \item{mesopelagic }{bb  }
  \item{bent }{bb  }
  \item{visual }{bb  }
  \item{nStages }{bb  }
  \item{names }{the names of the functional group - will be used to label the model output
  }
  \item{mort0 }{natural mortality [/year]. 
  }
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
Ken H Andersen
Karline Soetaert <karline.soetaert@nioz.nl>
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

