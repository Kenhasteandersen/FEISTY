\name{simulateFeisty}
\alias{simulateFeisty}
\alias{derivativesFeistyR}
\title{
dynamic application of the Feisty model
}
\description{

FEISTY is a size- and trait-based model of higher trophic level dynamics.

The Fishes in FEISTY are defined by their functional type, size, and maturity stage. 

The number of fish functional types as well as the number of size classes in each functional type can be chosen.

The biological rates in the model are based on allometry using the mass, w (g), of each size class.  

Fish size classes are defined using logarithmic size bins appropriate for the life history stage it represents
by holding the ratio of initial and final body sizes of each size-class, constant across all stages.

The stage-structured formulation is as in De Roos et al. (2008).  

Each fish functional type is modeled as a collection of size-classes, with smaller classes growing into larger size classes. 

The dynamics of the fish resources is modeled in a very simple manner.

The state variables are:


  \itemize{
   \item the \code{fish biomass density} in each size-class, and for each functional type, expressed in gWW/m2.
   \item the \code{biomass density of the food resources}, expressed in gWW/m2
  }

}

\usage{
simulateFeisty(p = setupBasic(), tEnd = 100, times = seq(from = 0, to = tEnd, by = 1), 
               yini = p$u0, USEdll = TRUE, Rmodel = derivativesFeistyR, simpleOutput = FALSE)
               
derivativesFeistyR(t, u,  p,  FullOutput=TRUE)
}
\arguments{
  \item{p }{the model parameters, as created with e.g. \link{setupBasic}
  }
  \item{tEnd, times }{the end time or the time sequence for which the model needs to be run.
  }
  \item{yini }{the initial condition of the resources and fish density biomass; if not provided, the default in the model parameter list will be used. 
  }
  \item{USEdll }{if TRUE, then the fortran version of the model will be used, else the R-version will be used. 
  }
  \item{simpleOutput }{if TRUE, then the output will be a list with resources (first item), the fish stages (second item), the yield and the Spawning Stock Biomass (SSB). 
  }
  \item{t, u }{the time and state variable input as required by the \link[deSolve]{ode} function from the \code{deSolve} package.
  }
  \item{FullOutput }{if TRUE: will also output the ordinary variables; if false only the state variables will be outputted
  } 
  \item{Rmodel }{Name of the R-function that calculates the model derivatives (as defined by the deSolve package) - only used if USEdll=FALSE}
}

\value{

  \code{simulateFeisty} either generates an object of class \code{deSolve} (if simplifyOutput=FALSE), or a list with resources, fish stage biomass, yield and spawning stoch biomass.

If simplifyOutput=FALSE, then the following output variables are generated for each stage (includes resources:fish stages):

\itemize{
 \item  f..,         [-],                          food limitation of each stage class
 \item  mortpred, [gWW/m2/time],          mortality rate due to predation of each stage

}
output variables for each fish stage are:
\itemize{
 \item  g..,         [/time],                       net growth rate of each stage class
 \item  Repro..,    [gWW/m2/time],                  reproduction rate of each stage class
 \item  Fin,  Fout,  [gWW/m2/time],                 input and output of each stage class
}

output variables for each fish functional type:
\itemize{
 \item  totMort, [gWW/m2/time],        total mortality of each fish functional type
 \item  totGrazing, [gWW/m2/time],     total grazing of each fish functional type
 \item  totLoss, [gWW/m2/time],        total loss rates of each fish functional type (respiration+faeces production)
 \item  totRepro, [gWW/m2/time],       total reproduction rates of each fish functional type 
 \item  totRecruit, [gWW/m2/time],     total recruitment into each fish functional type (a fraciton of totRepro)
 \item  totBiomass, [gWW/m2],          total biomass of each fish functional type 
}
  
}
\seealso{ 
\link{parametersFeisty} for a description of the parameter input file.
}

\author{
Ken H Andersen
Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{ 
#-----------------------------------------------
# run model with default parameter settings
#-----------------------------------------------

RR <- simulateFeisty()

colnames(RR)

plot(RR, which = 1:9) # plot first 9 state variables

plot(RR, which=c( "smallZoo",  "largeZoo", "smallBenthos", "largeBenthos",
   "totBiomass.smallPel", "totBiomass.largePel", "totBiomass.Demersals"))
   
matplot.0D(RR, type= "l", lty=1, ylab="g/m2", log="y", main="Large Pelagics",
  which=c("largePel_1", "largePel_2", "largePel_3", "totBiomass.largePel"))

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

STILL A LOT TO BE EXPLAINED HERE


}

\note{
There are two versions of the same model:

For efficiency reasons, the model was written in Fortran.

For model development, an R-version of the model is also available.
}
\keyword{misc}

