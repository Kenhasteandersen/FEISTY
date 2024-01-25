#'
#' 
#' @docType package
#' 
#' @name FEISTY
#' 
#' @title FEISTY: FishErIes Size and functional TYpe model.
#' 
#' @description 
#' FEISTY is a size- and trait-based model of marine higher trophic level dynamics (fish and cephalopods).
#' The model is described in the publications:  Petrik et al., 2019 and van Denderen et al., 2020.
#' This package is developed for FEISTY model simulations and development and includes modifications to 
#' the published setups. It also includes a Shiny app interface.
#' 
#' @details 
#' The package provides four setups for FEISTY, including two published setups and their revised versions.
#' Furthermore, this package allows users to customize their setups to simulate and develop their new ideas within the FEISTY model frame.
#' Functions of parameter setting are offered for customization.
#' 
#' The fish in FEISTY are defined by their functional type, size, and maturity stage. 
#' The number of fish functional types is defined by each of the predefined setups (see below), but can also be adapted by the user.
#' The physiological rates of each size class are based on allometry using the mass (gram wet weight) of each size class. 
#' The fish eat on resources of zooplankton or benthos, or, for larger size classes, on smaller fish. 
#' Fish size classes are defined using logarithmic size bins appropriate for the life history stage; 
#' they are represented by holding the ratio of initial and final body sizes of each size class, constant across all stages.
#' Each fish functional type is modeled as a collection of size classes, with smaller classes growing into larger size classes. 
#' The stage-structured formulation is as in de Roos et al. (2008).
#' 
#' The package provides two simulation approaches: either using a compiled Fortran code (very fast) or running everything in R (slow, but easier to modify).
#' Both rely on the \link{deSolve} package. 
#' 
#' The package operates with two main structures: the parameters structure and the simulation structure. \cr
#' The parameters are set with a call to one of four "setup" functions. The parameters structure is passed to the `simulateFEISTY`, which returns the results in a `sim` structure.
#' This structure contains all the information of the simulation and can be passed to routines for plotting. \cr
#' 
#' The simplest simulation example can be found below.
#' 
#' 
#' ... ... ...
#' 
#' @examples
#' # A simulation can be as simple as:
#'   p = setupBasic()
#'   sim = simulateFEISTY(p=p)
#'   plotSimulation(sim)
#' # To see what 'sim' contains, please check the help page of the function 'simulateFEISTY'.
#' 
#' @references
#' Petrik, C. M., Stock, C. A., Andersen, K. H., van Denderen, P. D., & Watson, J. R. (2019). Bottom-up drivers of global patterns of demersal, forage, and pelagic fishes. Progress in oceanography, 176, 102124.
#' 
#' van Denderen, P. D., Petrik, C. M., Stock, C. A., & Andersen, K. H. (2021). Emergent global biogeography of marine fish food webs. Global Ecology and Biogeography, 30(9), 1822-1834.
#' 
#' de Roos, A. M., Schellekens, T., Van Kooten, T., Van De Wolfshaar, K., Claessen, D., & Persson, L. (2008). Simplifying a physiologically structured population model to a stage-structured biomass model. Theoretical population biology, 73(1), 47-62.
#' 
#' Soetaert, K., Petzoldt, T., & Setzer, R. W. (2010). Solving differential equations in R: package deSolve. Journal of statistical software, 33, 1-25.
#' 
#' @author Ken H Andersen,
#'         Karline Soetaert <karline.soetaert@nioz.nl>,
#'         Yixin Zhao
#' 
#' @seealso 
#' Setups and simulations: \cr
#' \code{\link{setupBasic}} The setup following Petrik et al. (2019) \cr
#' \code{\link{setupBasic2}} A revised setup based on `setupBasic` \cr
#' \code{\link{setupVertical}} The setup following van Denderen et al. (2021) \cr
#' \code{\link{setupVertical2}} A revised setup based on `setupVertical` \cr
#' \code{\link{simulateFEISTY}} The main function to run FEISTY simulations \cr
#' \code{\link{webFEISTY}} A shiny interface for visualizing FEISTY model results
#' 
#' Plotting: \cr
#' \code{\link{plotSimulation}} Plot simulation results \cr
#' \code{\link{plotSSBtime}} Spawning stock biomass plot \cr
#' \code{\link{plotSpectra}} Biomass plot \cr
#' \code{\link{plotRates}} Plots for growth rate, mortality, and feeding level \cr
#' \code{\link{plotNetwork}} Food web plot \cr
#' \code{\link{plotDiet}} Diet plot
#' 
#' Fishing and diagnostics: \cr
#' \code{\link{calcSSB}} 	Spawning stock biomass calculation \cr
#' \code{\link{calcYield}} Yield calculation
#' 
#' Parameters: \cr
#' \code{\link{paramInit}} 	Initialize parameters for FEISTY \cr
#' \code{\link{paramAddResource}} 	Add resource parameters \cr
#' \code{\link{paramAddGroup}} 	Add parameters of one functional type \cr
#' \code{\link{paramAddPhysiology}} 	Add physiological parameters \cr
#' \code{\link{paramSizepref}} 	Size preference matrix calculation \cr
#' \code{\link{paramTeffect}} 	Add temperature effects \cr
#' \code{\link{setFishing}} 	Set fishing mortality \cr
#' 
#' \code{\link{derivativesFEISTYR}} The derivative function of FEISTY
#' 
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics axis box image legend lines mtext par points polygon
#' @importFrom stats na.omit simulate
#' @importFrom pracma erf size linspace
#' @importFrom latex2exp TeX
#' @import deSolve
#' @import shiny
#' @import shinyjs
#' @useDynLib FEISTY, .registration = TRUE
#' 
NULL # This 'NULL' is important in R package documentation file when using roxygen2.