#
# Code for calculating carbon fluxes, carbon injection, and carbon sequestration.
#
# Written by Julie Lemoine and Ken H Andersen.
# Sequestration calculations based on code by Andre W Visser
#
#library(FEISTY)
#library(data.table)
#library(pracma)
#library(Matrix)
#library(parallel)
#library(doParallel)
#library(foreach)
#library(ggplot2)
#library(patchwork) # To arrange two plots
#library(maps)

#' @import Matrix
#' @import doParallel
#' @import foreach
#' @import parallel
#' @importFrom pracma trapz
NULL

#
# Calculate the flux from carcasses, fecal pellets,  reproduction wastes, and respiration.
# The fluxes are calculated at the position of each size class in units of
# gWW/m2/year.
#
# In:
#  A simulation object
#
# Out:
#  A simulation object with additional fields:
#   fluxCarcass, fluxFecal,fluxRepro, and fluxRespiration
#  all in units: gWW/m2/year
#
#' @export
calcCarbonFluxes <- function(sim) {
  
  p = sim$p # parameters
  
  # Fluxes from grazing. sim$f and sim$B are time x fish-stage matrices, so
  # stage-specific rates must be applied by column. Plain vector * matrix
  # multiplication recycles column-wise in R and assigns rates to wrong stages.
  grazrate = sweep(sim$f, 2, p$Cmax[p$ixFish], "*") # grazing rate [yr^-1]
  graz = sim$B * grazrate # grazing flux (before assimilation) [gWW.m^-2.yr^-1]
  feces = graz * (1 - p$epsAssim)/2 # feces flux [gWW.m^-2.yr^-1] -> (1 - p$epsAssim)/2 A half of unassimilated food is feces.

  # Waste energy from the total energy invested into reproduction to eggs goes to respiration.
  # Dead eggs from eggs to larvae becomes detritus, so goes to fecal pellets.

  # sim$Repro already includes Fout of last stage for each functional type
  eps_egg   = 0.22 # 1-eps_egg is The fraction of reproductive invested used for respiration
  eps_R     = unique(p$epsRepro) / eps_egg # The fraction of eggs that survives
  rep2resp  = sim$Repro * (1 - eps_egg) # metabolic cost of egg production [gWW.m^-2.yr^-1]
  rep2feces = sim$Repro * eps_egg * (1 - eps_R) # dead eggs sink as feces; eps_egg = 0.22 from Andersen 2019 p47

  # Respiration components, kept separate for diagnostics:
  #   basal = allometric basal metabolism, p$metabolism[ixFish] * B (Petrik 2019: bM=-0.175, basal-only)
  #   SDA   = specific dynamic action, half of unassimilated food, (1-epsAssim)/2 * grazrate * B
  #   repro = metabolic cost of egg production, rep2resp = (1-eps_egg)*Repro
  resp_basal = sweep(sim$B, 2, p$metabolism[p$ixFish], "*") # [time x nFish] gWW.m^-2.yr^-1
  resp_sda   = graz * (1 - p$epsAssim)/2                   # [time x nFish] gWW.m^-2.yr^-1
  resp_repro = rep2resp                                     # [time x nFish] gWW.m^-2.yr^-1
  respiration = resp_basal + resp_sda + resp_repro

  # Get last 40% of timeseries
  etaTime <- 0.4
  ixTime  <- which(sim$t >= ((1 - etaTime) * sim$t[sim$nTime]))

  sim$fluxCarcass     <- colMeans(sweep(sim$B[ixTime,, drop = FALSE], 2, p$mort0[p$ixFish], "*")) # deadfalls flux [gWW.m^-2.yr^-1]
  sim$fluxFecal       <- colMeans(feces[ixTime,]) # total feces flux [gWW.m^-2.yr^-1]
  sim$fluxRepro       <- colMeans(rep2feces[ixTime,]) # feces flux from reproduction waste [gWW.m^-2.yr^-1]
  sim$fluxRespiration <- colMeans(respiration[ixTime,]) # total respiration flux [gWW.m^-2.yr^-1]
  # Decomposition of fluxRespiration into per-stage components (sum to fluxRespiration):
  sim$fluxRespirationBasal <- colMeans(resp_basal[ixTime,]) # basal metabolism only [gWW.m^-2.yr^-1]
  sim$fluxRespirationSDA   <- colMeans(resp_sda[ixTime,])   # SDA from grazing [gWW.m^-2.yr^-1]
  sim$fluxRespirationRepro <- colMeans(resp_repro[ixTime,]) # cost of egg production [gWW.m^-2.yr^-1]
  
  return(sim)
}

#
# Calculate injection fluxes calculated at all z-levels
# Units i gC/m3/year
#
#' @export
calcCarbonInjection = function(sim) {
  rho_gWW_gC = 9 # Gram wet weight per gram carbon
  #
  # Calculate the injection from sinking POC with velocity v:
  #
  calcPOCinject = function(J, v) { # v: sinking speed [m/day]

    Solve_Detritus_Euler <- function(z, alpha, zeta_X, v) {
      n <- length(z)
      dz <- c(diff(z), tail(diff(z), 1))

      DX <- numeric(n)
      DX[1] <- 0  # Boundary conditions at the surface -> no detritus

      for (i in 1:(n - 1)) {
        DX[i + 1] <- (DX[i] + dz[i] * zeta_X[i] / v) / (1 + dz[i] * alpha[i] / v)
      }

      return(DX)
    }

    # Multiply with the vertical probability distributions
    Jday <-  t(pDay) * J * 0.5 # Units: gWW/m^3/year
    Jnight <- t(pNight) * J * 0.5

    # remineralization at each depth
    z <- 0:p$bottom
    alpha <- rep(NA, length(z)) # bacterial degradation
    zeta_X <- rep(0, length(z)) # source = poc production at each depth

    # set a depth-dependent bacterial degradation
    alpha[z <= 100] <- rp
    alpha[z > 100 & z <= 1500] <- rm
    alpha[z > 1500] <- rb

    # sum total poc production
    Jtotal <- colSums(Jday + Jnight) # Units gWW/m3/year

    zeta_X[1:length(Jtotal)] <- Jtotal

    # calculate detritus production at each depth
    DX_Euler <- Solve_Detritus_Euler(z, alpha, zeta_X, v)
    # detritus remineralized
    Jinject = alpha * DX_Euler

    res = list()
    res$depth <- z
    res$Jbottom <- sum(Jtotal) - sum(Jinject) # gWW/m2/yr

    Jinject[length(Jinject)] <- Jinject[length(Jinject)] + res$Jbottom
    #Jinject_fecal[z < p$photic] <- 0 # carbon close to the surface is not sequestered

    res$inject <- Jinject

    return(res)
  }

  #
  # Per-stage version of calcPOCinject: keeps each fish stage as a separate
  # column so we can diagnose which stages contribute to below-euphotic flux.
  # The Euler step is linear in the source, so we can vectorize across stages.
  # Returns a matrix [nDepth x nFish] in gWW/m^3/yr.
  #
  calcPOCinject_perstage = function(J, v) {
    Solve_Detritus_Euler_vec <- function(z, alpha, zeta_X_mat, v) {
      n <- length(z)
      ns <- ncol(zeta_X_mat)
      dz <- c(diff(z), tail(diff(z), 1))
      DX <- matrix(0, n, ns)
      for (i in 1:(n - 1)) {
        DX[i + 1, ] <- (DX[i, ] + dz[i] * zeta_X_mat[i, ] / v) /
                       (1 + dz[i] * alpha[i] / v)
      }
      return(DX)
    }

    # Per-stage vertical source (do NOT sum across stages):
    Jday_mat   <- sweep(pDay,   2, J, "*") * 0.5  # [nDepth_pDay x nFish] gWW/m^3/yr
    Jnight_mat <- sweep(pNight, 2, J, "*") * 0.5
    Jtotal_mat <- Jday_mat + Jnight_mat

    z <- 0:p$bottom
    nz <- length(z)
    ns <- length(J)

    alpha <- rep(NA, nz)
    alpha[z <= 100] <- rp
    alpha[z > 100 & z <= 1500] <- rm
    alpha[z > 1500] <- rb

    zeta_X_mat <- matrix(0, nz, ns)
    src_len <- min(nrow(Jtotal_mat), nz)
    zeta_X_mat[1:src_len, ] <- Jtotal_mat[1:src_len, ]

    DX_mat <- Solve_Detritus_Euler_vec(z, alpha, zeta_X_mat, v)
    inject_mat <- sweep(DX_mat, 1, alpha, "*")  # [nz x ns] gWW/m^3/yr

    # Per-stage Jbottom spike: POC reaching seafloor unremineralized.
    Jbottom_per_stage <- colSums(zeta_X_mat) - colSums(inject_mat)
    inject_mat[nz, ] <- inject_mat[nz, ] + Jbottom_per_stage

    return(list(z = z, inject = inject_mat))
  }
  
  
  # parameters
  alpha0 = 0.25 # maximum bacterial degradation rate (from Pinti) [day^-1]
  Q10r = 2 # Q10 remin [-]
  
  p = sim$p
  # adjust bacterial degradation rate with temperature
  Tr=p$Tp # pelagic water temperature
  rfac = Q10r^((Tr-10)/10)
  rp = rfac*alpha0
  
  Tr=p$Tm # mid-water temperature
  rfac = Q10r^((Tr-10)/10)
  rm = rfac*alpha0
  
  Tr=p$Tb # bottom temperature
  rfac = Q10r^((Tr-10)/10)
  rb = rfac*alpha0
  
  col_indices <- p$ixFish
  #
  # Calculation the injection from respiration:
  #
  pDay = p$depthDay[, col_indices, drop = FALSE]
  pNight = p$depthNight[, col_indices, drop = FALSE]
  respiration = t(pDay+pNight) * sim$fluxRespiration/2 # For each size class

  # Per-stage respiration injection split by component [nDepth_pDay x nFish], gWW/m^3/yr.
  # Respiration is released at the fish's depth (DVM-weighted), no POC sinking.
  resp_basal_mat = sweep(pDay + pNight, 2, sim$fluxRespirationBasal, "*") / 2
  resp_sda_mat   = sweep(pDay + pNight, 2, sim$fluxRespirationSDA,   "*") / 2
  resp_repro_mat = sweep(pDay + pNight, 2, sim$fluxRespirationRepro, "*") / 2

  # poc production
  resFecal   = calcPOCinject( as.vector(sim$fluxFecal),    1000) # fecal pellet sinking speed [m/day]
  resCarcass = calcPOCinject( as.vector(sim$fluxCarcass),  2000) # carcass sinking speed [m/day]
  resRepro   = calcPOCinject( as.vector(sim$fluxRepro),    1000) # reproductive detritus sinking speed [m/day]

  # Per-stage POC injection profiles [nDepth x nFish], gWW/m^3/yr.
  resFecal_stage   = calcPOCinject_perstage( as.vector(sim$fluxFecal),   1000)
  resCarcass_stage = calcPOCinject_perstage( as.vector(sim$fluxCarcass), 2000)
  resRepro_stage   = calcPOCinject_perstage( as.vector(sim$fluxRepro),   1000)

  # Make list with injections as output and convert to carbon units:
  res = list()
  res$z = resFecal$depth
  res$Fecal = resFecal$inject / rho_gWW_gC
  res$Carcass = resCarcass$inject / rho_gWW_gC
  res$Repro = resRepro$inject / rho_gWW_gC
  res$Respiration = colSums( respiration ) / rho_gWW_gC

  res$total = res$Fecal + res$Carcass + res$Repro + res$Respiration

  # Per-stage injection profiles, gC/m^3/yr.  [nDepth x nFish]
  # Sums across stages reproduce res$Fecal, res$Carcass, res$Repro, res$Respiration
  # (modulo a tiny binning rounding for the bottom spike).
  res$Fecal_stage            = resFecal_stage$inject   / rho_gWW_gC
  res$Carcass_stage          = resCarcass_stage$inject / rho_gWW_gC
  res$Repro_stage            = resRepro_stage$inject   / rho_gWW_gC
  res$RespirationBasal_stage = resp_basal_mat          / rho_gWW_gC
  res$RespirationSDA_stage   = resp_sda_mat            / rho_gWW_gC
  res$RespirationRepro_stage = resp_repro_mat          / rho_gWW_gC
  res$Respiration_stage      = res$RespirationBasal_stage +
                               res$RespirationSDA_stage   +
                               res$RespirationRepro_stage

  return(res)
}

#
# Simulate a given position in the global map:

### PERHAPS MOVE TO MAIN FEISTY (including data file)
#
#' @export
simulatePosition = function(setup,
                            lat, lon,
                            Fmax=0, ixGroups=NULL, # Specification of fishing parameters sent to setFishing()
                            nStages=9, tEnd=500, glob=NULL) {
  # Output from COBALT

    pp = getParametersPosition(lat, lon, glob=glob)
    p = setup(szprod  = pp$szprod,
              lzprod  = pp$lzprod,
              dfbot   = pp$dfbot,
              depth   = pp$depth,
              Tp      = pp$Tp,
              Tm      = pp$Tm,
              Tb      = pp$Tb,
              photic  = pp$photic,
              nStages = nStages)
    #
    # Set fishing to: Fmax with a Q10=1.8 at depths < 1000 m,
    # and to 90% less at deeper than 1000 m.
    #
    if (length(ixGroups)>0) {
      # Per-group temperature for Q10=1.8 correction:
      #   1=smallPel(Tp) 2=mesoPel(Tm) 3=largePel(Tp) 4=mwpred(Tm) 5=dem(Tb)
      Tgroup <- c(pp$Tp, pp$Tm, pp$Tp, pp$Tm, pp$Tb)
      Fmax0  <- Fmax  # preserve input value across loop iterations

      for (g in ixGroups) {
        Fmax <- Fmax0 * 1.8^((Tgroup[g] - 15) / 10) # Q10=1.8 correction per group
        if (pp$depth > 1000)
          Fmax <- 0.1 * Fmax
        p <- setFishing(p, Fmax, groupidx = g)
      }
    }
    
    sim = simulateFEISTY(p = p, tEnd=tEnd) 
  return(sim)
}

#' @export
getParametersPosition = function(lat, lon, sFile="data/Cobalt global data.csv", glob=NULL) {
  if (is.null(glob))
    glob <- read.csv(sFile)
  
  if (lon<0)
    lon = 360+lon
  
  ix = which.min( (glob$lat-lat)^2 + (glob$lon-lon)^2 ) # Find the best fitting location

  return( list(
    szprod = glob[ix, "szprod"],        # small zooplankton production
    lzprod = glob[ix, "lzprod"],        # large zooplankton production
    dfbot  = glob[ix, "dfbot"],         # detrital flux reaching the bottom
    photic = glob[ix, "photic"],        # photic zone depth
    depth  = glob[ix, "depth"],         # water column depth
    Tp     = glob[ix, "Tp"],            # pelagic water temperature
    Tm     = glob[ix, "Tm"],            # mid-water temperature
    Tb     = glob[ix, "Tb"]
  ))
}

#
# Function to read the original matlab TM file.
#
# loadMatlabTransportMatrix = function(sFilename="data/CTL.mat") {
#   library("R.matlab")
#   #
#   # Transport matrix:
#   #
#   TM <- readMat(sFilename, sparseMatrixClass="Matrix")
# 
#   ## When loaded the names of every variables were missing in the CTL.mat object
#   # missing variable names in "CTL.mat"
#   var_names <- dimnames(TM)[[1]]
#   TM <- setNames(as.list(TM[,1,1]), var_names)
# 
#   # missing variable names in "msk"
#   msk_names <- dimnames(TM$msk)[[1]]
#   TM$msk <- setNames(as.list(TM$msk[,1,1]), msk_names)
# 
#   # missing variable names in "grid"
#   grid_names <- dimnames(TM$grid)[[1]]
#   TM$grid <- setNames(as.list(TM$grid[,1,1]), grid_names)
# 
#   # missing variable names in "MSKS"
#   MSKS_names <- dimnames(TM$MSKS)[[1]]
#   TM$MSKS <- setNames(as.list(TM$MSKS[,1,1]), MSKS_names)
# 
#   # Transfer only the needed variables:
#   T$TR = TM$TR
#   T$M3d = TM$M3d
#   T$msk = TM$msk
#   T$grid$xt = TM$grid$xt
#   T$grid$yt = TM$grid$yt
#   T$grid$zt = TM$grid$zt
#   T$grid$zw = TM$grid$zw
#   T$grid$dzt = TM$grid$dzt
#   T$grid$DXT3d = TM$grid$DXT3d
#   T$grid$DYT3d = TM$grid$DYT3d
#   T$grid$DZT3d = TM$grid$DZT3d
#   
#
# Calculation of A = TR - Sink:
#
  # m <- nrow(TM$TR)
  # sink <- rep(0,m)
  # sink[1:length(msk$hkeep)] <- 1e10 # a strong sink force (1e10) is attributed on surface cells only
  # SSINK <- sparseMatrix(i = 1:m, j = 1:m, x = sink) # sink vector on the diagonal of the SSINK matrix
  # A <- TM$TR - SSINK # calculation of A matrix
  # T$A = lu(A)
  
#   return(T)
# }

#' @export
loadTransportMatrix = function(sFilename=NULL, bLUdecompose=FALSE) {
  sLUfilename = 'LU factored TM.Rdata'
  
  # Always do LU decomposition if the LU-decomposed version does not exist:
  if (!file.exists(sLUfilename)) bLUdecompose=TRUE
  
  if (bLUdecompose) {
    cat("LU decomposing transport matrix.\nTakes time, but is only done once and then saved on disk for future use.\n")
    # Load the original transport matrix:
    if (is.null(sFilename)) {
      data(CTL)
    } else {
      load(sFilename)
    }
    
    # Calculation of A = TR - Sink:
    m <- nrow(TM$TR)
    sink <- rep(0,m)
    sink[1:length(TM$msk$hkeep)] <- 1e10 # a strong sink force (1e10) is attributed on surface cells only
    SSINK <- sparseMatrix(i = 1:m, j = 1:m, x = sink) # sink vector on the diagonal of the SSINK matrix
    A <- TM$TR - SSINK # calculation of A matrix
    TM$A = lu(A)
    
    # Save the LU-decomposed version:
    cat("Saving LU decomposed matrix on disk.\n")
    save(TM, file=sLUfilename, compression_level=9)
  }
  else
    load(sLUfilename)

  return(TM)
}

# simplifyTransportMatrix = function(TM) {
#   T$TR = TM$TR
#   T$M3d = TM$M3d
#   T$msk = TM$msk
#   T$grid$xt = TM$grid$xt
#   T$grid$yt = TM$grid$yt
#   T$grid$zt = TM$grid$zt
#   T$grid$zw = TM$grid$zw
#   T$grid$dzt = TM$grid$dzt
#   T$grid$DXT3d = TM$grid$DXT3d
#   T$grid$DYT3d = TM$grid$DYT3d
#   T$grid$DZT3d = TM$grid$DZT3d
#   
#   return(T)
# }

#
# Project the FEISTY/COBALT injection profile onto the OCIM TM grid.
#
# Output: list(inject = numeric(nz), ix = list(y, x))
#   inject[j] = column-integrated injection in OCIM layer j  (gC m^-2 yr^-1)
#
# FEISTY's discrete mass conservation is SUM-based: inject$total[i] is the
# mass injected within the 1 m bin ending at z[i], and the bottom bin
# (i = length(z)) carries an additional spike Jbottom for POC that reaches
# the seafloor without remineralizing in the water column. Trapezoidal
# integration would weight the bottom bin by 0.5 and silently drop ~half of
# Jbottom (a sizeable fraction of the total flux). We use right-closed
# sum-based binning instead, so every FEISTY meter — including the
# bottom-spike point at z = p$bottom — is counted exactly once.
#
# Bathymetry-mismatch repair: when M3d_col (the per-(i,j) OCIM wet/dry
# column) is supplied, any FEISTY mass that lands in a dry OCIM layer
# (below the OCIM seafloor at this lat/lon, or in a sandwiched dry layer)
# is redirected to the deepest wet OCIM layer. Mass-conserving.
#
project_injection_to_TM <- function(inject, lat, lon, grid, M3d_col = NULL) {
  nz <- length(grid$zt)
  integral <- numeric(nz)
  ix <- list(y = which.min((lat - grid$yt)^2),
             x = which.min((lon - grid$xt)^2))

  z <- inject$z
  f <- inject$total
  ok <- is.finite(z) & is.finite(f)
  z <- z[ok]; f <- f[ok]
  if (!length(z)) return(list(inject = integral, ix = ix))

  # SUM-based binning (NOT trapz): FEISTY stores Jbottom — the POC that
  # reached the seafloor without remineralizing — as a spike in the bottom
  # bin of inject$total. Trapezoidal integration would weight that endpoint
  # by 0.5 and drop ~half of it. Per-meter sum keeps every meter fully.
  dz    <- if (length(z) > 1) c(diff(z), tail(diff(z), 1)) else 1
  edges <- c(grid$zw, grid$zw[nz] + grid$dzt[nz])
  bin   <- findInterval(z, edges, rightmost.closed = TRUE, left.open = TRUE)
  bin   <- pmin(pmax(bin, 1L), nz)

  s <- tapply(f * dz, bin, sum)
  integral[as.integer(names(s))] <- as.numeric(s)

  # Bathymetry-mismatch repair: redirect mass landing in dry OCIM layers
  # (below the OCIM seafloor at this lat/lon, or any sandwiched dry layer)
  # to the deepest wet OCIM layer. Mass-conserving.
  if (!is.null(M3d_col)) {
    wet <- which(M3d_col == 1)
    if (!length(wet)) return(list(inject = numeric(nz), ix = ix))
    dry <- setdiff(seq_len(nz), wet)
    integral[max(wet)] <- integral[max(wet)] + sum(integral[dry])
    integral[dry] <- 0
  }

  list(inject = integral, ix = ix)
}

calc_per_area_sum = function(grid, matrix, depthUpper=0) {
  ix = grid$zt>depthUpper
  dz = grid$DZT3d
  dz[,,!ix] = 0
  return( apply(replace(matrix, is.na(matrix), 0) * dz, c(1,2), sum) )
}

#
# Calculate the amount of carbon sequestered and the sequstration time
#
#' @export
calcCarbonSequestration <- function(TM,  # Transport matrix
                                    matrixInject,  # injection matrix (lat, lon, depth) with same dimensions as TM$grid$M3d
                                    photic = 200   # euphotic zone depth [m]: scalar (default 200m)
                                                   # or 2D matrix (lat x lon) for spatially varying depth from Cobalt data
){
  
  project_to_TM = function(vector) {
    return( replace(array(NA, dim = dim(TM$M3d)), TM$msk$pkeep, vector) )
  }
  
  # Initialization of result list
  result <- NULL
  
  #
  # Setup grid from TM:
  #
  M3d <- TM$M3d                             # 3D array (lat x lon x depth) containing :
                                            #   1 = ocean
                                            #   0 = land
  grid <- TM$grid                           # grid metrics with coordinates and depth
  msk <- TM$msk                             # cells of interest masks
                                            #   hkeep = surface cells
                                            #   pkeep = ocean cells
                                            #   ckeep = interior ocean cells (ckeep = pkeep - hkeep)
  VOL = grid$DXT3d*grid$DYT3d*grid$DZT3d    # volume for each grid cell [m^3]
  V = VOL[msk$pkeep]                        # volume for each ocean grid cell in the transport matrix
  #
  # Calculation of A = TR - Sink:
  #
  #m <- nrow(TM$TR)
  #sink <- rep(0,m)
  #sink[1:length(msk$hkeep)] <- 1e10 # a strong sink force (1e10) is attributed on surface cells only
  #SSINK <- sparseMatrix(i = 1:m, j = 1:m, x = sink) # sink vector on the diagonal of the SSINK matrix
  #A <- TM$TR - SSINK # calculation of A matrix
  #
  # Injection
  #
  dz = grid$dzt  #  thickness of each layers
  Q <- array(0, dim = dim(M3d))   # latitude x longitude x depth
  
  index <- which(matrixInject > 0, arr.ind = TRUE)
  latindex <- index[, 1]  # latitude indices
  lonindex <- index[, 2]  # longitude indices
  depthindex <- index[, 3] # depth indices
  
  for(i in seq_along(latindex)) {
    lati <- latindex[i]
    loni <- lonindex[i]
    depthi <- depthindex[i]
    
    Q[lati, loni, depthi] <- matrixInject[lati, loni, depthi] / dz[depthi] # gC/m3/yr
  }

  # Carbon flux in the ocean cells [gC/m3/yr] -> a vector
  q_ocim <- Q[msk$pkeep]
  q_ocim[is.na(q_ocim)] <- 0
  
  # Carbon sequestration in each grid cell [gC/m3] -> a vector
  cseq <- Matrix::solve(TM$A, -q_ocim, sparse=TRUE)
  result$cseq <- cseq
  
  #
  # Calculate quantities in total, per area, and on the TM grid:
  #
  result$lon = TM$grid$xt
  result$lat = TM$grid$yt
  result$depth = TM$grid$zt
  result$inject_pr_vol = Q
  result$inject_per_area = matrixInject # Injection in each cell (gC/yr/m2)
  
  # Total carbon injection [PgC/yr]
  TotInject <- crossprod(V, q_ocim) / 1e15
  result$TotInject <- TotInject
  
  # Carbon injected below euphotic zone (gC/m2/yr):
  zt <- TM$grid$zt
  if (length(photic) == 1) {
    ixBelowEuphotic = zt > photic
    result$inject_below_euphotic = apply( result$inject_per_area[,,ixBelowEuphotic], c(1,2), sum)
  } else {
    nlat <- dim(matrixInject)[1]; nlon <- dim(matrixInject)[2]
    result$inject_below_euphotic <- matrix(0, nrow = nlat, ncol = nlon)
    for (i in 1:nlat)
      for (j in 1:nlon) {
        ix <- zt > photic[i, j]
        if (any(ix))
          result$inject_below_euphotic[i, j] <- sum(result$inject_per_area[i, j, ix])
      }
  }
  result$TotInject_below_euphotic = sum( result$inject_below_euphotic*grid$Areat ) / 1e15 #PgC/yr
  # Carbon sequestered (steady-state stock from solving A*c = -q):
  result$Cseq          = project_to_TM( cseq )                     # (gC/m3, stock)
  result$Cseq_per_area = calc_per_area_sum( TM$grid, result$Cseq ) # (gC/m2, depth-integrated stock)
  
  # Total carbon sequestered in the ocean [PgC]
  TotSeq <- crossprod(V, cseq) / 1e15
  result$TotSeq <- TotSeq
  
  # Sequestration time [year] on the TM grid
  local_export <- q_ocim
  local_export[local_export == 0] <- NA
  SeqTime <- cseq / local_export
  result$SeqTime <- project_to_TM( SeqTime )
  
  # Total sequestration time [year]
  # Defined as sequestered stock divided by below-euphotic injection,
  # following Pinti et al. (2023, Biogeosciences). Using full-column
  # TotInject would understate the timescale because carbon released
  # within the euphotic zone is removed by the surface sink almost
  # immediately and should not be counted in the denominator.
  TotSeqTime <- if (result$TotInject_below_euphotic > 0)
    TotSeq / result$TotInject_below_euphotic else NA
  result$TotSeqTime <- TotSeqTime
  
  #  df_long <- as.data.frame(cc) %>%
  #    setNames(unique(Cseq_JPOC$lon)) %>%
  #    mutate(lat = unique(Cseq_JPOC$lat)) %>%
  #    pivot_longer(cols = -lat, names_to = "lon", values_to = "cseq") %>%
  #    mutate(lon = as.numeric(lon), cseq = na_if(cseq, 0))
  
  return(result)
}

#
# Calculate carbon sequestration at a range of positions:
#
#' @export
calcGlobalCarbonSequestration = function(TM=loadTransportMatrix(),
                                         lon=c(0,360), lat=c(-90,90),  # Which latitudes to simulate over
                                         Fmax=0, ixGroups=NULL,        # Specification of fishing (set to setFishing())
                                         bPrintStatus=TRUE,
                                         nCores=detectCores()-2,
                                         tEnd=500) {

  tTotal = proc.time()

  # Initialize a matrix with all injections
  matrixInject = array(dim=dim(TM$M3d), data=0)
  
  # Make indices for the lat/lon range:
  lon[ lon<0 ] = 360 + lon[ lon<0 ]
  
  if (length(lon)==1) {
    ix_lon = which.min( (TM$grid$xt-lon)^2 ) # Find the best fitting location
  } else {
    ix_lon = which( TM$grid$xt>=lon[1] & TM$grid$xt<=lon[2] )
  }
  
  if (length(lat)==1) {
    ix_lat = which.min( (TM$grid$yt-lat)^2 ) # Find the best fitting location
  } else {
    ix_lat = which( TM$grid$yt>=lat[1] & TM$grid$yt<=lat[2] )
  }
  grid_idx <- expand.grid(i = ix_lat, j = ix_lon)

  # Filter out land points using the ocean mask (surface layer of M3d):
  is_ocean <- sapply(1:nrow(grid_idx), function(k) TM$M3d[grid_idx$i[k], grid_idx$j[k], 1] == 1)
  grid_idx <- grid_idx[is_ocean, ]

  grid = TM$grid
  glob = read.csv("data/Cobalt global data.csv")  # Read once, pass to all workers

  # Per-(i,j) M3d column for bathymetry-mismatch repair in project_injection_to_TM.
  m3d_cols <- vector("list", nrow(grid_idx))
  for (k in seq_len(nrow(grid_idx)))
    m3d_cols[[k]] <- TM$M3d[ grid_idx$i[k], grid_idx$j[k], ]

  # Setup parallel backend:
  cl <- makeCluster( nCores )
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)

  # Loop over ocean grid points only:
  cat("Simulating FEISTY to calculate injections at",dim(grid_idx)[1], "ocean position(s).\n")
  tStart = proc.time()

  injectTM = foreach(i = 1:dim(grid_idx)[1],
                     .packages = c("FEISTY","pracma"),
                     .verbose = FALSE) %dopar%
    {
      sim = simulatePosition(setupVertical2,
                             grid$yt[ grid_idx$i[i] ],
                             grid$xt[ grid_idx$j[i] ],
                             Fmax=Fmax, ixGroups=ixGroups,
                             tEnd=tEnd, glob=glob)

      # Calculate carbon fluxes at the position of the fish:
      sim = calcCarbonFluxes(sim)
      #totalFlux = sim$fluxCarcass + sim$fluxFecal + sim$fluxRepro + sim$fluxRespiration

      # Calculate the injection
      inject = calcCarbonInjection(sim)

      # Calculate injection on TM grid:
      injectTM = project_injection_to_TM(inject,
                                         grid$yt[ grid_idx$i[i] ],
                                         grid$xt[ grid_idx$j[i] ], grid,
                                         M3d_col = m3d_cols[[i]])

      ix = sim$t >= (1 - 0.4) * max(sim$t) # Last 40% of timeseries, consistent with calcCarbonFluxes

      # Per-stage, per-pathway injection below the local euphotic zone.
      # Units: gC/m^2/yr per stage (integrated over depth on FEISTY's 1 m grid).
      ix_below <- which(inject$z > sim$p$photic)
      colsum_below <- function(M) {
        if (length(ix_below)) colSums(M[ix_below, , drop = FALSE])
        else rep(0, ncol(M))
      }
      inj_be_Fecal_stage         <- colsum_below(inject$Fecal_stage)
      inj_be_Carcass_stage       <- colsum_below(inject$Carcass_stage)
      inj_be_Repro_stage         <- colsum_below(inject$Repro_stage)
      inj_be_RespBasal_stage     <- colsum_below(inject$RespirationBasal_stage)
      inj_be_RespSDA_stage       <- colsum_below(inject$RespirationSDA_stage)
      inj_be_RespRepro_stage     <- colsum_below(inject$RespirationRepro_stage)

      #injectTM$inject
      list( inject=injectTM$inject,
            Biomass=colMeans(sim$totBiomass[ix,]),
            SSB=colMeans(sim$SSB[ix,]),
            Y=colMeans(sim$yield[ix,]),
            B_per_stage=colMeans(sim$B[ix,]),
            fluxFecal_stage              = sim$fluxFecal,
            fluxCarcass_stage            = sim$fluxCarcass,
            fluxRepro_stage              = sim$fluxRepro,
            fluxRespiration_stage        = sim$fluxRespiration,
            # Per-stage source-flux respiration components (gWW/m^2/yr)
            fluxRespirationBasal_stage   = sim$fluxRespirationBasal,
            fluxRespirationSDA_stage     = sim$fluxRespirationSDA,
            fluxRespirationRepro_stage   = sim$fluxRespirationRepro,
            # Per-stage below-euphotic injection per pathway (gC/m^2/yr)
            inj_be_Fecal_stage           = inj_be_Fecal_stage,
            inj_be_Carcass_stage         = inj_be_Carcass_stage,
            inj_be_Repro_stage           = inj_be_Repro_stage,
            inj_be_RespBasal_stage       = inj_be_RespBasal_stage,
            inj_be_RespSDA_stage         = inj_be_RespSDA_stage,
            inj_be_RespRepro_stage       = inj_be_RespRepro_stage,
            photic=sim$p$photic )
    }
  tSim = (proc.time() - tStart)[3]
  cat("FEISTY simulations completed in", round(tSim, 1), "seconds.\n")

  # Put into the injection matrix:
  Biomass = array(data=0, c(dim(matrixInject)[1:2], 5))
  SSB = array(data=0, c(dim(matrixInject)[1:2], 5))
  Yield = array(data=0, c(dim(matrixInject)[1:2], 5))
  nStages_fish = length(injectTM[[1]]$B_per_stage)
  B_per_stage = array(data=0, c(dim(matrixInject)[1:2], nStages_fish))
  # Per-stage carbon fluxes per pathway [lat × lon × stage], units gWW/m^2/yr
  fluxFecal_stage              = array(data=0, c(dim(matrixInject)[1:2], nStages_fish))
  fluxCarcass_stage            = array(data=0, c(dim(matrixInject)[1:2], nStages_fish))
  fluxRepro_stage              = array(data=0, c(dim(matrixInject)[1:2], nStages_fish))
  fluxRespiration_stage        = array(data=0, c(dim(matrixInject)[1:2], nStages_fish))
  fluxRespirationBasal_stage   = array(data=0, c(dim(matrixInject)[1:2], nStages_fish))
  fluxRespirationSDA_stage     = array(data=0, c(dim(matrixInject)[1:2], nStages_fish))
  fluxRespirationRepro_stage   = array(data=0, c(dim(matrixInject)[1:2], nStages_fish))
  # Per-stage below-euphotic injection per pathway [lat × lon × stage], units gC/m^2/yr
  inj_be_Fecal_stage           = array(data=0, c(dim(matrixInject)[1:2], nStages_fish))
  inj_be_Carcass_stage         = array(data=0, c(dim(matrixInject)[1:2], nStages_fish))
  inj_be_Repro_stage           = array(data=0, c(dim(matrixInject)[1:2], nStages_fish))
  inj_be_RespBasal_stage       = array(data=0, c(dim(matrixInject)[1:2], nStages_fish))
  inj_be_RespSDA_stage         = array(data=0, c(dim(matrixInject)[1:2], nStages_fish))
  inj_be_RespRepro_stage       = array(data=0, c(dim(matrixInject)[1:2], nStages_fish))
  photic_map = matrix(200, nrow = length(TM$grid$yt), ncol = length(TM$grid$xt))
  for (i in seq_len(dim(grid_idx)[1])) {
    matrixInject[ grid_idx$i[i], grid_idx$j[i],] = injectTM[[i]]$inject
    Biomass[ grid_idx$i[i], grid_idx$j[i],] = injectTM[[i]]$Biomass
    SSB[ grid_idx$i[i], grid_idx$j[i],] = injectTM[[i]]$SSB
    Yield[grid_idx$i[i], grid_idx$j[i],] = injectTM[[i]]$Y
    B_per_stage[ grid_idx$i[i], grid_idx$j[i],] = injectTM[[i]]$B_per_stage
    fluxFecal_stage[             grid_idx$i[i], grid_idx$j[i],] = injectTM[[i]]$fluxFecal_stage
    fluxCarcass_stage[           grid_idx$i[i], grid_idx$j[i],] = injectTM[[i]]$fluxCarcass_stage
    fluxRepro_stage[             grid_idx$i[i], grid_idx$j[i],] = injectTM[[i]]$fluxRepro_stage
    fluxRespiration_stage[       grid_idx$i[i], grid_idx$j[i],] = injectTM[[i]]$fluxRespiration_stage
    fluxRespirationBasal_stage[  grid_idx$i[i], grid_idx$j[i],] = injectTM[[i]]$fluxRespirationBasal_stage
    fluxRespirationSDA_stage[    grid_idx$i[i], grid_idx$j[i],] = injectTM[[i]]$fluxRespirationSDA_stage
    fluxRespirationRepro_stage[  grid_idx$i[i], grid_idx$j[i],] = injectTM[[i]]$fluxRespirationRepro_stage
    inj_be_Fecal_stage[          grid_idx$i[i], grid_idx$j[i],] = injectTM[[i]]$inj_be_Fecal_stage
    inj_be_Carcass_stage[        grid_idx$i[i], grid_idx$j[i],] = injectTM[[i]]$inj_be_Carcass_stage
    inj_be_Repro_stage[          grid_idx$i[i], grid_idx$j[i],] = injectTM[[i]]$inj_be_Repro_stage
    inj_be_RespBasal_stage[      grid_idx$i[i], grid_idx$j[i],] = injectTM[[i]]$inj_be_RespBasal_stage
    inj_be_RespSDA_stage[        grid_idx$i[i], grid_idx$j[i],] = injectTM[[i]]$inj_be_RespSDA_stage
    inj_be_RespRepro_stage[      grid_idx$i[i], grid_idx$j[i],] = injectTM[[i]]$inj_be_RespRepro_stage
    photic_map[ grid_idx$i[i], grid_idx$j[i] ] = injectTM[[i]]$photic
  }

  # Solve the transport matrix to get sequestration etc.:
  cat("Calculating carbon sequestration...\n")
  tStart = proc.time()
  sequestration = calcCarbonSequestration(TM, matrixInject, photic = photic_map)
  tSeq = (proc.time() - tStart)[3]
  cat("Carbon sequestration completed in", round(tSeq, 1), "seconds.\n")
  #
  # Add results from simulations:
  #
  sequestration$matrixInject = matrixInject
  sequestration$Biomass = Biomass
  sequestration$SSB = SSB
  sequestration$Yield = Yield
  sequestration$B_per_stage = B_per_stage
  sequestration$fluxFecal_stage              = fluxFecal_stage
  sequestration$fluxCarcass_stage            = fluxCarcass_stage
  sequestration$fluxRepro_stage              = fluxRepro_stage
  sequestration$fluxRespiration_stage        = fluxRespiration_stage
  sequestration$fluxRespirationBasal_stage   = fluxRespirationBasal_stage
  sequestration$fluxRespirationSDA_stage     = fluxRespirationSDA_stage
  sequestration$fluxRespirationRepro_stage   = fluxRespirationRepro_stage
  sequestration$inj_be_Fecal_stage           = inj_be_Fecal_stage
  sequestration$inj_be_Carcass_stage         = inj_be_Carcass_stage
  sequestration$inj_be_Repro_stage           = inj_be_Repro_stage
  sequestration$inj_be_RespBasal_stage       = inj_be_RespBasal_stage
  sequestration$inj_be_RespSDA_stage         = inj_be_RespSDA_stage
  sequestration$inj_be_RespRepro_stage       = inj_be_RespRepro_stage

  # Calculate the per-area sequestration for the cells which are simulated:
  area = 0 # Area of all simulated cells
  for (i in seq_len(dim(grid_idx)[1]))
    area = area + TM$grid$Areat[grid_idx$i[i], grid_idx$j[i]]
  sequestration$TotSeq_per_area = sequestration$TotSeq / area * 1e15 # gC/m2
  #
  # Calc total biomass, spawning stock biomass, and yield:
  #
  for (i in 1:dim(sequestration$SSB)[3]) {
    sequestration$TotBiomass[i] = sum( TM$grid$Areat*sequestration$Biomass[,,i], na.rm=TRUE ) / 1e15 # Pg_WW
    sequestration$TotSSB[i] = sum( TM$grid$Areat*sequestration$SSB[,,i], na.rm=TRUE ) / 1e15 # Pg_WW
    sequestration$TotYield[i] = sum( TM$grid$Areat*sequestration$Yield[,,i], na.rm=TRUE ) / 1e15 # Pg_WW/yr
  }
  # Per-stage total biomass [Pg_WW] — one entry per fish size class
  sequestration$TotB_per_stage = sapply(seq_len(dim(sequestration$B_per_stage)[3]), function(s)
    sum( TM$grid$Areat*sequestration$B_per_stage[,,s], na.rm=TRUE ) / 1e15)
  # Per-stage global carbon fluxes per pathway [PgWW/yr]
  sequestration$TotFluxFecal_stage              = sapply(seq_len(nStages_fish), function(s)
    sum( TM$grid$Areat*sequestration$fluxFecal_stage[,,s],              na.rm=TRUE ) / 1e15)
  sequestration$TotFluxCarcass_stage            = sapply(seq_len(nStages_fish), function(s)
    sum( TM$grid$Areat*sequestration$fluxCarcass_stage[,,s],            na.rm=TRUE ) / 1e15)
  sequestration$TotFluxRepro_stage              = sapply(seq_len(nStages_fish), function(s)
    sum( TM$grid$Areat*sequestration$fluxRepro_stage[,,s],              na.rm=TRUE ) / 1e15)
  sequestration$TotFluxRespiration_stage        = sapply(seq_len(nStages_fish), function(s)
    sum( TM$grid$Areat*sequestration$fluxRespiration_stage[,,s],        na.rm=TRUE ) / 1e15)
  sequestration$TotFluxRespirationBasal_stage   = sapply(seq_len(nStages_fish), function(s)
    sum( TM$grid$Areat*sequestration$fluxRespirationBasal_stage[,,s],   na.rm=TRUE ) / 1e15)
  sequestration$TotFluxRespirationSDA_stage     = sapply(seq_len(nStages_fish), function(s)
    sum( TM$grid$Areat*sequestration$fluxRespirationSDA_stage[,,s],     na.rm=TRUE ) / 1e15)
  sequestration$TotFluxRespirationRepro_stage   = sapply(seq_len(nStages_fish), function(s)
    sum( TM$grid$Areat*sequestration$fluxRespirationRepro_stage[,,s],   na.rm=TRUE ) / 1e15)
  # Per-stage global below-euphotic injection per pathway [PgC/yr]
  sequestration$TotInjBE_Fecal_stage            = sapply(seq_len(nStages_fish), function(s)
    sum( TM$grid$Areat*sequestration$inj_be_Fecal_stage[,,s],           na.rm=TRUE ) / 1e15)
  sequestration$TotInjBE_Carcass_stage          = sapply(seq_len(nStages_fish), function(s)
    sum( TM$grid$Areat*sequestration$inj_be_Carcass_stage[,,s],         na.rm=TRUE ) / 1e15)
  sequestration$TotInjBE_Repro_stage            = sapply(seq_len(nStages_fish), function(s)
    sum( TM$grid$Areat*sequestration$inj_be_Repro_stage[,,s],           na.rm=TRUE ) / 1e15)
  sequestration$TotInjBE_RespBasal_stage        = sapply(seq_len(nStages_fish), function(s)
    sum( TM$grid$Areat*sequestration$inj_be_RespBasal_stage[,,s],       na.rm=TRUE ) / 1e15)
  sequestration$TotInjBE_RespSDA_stage          = sapply(seq_len(nStages_fish), function(s)
    sum( TM$grid$Areat*sequestration$inj_be_RespSDA_stage[,,s],         na.rm=TRUE ) / 1e15)
  sequestration$TotInjBE_RespRepro_stage        = sapply(seq_len(nStages_fish), function(s)
    sum( TM$grid$Areat*sequestration$inj_be_RespRepro_stage[,,s],       na.rm=TRUE ) / 1e15)
  
  sequestration$ix_lat = ix_lat
  sequestration$ix_lon = ix_lon
  #
  # Print summary and make plots:
  #
  if (bPrintStatus) {
    cat( c("Total fish biomass: ", format(sum(sequestration$TotBiomass),digits=3), "PgWW \n") )
    cat( c("Total spawning stock biomass: ", format(sum(sequestration$TotSSB),digits=3), "PgWW \n") )
    cat( c("Total fish yield: ", format(sum(sequestration$TotYield),digits=3), "PgWW/yr \n") )
    cat( c("Total carbon injected: ", format(sequestration$TotInject,digits=3), "PgC/yr \n"))
    cat( c("Total carbon injected below euphotic: ", format(sequestration$TotInject_below_euphotic,digits=3), "PgC/yr \n"))
    cat( c("Total carbon sequestered: ", format(sequestration$TotSeq,digits=3), 'PgC \n') )
    cat( c("Average sequestered per area: ", format(sequestration$TotSeq_per_area,digits=3), 'gC/m2 \n') )
    cat( c("Average sequestration time: ", format(sequestration$TotSeqTime,digits=3), 'yr \n') )
    
    p1= plotGlobal( sequestration$lon, sequestration$lat, 
                    sequestration$inject_below_euphotic,
                    sTitle="Injection below euphotic zone", "gC/m2/yr")

    p2 = plotGlobal( sequestration$lon, sequestration$lat, sequestration$Cseq_per_area, 
                     sTitle="Carbon sequestered", "gC/m2")
    
    #p3 = plotGlobal( sequestration$lon, sequestration$lat,
    #                 sequestration$SeqTime[,,6],
    #                 sTitle="Sequestration time at 246 m", "yr")
    
    combined <- plot_grid(p1,p2,ncol=2)
    print(combined)
  }
  
  cat("Total time:", round((proc.time() - tTotal)[3], 1), "seconds.\n")
  return( sequestration )
}

#' @export
plotGlobal = function(lon, lat, data, sTitle="", sUnits="") {
  # Fix range of longitudes if needed:
  ix = which(lon>180)
  lon[ix] = lon[ix]-360
  
  # Transform the data into a data frame:
  df <- expand.grid(lon = lon, lat = lat)
  df$value <- as.vector(t(data))
  
  # Plot:
  ggplot(
    df, 
    aes(x = lon, y = lat, fill = value)) +
    geom_tile() +
    annotation_borders("world", colour = "black") +
    scale_fill_viridis_c() +
    coord_fixed(ratio = 1.3) +
    theme_minimal() +
    labs(x="Longitude", y="Latitude", fill=sUnits, title=sTitle)
  
}

#
# Test carbon calculations at a single position
#
#' @export
testCarbonCalculations_one_position = function(lat=60, lon=-15) {
  
  calcGlobalCarbonSequestration(lon=lon,lat=lat)
  
  # Simulate the position 60, -15 using Cobalt output:
  # sim = simulatePosition(setupVertical2, lat, lon)
  # 
  # # Calculate carbon fluxes at the position of the fish:
  # sim = calcCarbonFluxes(sim) 
  # totalFlux = sim$fluxCarcass + sim$fluxFecal + sim$fluxRepro + sim$fluxRespiration
  # barplot(totalFlux, xlab="Size class", ylab="Flux (gWW/m2/yr)")
  # 
  # # Calculate the injection
  # inject = calcCarbonInjection(sim)
  # 
  # z = -inject$z
  # plot( inject$total, z, type="l", lwd=3, 
  #       xlim=c(0,max(inject$total[1:length(inject$total)-1])),
  #       xlab="Injection (gC/m3/yr)",
  #       ylab="Depth (m)")
  # lines( inject$Fecal, z, col="brown" )
  # lines( inject$Carcass, z, col="grey" )
  # lines( inject$Repro, z, col="darkgreen")
  # lines( inject$Respiration, z, col="darkred")
  # legend("top",
  #        c("Total","Fecal pellets","Carcasses","Reproduction","Respiration"),
  #        lwd=c(3,1,1,1,1),
  #        col=c("black","brown","grey","darkgreen","darkred")
  # )
  # 
  # # Calculate injection on TM grid:
  # TM = loadTransportMatrix()
  # long=lon
  # if (lon<0)
  #   long = 360+lon
  # 
  # injectTM = project_injection_to_TM(inject, lat, long, TM) 
  # plot( injectTM$inject, -TM$grid$zt, ylim=c(2*min(z),0) )
  # 
  # # Assemble a matrix with all injections
  # matrixInject = array(dim=dim(TM$M3d), data=0)
  # matrixInject[injectTM$ix$y, injectTM$ix$x, ] = injectTM$inject
  # 
  # # Solve the transport matrix to get sequestration etc.:
  # sequestration = calc_CarbonSequestration(TM, matrixInject)
  # 
  # #
  # # Plots:
  # #
  # sequestration$lat = grid$xt
  # sequestration$lon = grid$yt
  # 
  # dat = as.data.frame( sequestration$Cseq_per_area )
  # world <- map_data("world2")
  # 
  # image( x=c(grid$xt[1]-1,grid$xt), y=c(-90,grid$yt), z=log10(t(sequestration$Cseq_per_area)))
}
