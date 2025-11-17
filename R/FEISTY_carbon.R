#
# Code for calculating carbon fluxes, carbon injection, and carbon sequestration.
#
# Written by Julie Lemoine and Ken H Andersen.
# Sequestration calculations based on code by Andre W Visser
#
#library(FEISTY)
library(data.table)
library(pracma)
library(Matrix)
library(doParallel)
library(foreach)
library(ggplot2)
library(patchwork) # To arrange two plots
library(maps)

# 
# Correct a longitude from the range -180:180 to 0:360:
#
longitude_correction = function(lon) {
  if (lon < 0)
    lon = 360 + lon
  return(lon)
}

inverse_longitude_correction = function(lon) {
  if (lon > 180)
    lon = lon-360
  return(lon)
}
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
calcCarbonFluxes <- function(sim) {
  
  p = sim$p # parameters
  
  # Fluxes from grazing
  grazrate = p$Cmax[p$ixFish] * sim$f # grazing rate [yr^-1]
  graz = sim$B * grazrate # grazing flux (before assimilation) [gWW.m^-2.yr^-1]
  feces = graz * 0.15 # feces flux [gWW.m^-2.yr^-1] -> (1 - p$epsAssim)
  
  # Reproduction waste goes to fecal pellets
  rep2feces = sim$Repro
  rep2feces[, sapply(p$ix, tail, n = 1) - p$nResources] =
    sim$Repro[, sapply(p$ix, tail, n = 1) - p$nResources] +
    sim$Fout[, sapply(p$ix, tail, n = 1) - p$nResources]
  rep2feces = rep2feces * (1 - unique(p$epsRepro)) * unique(p$epsRepro / 0.22) # eps_egg = 0.22 from Andersen 2019 p47
  
  feces = feces + rep2feces # feces flux [gWW.m^-2.yr^-1]
  
  # Flux from respiration
  resprate = p$metabolism[p$ixFish] + 0.15 * grazrate # respiration rate [yr-1]
  respiration = sim$B * resprate # respiration flux [gWW.m^-2.yr^-1]
  
  # Get last 40% of timeseries
  etaTime <- 0.4 
  ixTime  <- which(sim$t >= ((1 - etaTime) * sim$t[sim$nTime]))
  
  sim$fluxCarcass   <- colMeans(sim$B[ixTime,] * p$mort0[-c(1:p$nResources)]) # deadfalls flux [gWW.m^-2.yr^-1]
  sim$fluxFecal <- colMeans(feces[ixTime,]) # total feces flux [gWW.m^-2.yr^-1]
  sim$fluxRepro  <- colMeans(rep2feces[ixTime,]) # feces flux from reproduction waste [gWW.m^-2.yr^-1]
  sim$fluxRespiration <- colMeans(respiration[ixTime,]) # respiration flux [gWW.m^-2.yr^-1]
  
  return(sim)
}

#
# Calculate injection fluxes calculated at all z-levels
# Units i gC/m3/year
#
calcCarbonInjection = function(sim) {
  rho_gC_gWW = 9 # Gram carbon per gram wet weight
  #
  # Calculate the injection from sinking POC with velocity v:
  #
  calcPOCinject = function(J, v) {
    
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
  
  
  # parameters
  alpha0 = 0.25 # maximum bacterial degradation rate (from Pinti) [day^-1]
  Q10r = 2 # Q10 remin [-]
  
  p = sim$p
  functy = "all"
  
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
  
  # depth indices
  ix100 = ifelse(test = p$bottom >= 100, yes = 100 + 1, no = 0)
  ixbottom = p$bottom + 1
  
  col_indices <- switch(functy,
                        "smallPel" = 5:10,
                        "mesoPel" = 11:16,
                        "largePel" = 17:25,
                        "mwpred" = 26:34,
                        "dem" = 35:43,
                        "all" = 5:43,
                        stop("Unknown functy value : ", functy))
  #
  # Calculation the injection from respiration:
  #
  pDay = p$depthDay[, col_indices]  
  pNight = p$depthNight[, col_indices]
  respiration = t(pDay+pNight) * sim$fluxRespiration/2 # For each size class
  
  # poc production
  resFecal = calcPOCinject( as.vector(sim$fluxFecal), 500)
  resCarcass = calcPOCinject( as.vector(sim$fluxCarcass), 700)
  resRepro = calcPOCinject( as.vector(sim$fluxRepro), 500)
  
  # Make list with injections as output and convert to carbon units:
  res = list()
  res$z = resFecal$depth
  res$Fecal = resFecal$inject / rho_gC_gWW
  res$Carcass = resCarcass$inject / rho_gC_gWW
  res$Repro = resRepro$inject / rho_gC_gWW
  res$Respiration = colSums( respiration ) / rho_gC_gWW
  
  res$total = res$Fecal + res$Carcass + res$Repro + res$Respiration
  
  return(res)
}

#
# Simulate a given position in the global map:

### PERHAPS MOVE TO MAIN FEISTY (including data file)
#
#' @export
simulatePosition = function(setup, lat, lon, nStages=9, tEnd=200) {
  # Output from COBALT
    
    pp = getParametersPosition(lat,lon)
    p = setup(szprod = pp$szprod,
              lzprod = pp$lzprod,
              dfpho  = pp$dfbot,
              depth  = pp$depth,
              Tp     = pp$Tp,
              Tm     = pp$Tm,
              Tb     = pp$Tb,
              nStages = nStages)
  
    sim = simulateFEISTY(p = p, tEnd = 10) 
  return(sim)
}

getParametersPosition = function(lat, lon, sFile="data/Cobalt global data.csv") {
  glob <- read.csv(sFile)
  
  if (lon<0)
    lon = 360+lon
  
  ix = which.min( (glob$lat-lat)^2 + (glob$lon-lon)^2 ) # Find the best fitting location
  
  if ( min((glob$lat-lat)^2 + (glob$lon-lon)^2) < 10)
  {
      szprod = glob[ix, "szprod"]        # small zooplankton production
      lzprod = glob[ix, "lzprod"]        # large zooplankton production
      dfbot  = glob[ix, "dfbot"]         # detrital flux reaching the bottom
      photic = glob[ix, "photic"]        # photic zone depth
      depth  = glob[ix, "depth"]         # water column depth
      Tp     = glob[ix, "Tp"]            # pelagic water temperature
      Tm     = glob[ix, "Tm"]            # mid-water temperature
      Tb     = glob[ix, "Tb"]            # bottom water temperature
  }
  else
  {
    # No need to simulate land points for long:
      szprod = 0        # small zooplankton production
      lzprod = 0        # large zooplankton production
      dfbot  = 0         # detrital flux reaching the bottom
      photic = 200        # photic zone depth
      depth  = 10         # water column depth
      Tp     = 10            # pelagic water temperature
      Tm     = 10            # mid-water temperature
      Tb     = 10            # bottom water temperature
  }
  
  return( list(
    szprod = szprod,        # small zooplankton production
    lzprod = lzprod,        # large zooplankton production
    dfbot  = dfbot,         # detrital flux reaching the bottom
    photic = photic,        # photic zone depth
    depth  = depth,         # water column depth
    Tp     = Tp,            # pelagic water temperature
    Tm     = Tm,            # mid-water temperature
    Tb     = Tb  
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
loadTransportMatrix = function(sFilename="data/CTL.Rdata", bLUdecompose=FALSE) {
  sLUfilename = 'data/LU decomposed TM.Rdata'
  
  # Always do LU decomposition if the LU-decomposed version does not exist:
  if (!file.exists()) bLUdecompose=TRUE
  
  if (bLUdecompose) {
    # Load the original transport matrix:
    load(sFilename)
    
    # Calculation of A = TR - Sink:
    m <- nrow(TM$TR)
    sink <- rep(0,m)
    sink[1:length(TM$msk$hkeep)] <- 1e10 # a strong sink force (1e10) is attributed on surface cells only
    SSINK <- sparseMatrix(i = 1:m, j = 1:m, x = sink) # sink vector on the diagonal of the SSINK matrix
    A <- TM$TR - SSINK # calculation of A matrix
    TM$A = lu(A)
    
    # Save the LU-decomposed version:
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
# Project the injection calculations onto the TM grid by integrating
# over the entire vertical cell
#
# gC/m2/yr for each cell
#
project_injection_to_TM <- function(inject, lat,lon, TM) {
  integral = 0*unique(TM$grid$zt)
  # Find closest grid point:
  ix = list( 
    y = which.min( (lat-TM$grid$yt)^2 ),
    x = which.min( (lon-TM$grid$xt)^2 ))
  # Integrate along the depth:
  for (j in 1:length(TM$grid$zt)) {
    idx = ( (inject$z > TM$grid$zw[j]) 
            & (inject$z <= (TM$grid$zw[j] + TM$grid$dzt[j])))
    integral[j] = trapz( inject$z[idx], inject$total[idx])
    
  }
  return(list(inject=integral, ix=ix))
}

calc_per_area_sum = function(grid, matrix, depthUpper=0) {
  ix = grid$zt>depthUpper
  dz = grid$DZT3d
  dz[!ix] = 0
  return( apply(replace(matrix, is.na(matrix), 0) * TM$grid$DZT3d, c(1,2), sum) )
}

#
# Calculate the amount of carbon sequestered and the sequstration time
#
#' @export
calcCarbonSequestration <- function(TM,  # Transport matrix
                                    matrixInject  # injection matrix (lon, lat, depth) with same dimensions as TM$grid$M3d
){
  
  project_to_TM = function(vector) {
    return( replace(array(NA, dim = dim(TM$M3d)), TM$msk$pkeep, vector) )
  }
  
  # Initialization of result list
  result <- NULL
  
  #
  # Setup grid from TM:
  #
  M3d <- TM$M3d                             # 3D array (lon x lat x depth) containing :
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
  m <- nrow(TM$TR)
  sink <- rep(0,m)
  sink[1:length(msk$hkeep)] <- 1e10 # a strong sink force (1e10) is attributed on surface cells only
  SSINK <- sparseMatrix(i = 1:m, j = 1:m, x = sink) # sink vector on the diagonal of the SSINK matrix
  A <- TM$TR - SSINK # calculation of A matrix
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
  cseq <- solve(A, -q_ocim, sparse=TRUE)
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
  ixBelowEuphotic = TM$grid$zt > 200
  result$inject_below_euphotic = apply( result$inject_per_area[,,ixBelowEuphotic], c(1,2), sum)
  

  # Carbon sequestered on the grid and per area (gC/yr/m3):
  result$Cseq = project_to_TM( cseq )
  result$Cseq_per_area = calc_per_area_sum( TM$grid, result$Cseq ) # (gC/yr/m2)
  
  # Total carbon sequestered in the ocean [PgC]
  TotSeq <- crossprod(V, cseq) / 1e15
  result$TotSeq <- TotSeq
  
  # Sequestration time [year] on the TM grid
  local_export <- q_ocim
  local_export[local_export == 0] <- NA
  SeqTime <- cseq / local_export
  result$SeqTime <- project_to_TM( SeqTime )
  
  # Total sequestration time [year]
  TotSeqTime <- TotSeq / TotInject
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
calcGlobalCarbonSequestration = function(TM=loadTransportMatrix(), lon=c(0,360), lat=c(-90,90), bPrintStatus=TRUE) {

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
  
  # Setup parallel backend:
  cl <- makeCluster(detectCores()-1)
  registerDoParallel(cl)
  
  # Loop over all grid points in the lat/lon range:
  cat("Simulating FEISTY to calculate injections at",dim(grid_idx)[1], "position(s).\n")
  injectTM = foreach(i = 1:dim(grid_idx)[1],
                     .packages = c("FEISTY","pracma"),
                     .verbose = FALSE) %dopar% 
    {
      sim = simulatePosition(setupVertical2, 
                             TM$grid$yt[ grid_idx$i[i] ], 
                             TM$grid$xt[ grid_idx$j[i] ] )
      
      # Calculate carbon fluxes at the position of the fish:
      sim = calcCarbonFluxes(sim) 
      #totalFlux = sim$fluxCarcass + sim$fluxFecal + sim$fluxRepro + sim$fluxRespiration
      
      # Calculate the injection
      inject = calcCarbonInjection(sim)
      
      # Calculate injection on TM grid:
      injectTM = project_injection_to_TM(inject, 
                                         TM$grid$yt[ grid_idx$i[i] ], 
                                         TM$grid$xt[ grid_idx$j[i] ], TM) 
      injectTM$inject
    } 
  stopCluster(cl)
  
  # Put into the injection matrix:
  for (i in 1:dim(grid_idx)[1])
    matrixInject[ grid_idx$i[i], grid_idx$j[i],] = injectTM[[i]]
  
  # Solve the transport matrix to get sequestration etc.:
  cat("Calculating carbon sequestration\n")
  sequestration = calcCarbonSequestration(TM, matrixInject)
  
  # Calculate the per-area sequestration for the cells which are simulated:
  area = 0 # Area of all simulated cells
  for (i in 1:dim(grid_idx)[1])
    area = area + TM$grid$Areat[grid_idx$i[i], grid_idx$j[i]]
  sequestration$TotSeq_per_area = sequestration$TotSeq / area * 1e15 # gC/m2
  
  sequestration$ix_lat = ix_lat
  sequestration$ix_lon = ix_lon
  #
  # Print summary and make plots:
  #
  if (bPrintStatus) {
    cat( c("Total carbon injected: ", format(sequestration$TotInject,digits=3), "pgC/yr \n"))
    cat( c("Total carbon sequestered: ", format(sequestration$TotSeq,digits=3), 'pgC \n') )
    cat( c("Average sequestered per area: ", format(sequestration$TotSeq_per_area,digits=3), 'gC/m2 \n') )
    cat( c("Average sequestration time: ", format(sequestration$TotSeqTime,digits=3), 'yr \n') )
    
    p1 = plotGlobal( sequestration$lon, sequestration$lat, sequestration$Cseq_per_area, 
                     sTitle="Carbon sequestered", "gC/m2")
    
    p2= plotGlobal( sequestration$lon, sequestration$lat, 
                    sequestration$inject_below_euphotic,
                    sTitle="Injection below 200 m", "gC/m2/yr")
    
    p3 = plotGlobal( sequestration$lon, sequestration$lat,
                     sequestration$SeqTime[,,2],
                     sTitle="Sequestration time", "yr")
    
    p1 + p2 + p3
  }
  
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
testCarbonCalculations_one_position = function(lat=60, lon=-15) {
  # Simulate the position 60, -15 using Cobalt output:
  sim = simulatePosition(setupVertical2, lat, lon)
  
  # Calculate carbon fluxes at the position of the fish:
  sim = calcCarbonFluxes(sim) 
  totalFlux = sim$fluxCarcass + sim$fluxFecal + sim$fluxRepro + sim$fluxRespiration
  barplot(totalFlux, xlab="Size class", ylab="Flux (gWW/m2/yr)")
  
  # Calculate the injection
  inject = calcCarbonInjection(sim)
  
  z = -inject$z
  plot( inject$total, z, type="l", lwd=3, 
        xlim=c(0,max(inject$total[1:length(inject$total)-1])),
        xlab="Injection (gC/m3/yr)",
        ylab="Depth (m)")
  lines( inject$Fecal, z, col="brown" )
  lines( inject$Carcass, z, col="grey" )
  lines( inject$Repro, z, col="darkgreen")
  lines( inject$Respiration, z, col="darkred")
  legend("top",
         c("Total","Fecal pellets","Carcasses","Reproduction","Respiration"),
         lwd=c(3,1,1,1,1),
         col=c("black","brown","grey","darkgreen","darkred")
  )
  
  # Calculate injection on TM grid:
  TM = loadTransportMatrix()
  long=lon
  if (lon<0)
    long = 360+lon
  
  injectTM = project_injection_to_TM(inject, lat, long, TM) 
  plot( injectTM$inject, -TM$grid$zt, ylim=c(2*min(z),0) )
  
  # Assemble a matrix with all injections
  matrixInject = array(dim=dim(TM$M3d), data=0)
  matrixInject[injectTM$ix$y, injectTM$ix$x, ] = injectTM$inject
  
  # Solve the transport matrix to get sequestration etc.:
  sequestration = calc_CarbonSequestration(TM, matrixInject)
  
  #
  # Plots:
  #
  sequestration$lat = grid$xt
  sequestration$lon = grid$yt
  
  dat = as.data.frame( sequestration$Cseq_per_area )
  world <- map_data("world2")
  
  image( x=c(grid$xt[1]-1,grid$xt), y=c(-90,grid$yt), z=log10(t(sequestration$Cseq_per_area)))
}
