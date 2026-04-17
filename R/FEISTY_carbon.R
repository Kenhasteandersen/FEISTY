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
  
  # Fluxes from grazing
  grazrate = p$Cmax[p$ixFish] * sim$f # grazing rate [yr^-1]
  graz = sim$B * grazrate # grazing flux (before assimilation) [gWW.m^-2.yr^-1]
  feces = graz * (1 - p$epsAssim)/2 # feces flux [gWW.m^-2.yr^-1] -> (1 - p$epsAssim)/2 A half of unassimilated food is feces.
  
  # Waste energy from the total energy invested into reproduction to eggs goes to respiration. 
  # Dead eggs from eggs to larvae becomes detritus, so goes to fecal pellets.
  
  # sim$Repro already includes Fout of last stage for each functional type
  eps_egg   = 0.22 # 1-eps_egg is The fraction of reproductive invested used for respiration
  eps_R     = unique(p$epsRepro) / eps_egg # The fraction of eggs that survives
  rep2resp  = sim$Repro * (1 - eps_egg) # metabolic cost of egg production [gWW.m^-2.yr^-1]
  rep2feces = sim$Repro * eps_egg * (1 - eps_R) # dead eggs sink as feces; eps_egg = 0.22 from Andersen 2019 p47
  
  # Flux from respiration
  resprate = p$metabolism[p$ixFish] + (1 - p$epsAssim)/2 * grazrate # respiration rate [yr-1] A half of unassimilated food is specific dynamic action.
  respiration = sim$B * resprate + rep2resp # metabolic cost of egg production [gWW.m^-2.yr^-1]
  
  # Get last 40% of timeseries
  etaTime <- 0.4 
  ixTime  <- which(sim$t >= ((1 - etaTime) * sim$t[sim$nTime]))
  
  sim$fluxCarcass     <- colMeans(sim$B[ixTime,] * p$mort0[-c(1:p$nResources)]) # deadfalls flux [gWW.m^-2.yr^-1]
  sim$fluxFecal       <- colMeans(feces[ixTime,]) # total feces flux [gWW.m^-2.yr^-1]
  sim$fluxRepro       <- colMeans(rep2feces[ixTime,]) # feces flux from reproduction waste [gWW.m^-2.yr^-1]
  sim$fluxRespiration <- colMeans(respiration[ixTime,]) # respiration flux [gWW.m^-2.yr^-1]
  
  return(sim)
}

#
# Calculate injection fluxes calculated at all z-levels
# Units i gC/m3/year
#
#' @export
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
    zeta_X <- rep(0, length(z)) # soa\q 1urce = poc production at each depth
    
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
simulatePosition = function(setup,
                            lat, lon,
                            Fmax=0, ixGroups=NULL, # Specification of fishing parameters sent to setFishing()
                            nStages=9, tEnd=200, glob=NULL) {
  # Output from COBALT

    pp = getParametersPosition(lat, lon, glob=glob)
    p = setup(szprod = pp$szprod,
              lzprod = pp$lzprod,
              dfbot  = pp$dfbot,
              depth  = pp$depth,
              Tp     = pp$Tp,
              Tm     = pp$Tm,
              Tb     = pp$Tb,
              nStages = nStages)
    #
    # Set fishing to: Fmax with a Q10=1.8 at depths < 1000 m,
    # and to 90% less at deeper than 1000 m.
    #
    if (length(ixGroups)>0) {
      Fmax = Fmax * 1.8^((pp$Tp-15)/10) # Q10=1.8 correction
      if (pp$depth>1000)
        Fmax = 0.1*Fmax
        
      p = setFishing(p,Fmax, groupidx=ixGroups)
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
# Project the injection calculations onto the TM grid by integrating
# over the entire vertical cell
#
# gC/m2/yr for each cell
#
project_injection_to_TM <- function(inject, lat,lon, grid) {
  integral = 0*unique(grid$zt)
  # Find closest grid point:
  ix = list( 
    y = which.min( (lat-grid$yt)^2 ),
    x = which.min( (lon-grid$xt)^2 ))
  # Integrate along the depth:
  for (j in 1:length(grid$zt)) {
    idx = ( (inject$z > grid$zw[j]) 
            & (inject$z <= (grid$zw[j] + grid$dzt[j])))
    integral[j] = trapz( inject$z[idx], inject$total[idx])
    
  }
  return(list(inject=integral, ix=ix))
}

calc_per_area_sum = function(grid, matrix, depthUpper=0) {
  ix = grid$zt>depthUpper
  dz = grid$DZT3d
  dz[!ix] = 0
  return( apply(replace(matrix, is.na(matrix), 0) * grid$DZT3d, c(1,2), sum) )
}

#
# Calculate the amount of carbon sequestered and the sequstration time
#
#' @export
calcCarbonSequestration <- function(TM,  # Transport matrix
                                    matrixInject,  # injection matrix (lon, lat, depth) with same dimensions as TM$grid$M3d
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
calcGlobalCarbonSequestration = function(TM=loadTransportMatrix(),
                                         lon=c(0,360), lat=c(-90,90),  # Which latitudes to simulate over
                                         Fmax=0, ixGroups=NULL,        # Specification of fishing (set to setFishing())
                                         bPrintStatus=TRUE,
                                         nCores=detectCores()-2) {

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
  
  # Setup parallel backend:
  cl <- makeCluster( nCores )
  registerDoParallel(cl)

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
                             glob=glob)

      # Calculate carbon fluxes at the position of the fish:
      sim = calcCarbonFluxes(sim)
      #totalFlux = sim$fluxCarcass + sim$fluxFecal + sim$fluxRepro + sim$fluxRespiration

      # Calculate the injection
      inject = calcCarbonInjection(sim)

      # Calculate injection on TM grid:
      injectTM = project_injection_to_TM(inject,
                                         grid$yt[ grid_idx$i[i] ],
                                         grid$xt[ grid_idx$j[i] ], grid)

      ix = sim$t>0.5*max(sim$t) # Last half of the timeseries


      #injectTM$inject
      list( inject=injectTM$inject, SSB=colMeans(sim$SSB[ix,]), Y=colMeans(sim$yield[ix,]),
            photic=sim$p$photic )
    }
  stopCluster(cl)
  tSim = (proc.time() - tStart)[3]
  cat("FEISTY simulations completed in", round(tSim, 1), "seconds.\n")

  # Put into the injection matrix:
  SSB = array(data=0, c(dim(matrixInject)[1:2], 5))
  Yield = array(data=0, c(dim(matrixInject)[1:2], 5))
  photic_map = matrix(200, nrow = length(TM$grid$yt), ncol = length(TM$grid$xt))
  for (i in 1:dim(grid_idx)[1]) {
    matrixInject[ grid_idx$i[i], grid_idx$j[i],] = injectTM[[i]]$inject
    SSB[ grid_idx$i[i], grid_idx$j[i],] = injectTM[[i]]$SSB
    Yield[grid_idx$i[i], grid_idx$j[i],] = injectTM[[i]]$Y
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
  sequestration$SSB = SSB
  sequestration$Yield = Yield

  # Calculate the per-area sequestration for the cells which are simulated:
  area = 0 # Area of all simulated cells
  for (i in 1:dim(grid_idx)[1])
    area = area + TM$grid$Areat[grid_idx$i[i], grid_idx$j[i]]
  sequestration$TotSeq_per_area = sequestration$TotSeq / area * 1e15 # gC/m2
  #
  # Calc total biomass and yield:
  #
  for (i in 1:dim(sequestration$SSB)[3]) {
    sequestration$TotSSB[i] = sum( TM$grid$Areat*sequestration$SSB[,,i], na.rm=TRUE ) / 1e15 # Pg_WW
    sequestration$TotYield[i] = sum( TM$grid$Areat*sequestration$Yield[,,i], na.rm=TRUE ) / 1e15 # Pg_WW/yr
  }
  
  sequestration$ix_lat = ix_lat
  sequestration$ix_lon = ix_lon
  #
  # Print summary and make plots:
  #
  if (bPrintStatus) {
    cat( c("Total fish biomass: ", format(sum(sequestration$TotSSB),digits=3), "PgWW \n") )
    cat( c("Total fish yield: ", format(sum(sequestration$TotYield),digits=3), "PgWW/yr \n") )
    cat( c("Total carbon injected: ", format(sequestration$TotInject,digits=3), "pgC/yr \n"))
    cat( c("Total carbon injected below euphotic: ", format(sequestration$TotInject_below_euphotic,digits=3), "pgC/yr \n"))
    cat( c("Total carbon sequestered: ", format(sequestration$TotSeq,digits=3), 'pgC \n') )
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
