
#===============================================================================
# Functions to setup model applications
# (generate specific MODEL parameter lists)
#===============================================================================

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
                      bprod  = 5,   # benthos production?
                      depth  = 100, # water column depth [m]
                      Tp     = 10,  # pelagic layer averaged temperature [Celsius]
                      Tb     = 8)  # bottom layer depth [Celsius]
  {
  # Initialize the parameters:
  param = paramInit(depth=depth, szprod=szprod, lzprod=lzprod, bprod=bprod,Tp=Tp,Tb=Tb)
  
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
                        mortF=c(0,0.03,0.3), nStages=3, name="Demersals")
  
  # physiology of all fish stages
  param = paramAddPhysiology(param)
  
  param=paramTeffect(param, # only for setupbasic & 2
                           Tref=10,
                           Q10=1.88,
                           Q10m=2.35,
                           u=NA)
  # Tref=10
  # Q10=1.88
  # Q10m=2.35 #Petrik et al.,2019
  # 
  # fTemp=Q10^((Tp - Tref)/10)
  # fTempm=Q10m^((Tp - Tref)/10)
  # fTempdem=Q10^((Tb - Tref)/10)
  # fTempmdem=Q10m^((Tb - Tref)/10)
  # 
  # ix = c(param$ix[[1]],param$ix[[2]])
  # m = param$mc[ix]
  # param$Cmax[ix] = fTemp* param$Cmax[ix] # maximum consumption rate 
  # param$V[ix] = fTemp* param$V[ix] # clearance rate 
  # param$metabolism[ix] = fTempm* param$metabolism[ix] # 0.2*param$Cmax[ix] # standard metabolism
  # ix = param$ix[[3]]
  # m = param$mc[ix]
  # param$Cmax[ix] = fTempdem* param$Cmax[ix] # maximum consumption rate 
  # param$V[ix] = fTempdem* param$V[ix] # clearance rate 
  # param$metabolism[ix] = fTempmdem* param$metabolism[ix] # 0.2*param$Cmax[ix] # standard metabolism 
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
  param$theta["smallPel_2", "Demersals_1"] = 1

  # Large pelagics:
  param$theta["largePel_1", "smallZoo"] = 1    
  param$theta["largePel_2", "smallZoo"] = 0.25 
  param$theta["largePel_2", "largeZoo"] = 1   
  param$theta["largePel_2", c("smallPel_1", "largePel_1", "Demersals_1")] = 1 
  param$theta["largePel_3", "smallPel_2"] = 0.5 
  param$theta["largePel_3", "largePel_2"] = 1 

  # Demersals:
  param$theta["Demersals_1", "smallZoo"] = 1
  param$theta["Demersals_2", "smallBenthos"] = 1
  if (param$depth < 200){
  param$theta["Demersals_3", "smallPel_2"] = 1
  param$theta["Demersals_3", "smallPel_2"] = 0.75/2
  param$theta["Demersals_3", "largePel_2"] = 0.75 
  }
  param$theta["Demersals_3", "Demersals_2"] = 1 

  return(param)
}

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
                       bprod  = 5,   # benthos production?
                       depth  = 100, # water column depth [m]
                       Tp     = 10,  # pelagic layer averaged temperature [Celsius]
                       Tb     = 8,  # bottom layer depth [Celsius]
                       nStages=9,
                       etaMature=0.25) {
  
  # Initialize the parameters:
  param = paramInit(depth=depth, szprod=szprod, lzprod=lzprod, bprod=bprod,Tp=Tp,Tb=Tb,etaMature=etaMature)

  # Setup resource groups:
  param = paramAddResource(
        param, 
        names= c("smallZoo", "largeZoo", "smallBenthos", "largeBenthos"),
        K    = c(szprod, lzprod, bprod, 0),  # g ww/m2  - maximum resource concentration
        r    = c(1, 1, 1, 1),              # [/yr] nudging coefficient
        mc   = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500)))


# Add fish groups:
  nSmall = round(0.66*nStages)
  param = paramAddGroup(param, mMin=0.001, mMax=   250, mMature=etaMature*250, 
                        mortF=0,      nStages=nSmall, name="smallPel")
  
  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=etaMature*125000, 
                        mortF=0, nStages=nStages, name="largePel") 
  
  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=etaMature*125000, 
                        mortF=0, nStages=nStages, name="Demersals")

  # physiology of all fish stages
  param = paramAddPhysiology(param)
  
  param=paramTeffect(param, # only for setupbasic & 2
                     Tref=10,
                     Q10=1.88,
                     Q10m=2.35,
                     u=NA)
  
  # Tref=10
  # Q10=1.88
  # Q10m=2.35 #Petrik et al.,2019
  # 
  # fTemp=Q10^((Tp - Tref)/10)
  # fTempm=Q10m^((Tp - Tref)/10)
  # fTempdem=Q10^((Tb - Tref)/10)
  # fTempmdem=Q10m^((Tb - Tref)/10)
  # 
  # ix = c(param$ix[[1]],param$ix[[2]])
  # m = param$mc[ix]
  # param$Cmax[ix] = fTemp* param$Cmax[ix] # maximum consumption rate
  # param$V[ix] = fTemp* param$V[ix] # clearance rate
  # param$metabolism[ix] = fTempm* param$metabolism[ix] # 0.2*param$Cmax[ix] # standard metabolism
  # ix = param$ix[[3]]
  # m = param$mc[ix]
  # param$Cmax[ix] = fTempdem* param$Cmax[ix] # maximum consumption rate
  # param$V[ix] = fTempdem* param$V[ix] # clearance rate
  # param$metabolism[ix] = fTempmdem* param$metabolism[ix] # 0.2*param$Cmax[ix] # standard metabolism

  # Setup size interaction matrix:
  thetaA = 0.5  # Large fish pref for medium forage fish
  thetaD = 0.75 # Pref of large demersal on pelagic prey
  
  param$theta = matrix(nrow=param$nStages, ncol=param$nStages, data=0)
  rownames(param$theta) <- colnames(param$theta) <- param$stagenames

   # Size-based interactions:  
   beta  = 400
   sigma = 1.3

   for (i in param$ixFish) {
     param$theta[i,] = exp( -(log(param$mc[i]/(beta*param$mc)))^2 / (2*sigma)^2  )
     param$theta[i, param$mc > param$mc[i]] = 0
   }
   param$theta[is.na(param$theta)] = 0
  
   #
   # Setup interactions between groups and resources:
   #
   ixSmall = param$ix[[1]]
   ixLarge = param$ix[[2]]
   ixDem   = param$ix[[3]]
  
   mMedium = 0.5
   mLarge = 250
   ixSmallSizeDem = ixDem[ (param$mc[ixDem]<=mMedium) ]
   ixMediumSizeDem = ixDem[ (param$mc[ixDem]>mMedium) &
                            (param$mc[ixDem]<mLarge) ]
   
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
   
   # Large demersals feed have reduced feeding effiiency on pelagic species:
   param$theta[ixDem, ixSmall] = thetaA * thetaD * param$theta[ixDem,ixSmall] 
   param$theta[ixDem, ixLarge] = thetaD * param$theta[ixDem, ixLarge] 

  #
  # Mortality
  #
  param$mort0 = 0.1
  param$mortF[param$ixFish] = 0.#*c(0,1,0,0,1,0,0,1) # Fishing only on mature stages
  
  
  param$metabolism[is.na(param$metabolism)]=0
  param$mortF[is.na(param$mortF)]=0
  param$psiMature[is.na(param$psiMature)]=0
  param$z[is.na(param$z)]=0
  
  return(param)
}

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
                         bent = 150, # Detrital flux out of photic zone
                         nStages=6, # No. of size groups
                         region = 4, # Temperature profile regions: 1 Tropical, 2 Temperate, 3 Boreal, 4 Default 10 Celsius 
                         depth=1500, # Bottom depth
                         photic=150, # Photic zone depth
                         mesopelagic=250, # mesopelagic depth
                         visual=1.5,# >1 visual predation primarily during the day, = 1 equal day and night
                         etaMature = 0.25 # Size of matureation relative to
                                          # asymptotic size. Different from
                                          # van Denderen (2021), where it is 0.002
                               ) {
  
  #------------------  
  # Initialize the parameters:
  # habitat and small benthos
  #------------------  
  bprod=0.1*(bent*(depth/photic)^-0.86)
  if(bprod>=0.1*bent)bprod=0.1*bent

  param = paramInit(bottom=depth, szprod=szprod, lzprod=lzprod, photic=photic,
                    mesop=mesopelagic, visual=visual, bent=bent,etaMature=etaMature,region=region)
  
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
        u0     = 0.5)

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
                        mortF=0, nStages=nStages, name="Demersals")
  param$mortF[length(param$mortF)]=0.5
  #------------------  
  # Setup physiology:
  #------------------  
  param = paramAddPhysiology(param,am = 0.2*20) # 20% * Max. consumption coefficient)
    
  #------------------  
  #------------------  
  # theta (preferences):
  #------------------  
  #------------------  
  param$theta  = matrix(nrow=param$nStages, ncol=param$nStages, data=0)
  rownames(param$theta) <- colnames(param$theta) <- param$stagenames

  beta  = 400
  sigma = 1.3
  param$sizeprefer = matrix(nrow=param$nStages, ncol=param$nStages, data=0)
  param$vertover   = matrix(nrow=param$nStages, ncol=param$nStages, data=0)
  
  # calculate size-preference matrix
  for (i in param$ixFish[[1]]: param$nStages){
       for (j in 1: param$nStages){
            param$sizeprefer[i, j] = sqrt(pi/2)*sigma*(
              erf((log(param$mUpper[j]) - log(param$mc[i]/beta))/(sqrt(2)*sigma))
                  - erf((log(param$mLower[j]) - log(param$mc[i]/beta))/(sqrt(2)*sigma)))
  param$sizeprefer[i, j] = param$sizeprefer[i, j]/(log(param$mUpper[j]) - log(param$mLower[j]))
  }
}

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
  
  ixjuv = which.min(abs(param$mLower[param$ix[[5]]] - 0.5))
  ixadult = which.min(abs(param$mLower[param$ix[[5]]] - param$mMature[[5]]))

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

  xlocvec[ixadult:length(xlocvec)] = param$dvm ### or demmig???
  dem_d =  VertDist(sigmap[ix], xlocvec)

  #if shallower than euphotic depth, adult demersals feed across-habitats
  if (param$bottom <= param$photic) {
  dem_d = (dem_d + dem_n)/2
  dem_n = dem_d
  }
  
  # calculate overlap during day
  depthDay = matrix(nrow=length(xrange), ncol=param$nStages, data=0)
  test     = matrix(nrow=length(xrange), ncol=param$nStages, data=0)
  param$dayout = matrix(nrow=param$nStages, ncol=param$nStages, data=0)
  
  depthDay[, 1:2] = zp_d # resources
  depthDay[, 3:4] = bent_dn # resources
  depthDay[, param$ix[[1]]] = spel_dn
  depthDay[, param$ix[[2]]] = mpel_d
  depthDay[, param$ix[[3]]] = lpel_d
  depthDay[, param$ix[[4]]] = bpel_d
  depthDay[, param$ix[[5]]] = dem_d
  
  for (i in 1: param$nStages) {
    for ( j in 1: param$nStages) {
     test[, j] = pmin(depthDay[, i], depthDay[, j])
    }
  param$dayout[, i] = colSums(test)
  }
  
  # calculate overlap during night
  depthNight = matrix(nrow=length(xrange), ncol=param$nStages, data=0)
  # test will be overwritten
  param$nightout = matrix(nrow=param$nStages, ncol=param$nStages, data=0)
  
  depthNight[, 1:2] = zp_n # resources
  depthNight[, 3:4] = bent_dn # resources
  depthNight[, param$ix[[1]]] = spel_dn
  depthNight[, param$ix[[2]]] = mpel_n
  depthNight[, param$ix[[3]]] = lpel_n
  depthNight[, param$ix[[4]]] = bpel_n
  depthNight[, param$ix[[5]]] = dem_n
  
  for (i in 1: param$nStages) {
    for ( j in 1: param$nStages) {
      test[, j] = pmin(depthNight[, i], depthNight[, j])
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
  idx_be = param$ixFish[1]: (param$ix[[5]][1] + (ixjuv - 2)) # all pelagic and larval demersals
  param$theta[idx_be, 3:4] = 0 # all pelagic and larval demersals do not eat benthos,
                               # only juvenile & adult demersals eat benthos
  
  # small demersals are less preyed on
  idx_smd = (param$ix[[5]][1] + (ixjuv - 1)): (param$ix[[5]][1] + (ixadult - 2)) #
  param$theta[idx_be, idx_smd] = param$theta[idx_be, idx_smd]*0.25
  
  # juvenile & adult demersals do not eat zooplankton
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
  tempdata=read.table(system.file("data", "tempdata.dat", package = "FeistyR"), sep=',') #
  tempdata[,5]=10 #
  Q10=1.88
  Q10m=1.88
  
  dist=(depthDay+depthNight)/2
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


