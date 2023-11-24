#===============================================================================
# Global routines for the Shiny web application. 
# You can run the application by clicking  the 'Run App' button above, OR:.
# use webFeisty() to trigger it from anywhere (if the FeistyR package is loaded)
#===============================================================================

#-------------------------------------------------------------------------------
# Background tools for plotting
#-------------------------------------------------------------------------------

singlewidth <- 8/2.54 # panel width in cm
doublewidth <- 13/2.54
height <- 8/1.6/2.54 + 2.5*0.138889

if(.Platform$OS.type=="windows") {
  quartz<-function(width, height) 
    windows(width=width, height=height)
}

# Base R plotting tools =====================================================

TeX <- function(...) latex2exp::TeX(...)
bottom <- 1
left <- 2
top <- 3
right <- 4

ticklength <- 0.2 # tick mark length
omargin <- 0.7 # outer margins
cex <- 1#0.9
cexaxis <- 1#0.8 # Scaling of size of numbers on axes
axis.lwd <- 0.8

iPlot <- 1 # Static variable used for labels on plots

solid = 1
thick = 2
dashed <- 2
dotted <- 3
dashdotted = 4
dots <- 16
triangles <- 17
squares <- 15
stdgrey <- grey(0.4)
darkgrey <- grey(0.15)
lightgrey <- grey(0.7)
black <- grey(0)


#===============================================================================
# Global routines for the Shiny web application. 
# You can run the application by clicking  the 'Run App' button above, OR:.
# use webFeisty() to trigger it from anywhere (if the FeistyR package is loaded)
#===============================================================================

calcDerivativesF = function(t, y, p, FullOutput, ...) 
  simulateFeisty(p=p, times=NA, yini=y, USEdll=TRUE, ...)
  
calcDerivativesR = function(t, y, p, FullOutput, ...) 
  simulateFeisty(p=p, times=NA, yini=y, USEdll=FALSE, ...)

#-------------------------------------------------------------------------------
# A plot of the Biomass of all functional groups
# as a function time.
#-------------------------------------------------------------------------------
plotBiomasstime = function(sim, bPlot=TRUE) {
  if (bPlot)
    defaultplot()
  p = sim$p
  
  Btime = matrix(nrow=sim$nTime, ncol=p$nGroups)
  for(i in 1:p$nGroups)
    Btime[,i]=rowSums( sim$B[,p$ix[[i]]-p$nResources] )
  
  semilogypanel(xlim=sim$t, ylim=c(1E-3,max(max(sim$SSB)*10,max(sim$R)*10)),
                xlab="Time (yr)", ylab = "Biomass (gww  m-2)")
  
  colnames(Btime)=p$groupnames[-p$ixR]
  #
  # Plot fish
  #
  for (i in 1:p$nGroups)
    lines(sim$t, Btime[,i], col= sim$p$my_palette[attr(sim$p$my_palette,"name") %in% sim$p$groupnames[-sim$p$ixR]] [i],lwd = 3)
  #
  # Plot resources:
  #
  for (i in p$ixR)
    lines(sim$t, sim$R[,i], col=sim$p$my_palette[attr(sim$p$my_palette,"name") %in% sim$p$groupnames[sim$p$ixR]] [i],lwd = 3)
  
  # legend(x='bottomright',
  #        legend=c('Resources','Fish'),
  #        lty=c(1,1),
  #        col=c('blue','black'),
  #        bty='n')
  
  legend(x='bottom',
         legend=sim$p$my_names[attr(sim$p$my_names,"names") %in% sim$p$groupnames],
         lty=c(1,1),
         col=sim$p$my_palette[attr(sim$p$my_palette,"names") %in% sim$p$groupnames],
         bty='n',
         ncol = 6, cex = 1,
         lwd = 3)
  
}

#-------------------------------------------------------------------------------
# Makes a basic plot of the adult biomass (SSB) of all functional groups
# as a function time.
#-------------------------------------------------------------------------------

plotSSBtime = function(sim, bPlot=TRUE) {
  if (bPlot)
    defaultplot()
  p = sim$p
  
  semilogypanel(xlim=sim$t, ylim=c(1E-3,max(max(sim$SSB)*10,max(sim$R)*10)),
                xlab="Time (yr)", ylab = "SSB (gww  m-2)")
  #
  # Plot fish
  #
  for (i in 1:p$nGroups)
    lines(sim$t, sim$SSB[,i],col= sim$p$my_palette[attr(sim$p$my_palette,"name") %in% sim$p$groupnames[-sim$p$ixR]] [i],lwd = 3)
  #
  # Plot resources:
  #
  for (i in p$ixR)
    lines(sim$t, sim$R[,i], col=sim$p$my_palette[attr(sim$p$my_palette,"name") %in% sim$p$groupnames[sim$p$ixR]] [i],lwd = 3)

  # legend(x='bottomright',
  #        legend=c('Resources','Fish'),
  #        lty=c(1,1),
  #        col=c('blue','black'),
  #        bty='n')
  
  # legend(x='bottom',
  #        legend=sim$p$my_names[attr(sim$p$my_names,"names") %in% sim$p$groupnames],
  #        lty=c(1,1),
  #        col=sim$p$my_palette[attr(sim$p$my_palette,"names") %in% sim$p$groupnames],
  #        bty='n',
  #        ncol = 4, cex = 1,
  #        lwd = 3)
  
}


#-------------------------------------------------------------------------------
# Makes a basic plot of the Yield of all functional groups
# as a function time.
#-------------------------------------------------------------------------------

plotYieldtime = function(sim, bPlot=TRUE) {
  if (bPlot)
    defaultplot()
  p = sim$p
  colnames(sim$yield)=p$groupnames[-p$ixR]
  semilogypanel(xlim=sim$t, ylim=c(1E-3,max(max(sim$yiled)*10,max(sim$R)*10)),
                xlab="Time (yr)", ylab = "Yield (gww  m-2)")
  #
  # Plot fish
  #
  for (i in 1:p$nGroups)
    lines(sim$t, sim$yield[,i], col= sim$p$my_palette[attr(sim$p$my_palette,"name") %in% sim$p$groupnames[-sim$p$ixR]] [i],lwd = 3)
  
  # legend(x='bottomright',
  #        legend=c('Resources','Fish'),
  #        lty=c(1,1),
  #        col=c('blue','black'),
  #        bty='n')
  
  legend(x='bottom',
         legend=sim$p$my_names[attr(sim$p$my_names,"names") %in% sim$p$groupnames],
         lty=c(1,1),
         col=sim$p$my_palette[attr(sim$p$my_palette,"names") %in% sim$p$groupnames],
         bty='n',
         ncol = 6, cex = 1,
         lwd = 3)
  
}

#-------------------------------------------------------------------------------
# Plots the mortalities and feeding levels
#-------------------------------------------------------------------------------

plotRates = function(p, u=p$u0, bPlot=TRUE) {
#  if (p$USEdll) {
#    rates = calcDerivativesF(0, u, p, FullOutput = TRUE)[[2]]
#  }else{
    rates = calcDerivativesR(0, u, p, FullOutput = TRUE)
#  }
    
  if (bPlot)
    defaultplot(mfcol=c(3,1))
  xlim = range(p$mc[p$ixFish])
  
  #
  # Growth rate
  # 
  loglogpanel(xlim = xlim, ylim=c(1e-1,max(rates$g)*10),
                ylab="Growth rate (1/year)", xlab="-", xaxis = FALSE)
  for (i in 1:p$nGroups) {
    lines(p$mc[p$ix[[i]]], rates$g[p$ix[[i]]-length(p$ixR)], lwd=i,
          col=p$my_palette[attr(p$my_palette,"name") %in% p$groupnames[-p$ixR]] [i])
  }

  #
  # Mortalities:
  #
  loglogpanel(xlim=xlim, ylim=c(1E-3,max(rates$mortpred)+10),
              xlab="-", ylab="mort (1/year)", xaxis = FALSE)
  for (i in 1:p$nGroups) {
    lines(p$mc[p$ix[[i]]], rates$mortpred[p$ix[[i]]], lwd=i, col='red')
    lines(p$mc[p$ix[[i]]], p$mortF[p$ix[[i]]], lwd=i, col='blue')
  }
  hline(p$mort0[p$ix[[i]]][1])
  
  legend(x='bottomleft',
         legend=c('Predation','Fishing','Background'),
         lty=c(1,1,dotted),
         col=c('red','blue','black'),
         bty='n')
  #
  # Feeding level
  # 
  semilogxpanel(xlim = xlim, ylim=c(0,1),
                ylab="Feeding level, f", xlab="Mass (gww)")
  for (i in 1:p$nGroups) {
    ix = p$ix[[i]]
    lines(p$mc[ix], rates$f[ix], lwd=i,
          col=p$my_palette[attr(p$my_palette,"name") %in% p$groupnames[-p$ixR]] [i])
    #lines(p$mc[ix], p$metabolism[ix]/(p$epsAssim*p$Cmax[ix]), lwd=i, lty=dotted) # Critical feeding level
  }
  
  i = which.max(sapply(p$ix,"length"))
  ix = p$ix[[i]]
  lines(p$mc[ix], p$metabolism[ix]/(p$epsAssim*p$Cmax[ix]), lwd=1, lty=dotted) # Critical feeding level
  
  legend(x='bottomright',
         legend=c('Critical feeding level'),
         lty=c(dotted),
         col=c('black'),
         bty='n')
  
  return(rates)
}

#-------------------------------------------------------------------------------
# Plot the biomasses of all groups
#-------------------------------------------------------------------------------

plotSpectra = function(sim, iTime=sim$nTime, bPlot=TRUE) {
  if (bPlot)
    defaultplot()
  p = sim$p
  
  loglogpanel(xlim=p$mc[p$ixFish], ylim=c(1e-3,max(colMeans(sim$B[round(0.6*iTime):iTime,])*10)),
               ylab="Biomass (gww m-2)", xlab="-", xaxis = FALSE)
  
  # Total biomass spectra
  totBspec= matrix(nrow=1, ncol=max(sapply(p$ix, length)),data=0)
  for (i in 1:p$nGroups) {
    totBspec[1:length(p$ix[[i]])] =totBspec[1:length(p$ix[[i]])]+ colMeans(sim$B[round(0.6*iTime):iTime,p$ix[[i]]-p$ixFish[1]+1])
  }
  
  lines(p$mc[p$ix[[which.max(sapply(p$ix, length))]]], totBspec,
        col=black,
        lwd=3)
  
  # mean biomass of last 40% of simulation time  
  for (i in 1:p$nGroups) {
    lines(p$mc[p$ix[[i]]], colMeans(sim$B[round(0.6*iTime):iTime,p$ix[[i]]-p$ixFish[1]+1]),
          col=sim$p$my_palette[attr(sim$p$my_palette,"name") %in% sim$p$groupnames[-sim$p$ixR]] [i],
          lwd=i)
  }
  
}

#-------------------------------------------------------------------------------
# Add legends in an independent plot panel
#-------------------------------------------------------------------------------
addLegends=function(sim){
  defaultpanel(xlim=c(0,10),ylim=c(0,10),xlab='', ylab='',
             xaxis=FALSE, yaxis=FALSE, label=FALSE,bty="n")
  
  legend(x='top',
         legend=sim$p$my_names[attr(sim$p$my_names,"names") %in% sim$p$groupnames],
         lty=1,
         col=sim$p$my_palette[attr(sim$p$my_palette,"names") %in% sim$p$groupnames],
         bty='n',
         ncol = 4, cex = 1.2,
         lwd = 3) 
}


#-------------------------------------------------------------------------------
# Make 4 panels of simulation
#-------------------------------------------------------------------------------

plotSimulation = function(sim) {
  defaultplot(mfcol=c(6,1))
  plotSSBtime(sim,bPlot=FALSE)
  plotSpectra(sim, bPlot=FALSE)
  rates = plotRates(sim$p, u=c( sim$R[sim$nTime,], sim$B[sim$nTime,]),bPlot=FALSE)
  addLegends(sim)
  return(rates)
}

#-------------------------------------------------------------------------------
# Plot the interaction matrix:
#-------------------------------------------------------------------------------

plotTheta = function(p) {
  par(mfcol=c(1,1))

  image( t(p$theta[p$ixFish,]), 
         x = seq(1,max(p$ixFish)), 
         y = seq(1, length(p$ixFish)),
         xlab="Prey group", ylab="Predator group" )
  #
  # Add a line for the resource groups
  #
  y = c(0,max(p$ixFish))
  lines(x = 0.5+c(1,1)*max(p$ixR), y = y, col="red")
  #
  # Add a line for each fish group:
  #
  for (iGroup in 1:p$nGroups) {
    lines(x = 0.5+c(1,1)*max(p$ix[[iGroup]]), y-max(p$ixR))
    lines(x=y, y=c(1,1)*max(p$ix[[iGroup]])-max(p$ixR)+0.5)
  }
}

#-------------------------------------------------------------------------------
# Network plot 
# Revised Based on work of Daniel Ottmann Riera and Solenne Roux
# u: sim$u simpleOutput=TRUE
#-------------------------------------------------------------------------------

plotNetwork <- function(p, u) {
  
  # Number of groups and biomass:
  ngroup <- p$ix[[length(p$ix)]][length(p$ix[[length(p$ix)]])] #ressources (4) + fish 
  biomass <-u# sim[,2:(p$ix[[length(p$ix)]][length(p$ix[[length(p$ix)]])]+1)] #ressources + poissons car se mangent entre eux (+1 car première colonne = temps)
  
  #Average of the biomass : 
  Bi <- colMeans(biomass[round(0.6*nrow(biomass), digits = 0):nrow(biomass),]) # mean value of the last 20% time 
  
  if (p$setup == "setupBasic"){
    Av_depth <- c(-1,-1,-4,-4,0,0,-2,-2,-2,-3,-3,-3)
    
    p$SpId <- c('smallPel','largePel', 'demersals')
    SpId <- c("smallZoo", "largeZoo", "smallBenthos", "largeBenthos", 
              rep(p$SpId[1], length(p$ix[[1]])),
              rep(p$SpId[2], length(p$ix[[2]])),
              rep(p$SpId[3], length(p$ix[[3]])))

  }  
  
  if (p$setup == "setupBasic2"){
    Av_depth <- c(-1,-1,-4,-4,rep(0, length(p$ix[[1]])),rep(-2, length(p$ix[[2]])),rep(-3, length(p$ix[[3]])))
    
    p$SpId <- c('smallPel','largePel', 'demersals')
    SpId <- c("smallZoo", "largeZoo", "smallBenthos", "largeBenthos", 
              rep(p$SpId[1], length(p$ix[[1]])),
              rep(p$SpId[2], length(p$ix[[2]])),
              rep(p$SpId[3], length(p$ix[[3]])))
    
  }  
  
  
  
  if (p$setup == "setupVertical" | p$setup == "setupVertical2"){
    
    #Calculate average depth day/night
    
    Av_depth_day <- 1 : p$nStages
    Av_depth_night <- 1 : p$nStages
    for (i in 1:p$nStages) {
      Av_depth_day[i] <- which.max(p$depthDay[ ,i])
      Av_depth_night[i] <- which.max(p$depthNight[ ,i]) 
      
    }
    #ce calcul prend l'indice de la matrice p$depthDay qui représente en réalité la profondeur et la valeur correspondant à l'indice est la probabilité de trouver x poisson à cette profondeur
    Av_depth <- -(Av_depth_day + Av_depth_night) / 2
    
    
    # Change a bit for visualization:
    Av_depth[p$ix[[1]][1]:p$ix[[1]][length(p$ix[[1]])]] <- Av_depth[p$ix[[1]][1]:p$ix[[1]][length(p$ix[[1]])]] + 0.1 * p$bottom
    Av_depth[p$ix[[3]][1]:p$ix[[3]][length(p$ix[[3]])]] <- Av_depth[p$ix[[3]][1]:p$ix[[3]][length(p$ix[[3]])]] - 0.1 * p$bottom
    
    # Create flux from interaction: 
    # Coordinates for lines between points
    # Select major interactions and scale sizes:
    # Set color palette 
    
    p$SpId <- c('smallPel','mesoPel','largePel', 'bathyPel', 'demersals')
    SpId <- c("smallZoo", "largeZoo", "smallBenthos", "largeBenthos", 
              rep(p$SpId[1], length(p$ix[[1]])),
              rep(p$SpId[2], length(p$ix[[2]])),
              rep(p$SpId[3], length(p$ix[[3]])),
              rep(p$SpId[4], length(p$ix[[4]])),
              rep(p$SpId[5], length(p$ix[[5]])))
    
  }
  
  # Marker size depends on biomass: 
  # Using real biomass yields bubles with too many orders of magnitude difference
  # Thus we group them by quantiles
  Msize <- Bi / max(Bi)
  Msize[Msize == 0] <- NA
  idxM <- quantile(Msize, prob = c(0.2, 0.4, 0.6, 0.8), na.rm = T) # get quantiles
  
  # Specify buble size for each quantile:
  Msize[Msize >= idxM[4] & !is.na(Msize)] <- 20
  Msize[Msize >= idxM[3] & Msize < idxM[4] & !is.na(Msize)] <- 15
  Msize[Msize >= idxM[2] & Msize < idxM[3] & !is.na(Msize)] <- 8
  Msize[Msize >= idxM[1] & Msize < idxM[2] & !is.na(Msize)] <- 3
  Msize[Msize < idxM[1] & !is.na(Msize)] <- .8
  
  # Create line width: 
  Mat <- rep(0, ngroup) 
  Mat[Bi != 0] <- 1
  Theta <- t(t(p$theta) * Bi) * Mat # flux equal the rate * the prey biomass (* 0 if pred <- 0)
  Theta <- c(Theta) 
  threshold <- 0.05 # min(tail(sort(Theta), 100)) # Alternatively, use 100 strongest relations regardless of absolute value of the threshold
  indx <- which(Theta >= threshold) # takes the x highest values of theta
  
  
  # Set values of each coordinate and put them together:
  coord_1 <- data.frame(index = 1:p$nStages^2,
                        mc = rep(p$mc[1:p$nStages], p$nStages), 
                        depth = rep(Av_depth[1:p$nStages], p$nStages), 
                        SpId = rep(SpId, p$nStages),
                        Msize = rep(Msize, p$nStages), 
                        LineWdth = Theta/max(Theta),
                        Alpha = Theta/max(Theta))
  
  coord_2 <- data.frame(index = 1:p$nStages^2, # Notice that here repetition ys grouped by "each" to change order
                        mc = rep(p$mc[1:p$nStages], each = p$nStages), 
                        depth = rep(Av_depth[1:p$nStages], each = p$nStages), 
                        SpId = rep(SpId, each = p$nStages),
                        Msize = rep(Msize, each = p$nStages),
                        LineWdth = Theta/max(Theta),
                        Alpha = Theta/max(Theta))
  
  df <- rbind(coord_1, coord_2)
  
  df <- df %>% filter(index %in% indx) %>%
    arrange(desc(Msize))
  
  if (length(p$ix)==3){
    p <- ggplot(data = df) +
      geom_line(aes(x = mc, y = depth, group = index, size = LineWdth, color = SpId, alpha = Alpha), show.legend = F) +
      geom_point(aes(x = mc, y = depth, color = SpId, size = Msize)) +
      scale_color_manual(values = p$my_palette[attr(p$my_palette, "names") %in% df$SpId], 
                         labels = p$my_names[attr(p$my_palette, "names") %in% df$SpId]) +
      scale_size_continuous(range = c(1, 15)) +
      scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x))) +
      annotation_logticks(sides = "b") +
      labs(x ="Weight (grams)", y = "", color = "Group") +
      theme_base() + 
      guides(size = "none") +
      theme(legend.position = "bottom",
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
  }
  
  if (length(p$ix)==5){
    p <- ggplot(data = df) +
      geom_line(aes(x = mc, y = depth, group = index, size = LineWdth, color = SpId, alpha = Alpha), show.legend = F) +
      geom_point(aes(x = mc, y = depth, color = SpId, size = Msize)) +
      scale_color_manual(values = p$my_palette[attr(p$my_palette, "names") %in% df$SpId], 
                         labels = p$my_names[attr(p$my_palette, "names") %in% df$SpId]) +
      scale_size_continuous(range = c(1, 15)) +
      scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x))) +
      annotation_logticks(sides = "b") +
      labs(x ="Weight (grams)", y = "Depth (m)", color = "Group") +
      theme_base() + 
      guides(size = "none") +
      theme(legend.position = "bottom")
    
  }
  # ggsave("plot_network.png", p, height = 45 , width = 80, units = "mm", scale = 3)
  
  return(p)
}


#-------------------------------------------------------------------------------
# Diet plot 
# Revised Based on work of Daniel Ottmann Riera and Solenne Roux
# u: sim$u simpleOutput=TRUE
#-------------------------------------------------------------------------------

plotDiet <- function(p, u) {
  p$nstage <-lengths <- max(sapply(p$ix, length)) #maximum number of stages for one group
  biomass <- u#sim [,2:(p$ix[[length(p$ix)]][length(p$ix[[length(p$ix)]])]+1)]
  Bin <- round(0.6 * nrow(biomass), digits = 0)
  biomassend <- colMeans(biomass[Bin:nrow(biomass),])
  biomassstage <- p$ixFish[length(p$ixFish)]
  biomasssmall <- p$nstage - round(2/3*p$nstage, digits = 0)
  Enc = p$V * (p$theta %*% biomassend)
  f   = Enc / (p$Cmax + Enc)
  f[is.na(f)] = 0  
  #f <- calcEncounter(biomassend, p)$f
  
  bom <- t(t(p$theta[5:biomassstage, ]) * colMeans(biomass[Bin:nrow(biomass),])) 
  fbom <- f[5:biomassstage] / rowSums(bom)
  output <- bom * fbom
  
  if (length(p$ix)==5){
    
    p$SpId <- c('smallPel','mesoPel','largePel', 'bathyPel', 'demersals')
    SpId <- c("smallZoo", "largeZoo", "smallBenthos", "largeBenthos", 
              rep(p$SpId[1], length(p$ix[[1]])),
              rep(p$SpId[2], length(p$ix[[2]])),
              rep(p$SpId[3], length(p$ix[[3]])),
              rep(p$SpId[4], length(p$ix[[4]])),
              rep(p$SpId[5], length(p$ix[[5]])))
    
    p$RSpName <- c("Small zooplankton", "Big zooplankton", "Benthos", "Small pelagics",
                   "Mesopelagics", "Large pelagics", "Bathypelagics", "Demersals")
    
  } else {
    
    p$SpId <- c('smallPel','largePel', 'demersals')
    SpId <- c("smallZoo", "largeZoo", "smallBenthos", "largeBenthos", 
              rep(p$SpId[1], length(p$ix[[1]])),
              rep(p$SpId[2], length(p$ix[[2]])),
              rep(p$SpId[3], length(p$ix[[3]])))
    
    p$RSpName <- c("Small zooplankton", "Large zooplankton", "Benthos", "Small pelagics",
                   "Large pelagics", "Demersals")   
    
  }
  
  
  # small pelagics: ---------------------------------------------------------
  small_pel <- output[(p$ix[[1]][1] - length(p$ixR)):(p$ix[[1]][length(p$ix[[1]])] - length(p$ixR)), ] 
  small_pel <- t(rbind(small_pel, matrix(0, biomasssmall, biomassstage)))
  small_pel <- data.frame(val = c(small_pel), 
                          stage = rep(1:p$nstage, each = nrow(small_pel)), 
                          SpId = as.factor(rep(SpId, p$nstage)))
  
  p1 <- ggplot(data = small_pel) +
    geom_bar(aes(x = stage, y = val, fill = SpId), stat = "identity") +
    scale_fill_manual(values = p$my_palette) +
    theme_base() +
    labs(x = "", y = "Fraction of stomach", fill = "Prey group", title = p$RSpName[4]) +
    theme(legend.position = "none")
  
  
  if (length(p$ix)==5){
    
    # Large pelagics: --------------------------------------------------------- 
    large_pel <- t(output[(p$ix[[3]][1] - length(p$ixR)):(p$ix[[3]][length(p$ix[[3]])] - length(p$ixR)), ])
    large_pel <- data.frame(val = c(large_pel), 
                            stage = rep(1:p$nstage, each = nrow(large_pel)), 
                            SpId = as.factor(rep(SpId, p$nstage)))
    
    p2 <- ggplot(data = large_pel) +
      geom_bar(aes(x = stage, y = val, fill = SpId), stat = "identity") +
      scale_fill_manual(values = p$my_palette) +
      theme_base() +
      labs(x ="size-class", y = "", fill = "Prey group", title = p$RSpName[6]) +
      theme(legend.position = "none") 
    
    
    # Demersal: ---------------------------------------------------------------
    demers <- t(output[(p$ix[[5]][1] - length(p$ixR)):(p$ix[[5]][length(p$ix[[5]])] - length(p$ixR)), ])
    demers <- data.frame(val = c(demers), 
                         stage = rep(1:p$nstage, each = nrow(demers)), 
                         SpId = as.factor(rep(SpId, p$nstage)))
    
    p3 <- ggplot(data = demers) +
      geom_bar(aes(x = stage, y = val, fill = SpId), stat = "identity") +
      scale_fill_manual(values = p$my_palette) +
      theme_base() +
      labs(x ="size-class", y = "", fill = "Prey group", title = p$RSpName[8]) +
      theme(legend.position = "none")
    
    
    # Mesopelagics: ----------------------------------------------------------- 
    meso_pel <- output[(p$ix[[2]][1] - length(p$ixR)):(p$ix[[2]][length(p$ix[[2]])] - length(p$ixR)), ]
    meso_pel <- t(rbind(meso_pel, matrix(0, biomasssmall, biomassstage)))
    meso_pel <- data.frame(val = c(meso_pel), 
                           stage = rep(1:p$nstage, each = nrow(meso_pel)), 
                           SpId = as.factor(rep(SpId, p$nstage)))
    
    p5 <- ggplot(data = meso_pel) +
      geom_bar(aes(x = stage, y = val, fill = SpId), stat = "identity") +
      scale_fill_manual(values = p$my_palette) +
      theme_base() +
      labs(x ="size-class", y = "", fill = "Prey group", title = p$RSpName[5]) +
      theme(legend.position = "none") 
    
    # Bathypelagics: ---------------------------------------------------------- 
    bathy_pel <- t(output[(p$ix[[4]][1] - length(p$ixR)):(p$ix[[4]][length(p$ix[[4]])] - length(p$ixR)), ])
    bathy_pel <- data.frame(val = c(bathy_pel), 
                            stage = rep(1:p$nstage, each = nrow(bathy_pel)), 
                            SpId = as.factor(rep(SpId, p$nstage)))
    
    p6 <- ggplot(data = bathy_pel) +
      geom_bar(aes(x = stage, y = val, fill = SpId), stat = "identity") +
      scale_fill_manual(values = p$my_palette) +
      theme_base() +
      labs(x ="size-class", y = "", fill = "Prey group", title = p$RSpName[7]) +
      theme(legend.position = "none") 
    
    
  } else {
    
    # Large pelagics: --------------------------------------------------------- 
    large_pel <- t(output[(p$ix[[2]][1] - length(p$ixR)):(p$ix[[2]][length(p$ix[[2]])] - length(p$ixR)), ])
    large_pel <- data.frame(val = c(large_pel), 
                            stage = rep(1:p$nstage, each = nrow(large_pel)), 
                            SpId = as.factor(rep(SpId, p$nstage)))
    
    p2 <- ggplot(data = large_pel) +
      geom_bar(aes(x = stage, y = val, fill = SpId), stat = "identity") +
      scale_fill_manual(values = p$my_palette) +
      theme_base() +
      labs(x ="size-class", y = "", fill = "Prey group", title = p$RSpName[5]) +
      theme(legend.position = "none") 
    
    
    # Demersal: ---------------------------------------------------------------
    demers <- t(output[(p$ix[[3]][1] - length(p$ixR)):(p$ix[[3]][length(p$ix[[3]])] - length(p$ixR)), ])
    demers <- data.frame(val = c(demers), 
                         stage = rep(1:p$nstage, each = nrow(demers)), 
                         SpId = as.factor(rep(SpId, p$nstage)))
    
    p3 <- ggplot(data = demers) +
      geom_bar(aes(x = stage, y = val, fill = SpId), stat = "identity") +
      scale_fill_manual(values = p$my_palette) +
      theme_base() +
      labs(x ="size-class", y = "", fill = "Prey group", title = p$RSpName[6]) +
      theme(legend.position = "none")
    
  }
  
  
  # Legend: ----------------------------------------------------------------- 
  legend_colors <- data.frame(val = 0, stage = 0, SpId = unique(large_pel$SpId))
  
  p7 <- ggplot(data = legend_colors) +
    geom_bar(aes(x = stage, y = val, fill = SpId), stat = "identity") +
    scale_fill_manual(values = p$my_palette[attr(p$my_palette, "names") %in% legend_colors$SpId]) +
    labs(fill = "Prey group") +
    theme_void() +
    theme(legend.position = c(.45,.4),
          legend.direction = "horizontal",
          legend.text = element_text(size = 12))+
    guides(fill = guide_legend(
      nrow = 3,  # Spécifiez le nombre de lignes pour la légende
      byrow = TRUE,  # Indique que les étiquettes de légende doivent être disposées par ligne
      title.position = "top",  # Positionne le titre de la légende en haut
    ))
  
  
  # Put all panels together:
  
  if (length(p$ix)==5){
    
    if(p$bottom > p$mesop) {
    plots <- p1 + p2 + p3 + p5 + p6 +p7 & 
      theme(plot.background = element_blank())
    }else{
    plots <- p1 + p2 + p3 +p7 & 
      theme(plot.background = element_blank())
    }
    
  } else {
    
    plots <- p1 + p2 + p3 +p7 & 
      theme(plot.background = element_blank())
  }
  
  
  
  return(plots)
}


#-------------------------------------------------------------------------------
# Make a basic run:
#-------------------------------------------------------------------------------

baserun = function(USEdll=TRUE) {
  p = setupBasic()
  sim = simulate(p,tEnd = 100,USEdll)
  plotSimulation(sim)
  return(sim)
}


defaultplot <- function(
  mfcol=c(1,1), 
  oma=c(0, 0.4, omargin, omargin), # Outer margins (bottom, left, top, right)
  mar=c(2.1,2.3,0,0), # Margins
  ps=10, # Point size
  tcl=ticklength, # Tick mark length
  mgp=c(1.1,0,0), ...)  # Margins for axis title, axis labels, and axis line
{
  par(ps=ps)
  par(mfcol=mfcol, oma=oma, mar=mar, tcl=tcl, mgp=mgp, cex=cex, cex.axis=cexaxis, ...) 
  assign("iPlot", 1, envir = .GlobalEnv)
}

defaultplothorizontal <- function(
  nPanels=1,
  mfcol=c(1,nPanels), 
  oma=c(2.1, 2.7, omargin, omargin), # Outer margins (bottom, left, top, right)
  mar=c(0,0,0,0.3), # Margins
#  oma=c(0,0.4, omargin, omargin), # Outer margins (bottom, left, top, right)
#  mar=c(2.7,2.5,0,0.3), # Margins
  ps=10, # Point size
  tcl=ticklength, # Tick mark length
  mgp=c(1.1,0,0), ...)  # Margins for axis title, axis labels, and axis line
{
  par(mfcol=mfcol, oma=oma, mar=mar, ps=ps, tcl=tcl, mgp=mgp, cex=cex, cex.axis=cexaxis, ...) 
  assign("iPlot", 1, envir = .GlobalEnv)
}

defaultplotvertical <- function(
  nPanels=1,
  mfcol=c(nPanels,1), 
  oma=c(2.1, 2.5, omargin, omargin), # Outer margins (bottom, left, top, right)
  mar=c(0.3,0,0,0), # Margins
  ps=10, # Point size
  tcl=ticklength, # Tick mark length
  mgp=c(2,0,0), ...)  # Margins for axis title, axis labels, and axis line
{
  par(mfcol=mfcol, oma=oma, mar=mar, ps=ps, tcl=tcl, mgp=mgp, cex=cex, cex.axis=cexaxis, ...) 
  assign("iPlot", 1, envir = .GlobalEnv)
}

test <- function() {
  defaultplot()
  defaultpanel(xlab="x", ylab="y",xlim=c(0.1,10), ylim=c(0.1,10))
  points(1,1)
  
  defaultplot()
  loglogpanel(xlab="x", ylab="y",xlim=c(0.1,10), ylim=c(0.1,10))
  points(1,1)

  defaultplothorizontal(nPanels = 3)
  loglogpanel(xlab="x", ylab="y",xlim=c(0.1,10), ylim=c(0.1,10))
  points(1,1)
  loglogpanel(xlab="x", ylab="",xlim=c(0.1,10), ylim=c(0.1,10), yaxis = FALSE)
  points(1,1)
  loglogpanel(xlab="x", ylab="",xlim=c(0.1,10), ylim=c(0.1,10), yaxis = FALSE)
  points(1,1)
  
  defaultplothorizontal(nPanels = 2)
  defaultpanel(xlab="x", ylab="y",xlim=c(0.1,10), ylim=c(0.1,10))
  points(1,1)
  defaultpanel(xlab="x", ylab="",xlim=c(0.1,10), ylim=c(0.1,10), yaxis = FALSE)
  points(1,1)
  
  defaultplotvertical(nPanels = 2)
  loglogpanel(xlab="", ylab="y",xlim=c(0.1,10), ylim=c(0.1,10), xaxis=FALSE)
  points(1,1)
  loglogpanel(xlab="x", ylab="y",xlim=c(0.1,10), ylim=c(0.1,10))
  points(1,1)
  
}

defaultpanel <- function(xlim, ylim, 
                         xlab='', ylab='', 
                         xaxis=TRUE, yaxis=TRUE, label=FALSE, new=FALSE,
                         bty="o", xaxs="r",
                         main="") {
  plot(1, type='n', 
       ylim=range(ylim[!is.na(ylim)]), 
       xlim=range(xlim[!is.na(xlim)]), axes=FALSE, xlab='', ylab='', par(new=new),
       bty=bty, xaxs=xaxs, main=main)
  if (xaxis)
    mtext(side=bottom, line=1, TeX(xlab), cex=par()$cex)
  if (yaxis)
    mtext(side=left, line=1, TeX(ylab), cex=par()$cex)
  if (label) 
    makepanellabel()
  if (xaxis)
    axis(bottom, labels=xaxis, cex.axis=par()$cex.axis, lwd=axis.lwd, lwd.ticks=axis.lwd)
  if (yaxis)
    axis(left, labels=yaxis, cex.axis=par()$cex.axis, lwd=axis.lwd, lwd.ticks=axis.lwd)
  if (bty == "o")
    box(lwd=axis.lwd)

}

semilogxpanel <- function(xlim, ylim, xlab='', ylab='', 
                          xaxis=TRUE, yaxis=TRUE, label=FALSE, new=FALSE, 
                          bty="o", xaxs="r") {
  ylim <- range(na.omit(ylim))
  xlim <- range(na.omit(xlim))
  plot(1, type='n', log='x',
       ylim=ylim, 
       xlim=xlim, axes=FALSE, xlab='',ylab='', xaxs=xaxs, par(new=new),
       xaxs="i", yaxs="i")

  if (xaxis)
    mtext(side=bottom, line=1, TeX(xlab))
  if (yaxis)
    mtext(side=left, line=1, TeX(ylab))
  if (label) 
    makepanellabel()
#  if (xaxis)
  logaxes(bottom, lim=xlim, labels=xaxis)
#  if (yaxis)
  axis(left, labels=yaxis, lwd=axis.lwd, lwd.ticks=axis.lwd)
  if (bty=="o")
    box(lwd=axis.lwd)
}

semilogypanel <- function(xlim, ylim, xlab='', ylab='', 
                          xaxis=TRUE, yaxis=TRUE, label=FALSE,
                          bty="o") {
  plot(1, type='n', log='y',
       ylim=range(ylim[!is.na(ylim)]), 
       xlim=range(xlim[!is.na(xlim)]), axes=FALSE, xlab='',ylab='',yaxs="i",
       bty=bty)
  if (xaxis)
    mtext(side=bottom, line=1, TeX(xlab))
  if (yaxis) 
    mtext(side=left, line=1.5, TeX(ylab))
  if (label) 
    makepanellabel()
  if (xaxis)
    axis(bottom, labels=xaxis, lwd=axis.lwd, lwd.ticks=axis.lwd)
  logaxes(left, lim=ylim, labels=yaxis, drawaxis=yaxis)
  if (bty=="o")
    box(lwd=axis.lwd)
}

loglogpanel <- function(xlim, ylim, xlab='', ylab='', 
                        xaxis=TRUE, yaxis=TRUE, label=FALSE,
                        bExponential=TRUE, new=FALSE, powx=NA, powy=NA, xaxs="r") {
  ylim <- range(na.omit(ylim))
  xlim <- range(na.omit(xlim))
  plot(1, type='n', log='xy',
       ylim=ylim, 
       xlim=xlim, axes=FALSE, xlab='',ylab='', par(new=new), xaxs=xaxs,
       xaxs="i",yaxs="i")
  mtext(side=bottom, line=1.1, TeX(xlab))
  mtext(side=left, line=1.5, TeX(ylab))
  if (label) 
    makepanellabel()
  logaxes(bottom, lim=xlim, bExponential = bExponential, labels=xaxis, pow=powx)
  logaxes(left, lim=ylim, bExponential = bExponential, labels=yaxis, pow=powy)
  box(lwd=axis.lwd)
}


## side:  1 for x-axis, 2 for y-axis (3, 4 secondary x and y axes)
## labels: logical, if TRUE, adds them
## las: integer in [0,1,2,3] see ?par
logaxes <- function(side = bottom, 
                    lim, pow = NA, drawaxis=TRUE,
                    labels=TRUE, las = 1, col = 1, bExponential=TRUE) {
  poww <- pow
  if (is.na(pow[1]))
    poww = ceiling(min(log10(lim))):floor(max(log10(lim)))
  
  #axis(side, at = 10^poww, lwd=0, labels = FALSE, tcl=ticklength, lwd.ticks = axis.lwd, col = col)

  if (labels) 
    axis(side, at = 10^poww, lwd=0, labels = FALSE, tcl=-0., lwd.ticks = axis.lwd)

  if (side==bottom)
    line <- 0.2/(par()$oma[bottom] + par()$mar[bottom])
  else
    line <- 0.6/(par()$oma[left] + par()$mar[left])
  
  if (labels) {
    for(i in poww) {
      if (bExponential)
        mtext(side = side, at = 10^i, text=bquote(10^.(i)), line = line, 
              las = las, cex=cexaxis)
      else
        mtext(side = side, at = 10^i, text=10^i, line = 0.2, 
              las = las, cex=par()$cex.axis)
    }
  }
  # Minor ticks:
  if ((is.na(pow[1])) & drawaxis) {
    at <- as.vector(sapply(seq(range(poww)[1]-1, range(poww)[2]), function(p) (2:9)*10^p))
    axis(side, at = at, labels = FALSE, tcl=0.5*ticklength, lwd=0, lwd.ticks=axis.lwd, col = col)
  }
} 

makepanellabel <- function(line=-1.1) {
  mtext(letters[iPlot], side=top, line=line, adj=0.05)
  assign("iPlot", iPlot+1, envir = .GlobalEnv)
}

hline <- function(y=0, lty='dotted', lwd=1) {
  if (par("xlog"))
    lines(x=10^par("usr")[1:2], y=y*c(1,1), lty=lty, lwd=lwd)
  else
    lines(x=par("usr")[1:2], y=y*c(1,1), lty=lty, lwd=lwd)
}

vline <- function(x=0, lty='dotted', col="black") {
  if (par("ylog"))
    lines(x=x*c(1,1), y=10^par("usr")[3:4], lty=lty, col=col)
  else
    lines(x=x*c(1,1), y=par("usr")[3:4], lty=lty, col=col)
}

pdfplot <- function(filename, FUN, ..., width=singlewidth, height=height) {
  pdf(filename, width=width, height=height, useDingbats=FALSE)
  FUN(...)
  dev.off()  
}

addEpsPicture <- function(sName, x, y, width=1) {
  # Convert the picture
  sOutName = paste(sName,'.xml',sep='')
  PostScriptTrace(sName, sOutName)
  # Add it
  grid.picture(readPicture(sOutName), x=x, y=y, width=width)
}

ribbon <- function(x,ymin=NA,ymax,col=lightgrey) {
  x <- c(x, x[seq(length(x),1,by = -1)])
  polygon(x, c(ymin, ymax[seq(length(ymax),1,by = -1)]),
          col=col, border=NA)
}

tightaxes <- function() {
  par(xaxs="i", yaxs="i")
}
