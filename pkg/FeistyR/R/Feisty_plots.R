#===============================================================================
# Functions to creating plots to the Feisty model
#
# Slightly rewritten by Solenne Roux, based on code from Ken H. Andersen
#===============================================================================
library(tidyverse) #package which contains the functions to edit dataframes and generate plots
library(scales) #package for the functions to edit visualisation 
library(ggthemes) #package for the function "theme_base"
library(patchwork) #to combine multiple panels in a single plot
#-------------------------------------------------------------------------------
# Makes a basic plot of the adult biomass (SSB) of all functional groups
# as a function time.
#-------------------------------------------------------------------------------

plotSSBtime = function(sim, bPlot=TRUE) {
  if (bPlot)
    defaultplot()
  p = sim$p
  
  semilogypanel(xlim=sim$t, ylim=sim$SSB+1e-10,
                xlab="Time (yr)", ylab = "SSB (gww  m-3)")
  #
  # Plot fish
  #
  for (i in 1:p$nGroups)
    lines(sim$t, sim$SSB[,i])
  #
  # Plot resources:
  #
  for (i in p$ixR)
    lines(sim$t, sim$R[,i], col='blue')
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
  loglogpanel(xlim = xlim, ylim=rates$g+1e-10,
              ylab="Growth rate (1/year)", xlab="-", xaxis = FALSE)
  for (i in 1:p$nGroups) {
    lines(p$mc[p$ix[[i]]], rates$g[p$ix[[i]]-length(p$ixR)], lwd=i, col='black')
  }
  
  #
  # Mortalities:
  #
  loglogpanel(xlim=xlim, ylim=rates$mortpred+1e-10,
              xlab="-", ylab="mort (1/year)", xaxis = FALSE)
  for (i in 1:p$nGroups) {
    lines(p$mc[p$ix[[i]]], rates$mortpred[p$ix[[i]]], lwd=i, col='red')
    lines(p$mc[p$ix[[i]]], p$mortF[p$ix[[i]]], lwd=i, col='blue')
  }
  hline(p$mort0[p$ix[[i]]][1])
  
  legend(x='bottomleft',
         legend=c('Predation','Background'),
         lty=c(1,dotted),
         col=c('red','black'),
         bty='n')
  #
  # Feeding level
  # 
  semilogxpanel(xlim = xlim, ylim=c(0,1),
                ylab="Feeding level, f", xlab="Mass (gww)")
  for (i in 1:p$nGroups) {
    ix = p$ix[[i]]
    lines(p$mc[ix], rates$f[ix], lwd=i)
    lines(p$mc[ix], p$metabolism[ix]/(p$epsAssim*p$Cmax[ix]), lwd=i, lty=dotted) # Critical feeding level
  }
  
  return(rates)
}

#-------------------------------------------------------------------------------
# Plot the biomasses of all groups
#-------------------------------------------------------------------------------

plotSpectra = function(sim, iTime=sim$nTime, bPlot=TRUE) {
  if (bPlot)
    defaultplot()
  p = sim$p
  
  loglogpanel(xlim=p$mc[p$ixFish], ylim=pmax(1e-10,sim$B[iTime,]),
              xlab = "Mass (gww)", ylab="Biomass (gww m-3)")
  for (i in 1:p$nGroups) {
    lines(p$mc[p$ix[[i]]], sim$B[iTime, p$ix[[i]]-p$ixFish[1]+1], lwd=i)
  }
}

#-------------------------------------------------------------------------------
# Make 4 panels of simulation
#-------------------------------------------------------------------------------

plotSimulation = function(sim) {
  defaultplot(mfcol=c(5,1))
  plotSSBtime(sim,bPlot=FALSE)
  plotSpectra(sim, bPlot=FALSE)
  rates = plotRates(sim$p, u=c( sim$R[sim$nTime,], sim$B[sim$nTime,]),bPlot=FALSE)
  
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
# Bubbles plot :
#-------------------------------------------------------------------------------

plot_network <- function(p, sim) {
  
  # Number of groups and biomass:
  ngroup <- p$ix[[length(p$ix)]][length(p$ix[[length(p$ix)]])] #ressources (4) + fish 
  biomass <- sim[,2:(p$ix[[length(p$ix)]][length(p$ix[[length(p$ix)]])]+1)] #ressources + poissons car se mangent entre eux (+1 car première colonne = temps)
  
  #Average of the biomass : 
  Bi <- colMeans(biomass[(nrow(biomass) -50):nrow(biomass),])
  
  my_palette_1 <- c("smallZoo" = "#BBDEFB",
                    "largeZoo" = "#2196F3",
                    "smallBenthos" = "#0D47A1",
                    "largeBenthos" = "#795548",
                    "smallPel" = "#F57C0D",
                    "largePel" = "#FFEE58",
                    "Demersals" =  "#F9A825")
  
  my_names_1 <- c("smallZoo" = "Small zooplankton",
                  "largeZoo" = "Large zooplankton",
                  "smallBenthos" = "Small Benthos",
                  "largeBenthos" = "Large Benthos",
                  "smallPel" = "Small pelagics",
                  "largePel" = "Large pelagics",
                  "Demersals" =  "Demersals")
  
  my_palette_2 <- c("smallZoo" = "#BBDEFB",
                    "largeZoo" = "#2196F3",
                    "smallBenthos" = "#0D47A1",
                    "largeBenthos" = "#795548",
                    "smallPel" = "#F57C0D",
                    "largePel" = "#FFEE58",
                    "Demersals" =  "#F9A825")
  
  my_names_2 <- c("smallZoo" = "Small zooplankton",
                  "largeZoo" = "Large zooplankton",
                  "smallBenthos" = "Small Benthos",
                  "largeBenthos" = "Large Benthos",
                  "smallPel" = "Small pelagics",
                  "largePel" = "Large pelagics",
                  "Demersals" =  "Demersals")
  
  
  if (p$setup == "setupBasic"){
    Av_depth <- c(-1,-1,-4,-4,0,0,-2,-2,-2,-3,-3,-3)
    
    p$SpId <- c('smallPel','largePel', 'Demersals')
    SpId <- c("smallZoo", "largeZoo", "smallBenthos", "largeBenthos", 
              rep(p$SpId[1], length(p$ix[[1]])),
              rep(p$SpId[2], length(p$ix[[2]])),
              rep(p$SpId[3], length(p$ix[[3]])))
    
    p$my_palette <- my_palette_1
    
    p$my_names <- my_names_1
  }  
  
  if (p$setup == "setupBasic2"){
    Av_depth <- c(-1,-1,-4,-4,0,0,0,0,0,0,-2,-2,-2,-2,-2,-2,-2,-2,-2,-3,-3,-3,-3,-3,-3,-3,-3,-3)
    
    p$SpId <- c('smallPel','largePel', 'Demersals')
    SpId <- c("smallZoo", "largeZoo", "smallBenthos", "largeBenthos", 
              rep(p$SpId[1], length(p$ix[[1]])),
              rep(p$SpId[2], length(p$ix[[2]])),
              rep(p$SpId[3], length(p$ix[[3]])))
    
    p$my_palette <- my_palette_1
    
    p$my_names <- my_names_1
  }  
  
  
  
  if (p$setup == "setupVertical"){
    
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
    
    p$SpId <- c('smallPel','mesoPel','largePel', 'bathyPel', 'Demersals')
    SpId <- c("smallZoo", "largeZoo", "smallBenthos", "largeBenthos", 
              rep(p$SpId[1], length(p$ix[[1]])),
              rep(p$SpId[2], length(p$ix[[2]])),
              rep(p$SpId[3], length(p$ix[[3]])),
              rep(p$SpId[4], length(p$ix[[4]])),
              rep(p$SpId[5], length(p$ix[[5]])))
    
    p$my_palette <- my_palette_2
    
    p$my_names <- my_names_2
    
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
# Bubbles plot :
#-------------------------------------------------------------------------------

plotdiet <- function(p, sim) {
  p$nstage <-lengths <- max(sapply(p$ix, length)) #maximum number of stages for one group
  biomass <- sim [,2:(p$ix[[length(p$ix)]][length(p$ix[[length(p$ix)]])]+1)]
  Bin <- floor(0.8 * nrow(biomass))
  biomassend <- colMeans(biomass[Bin:nrow(biomass),])
  biomassstage <- p$ixFish[length(p$ixFish)]
  biomasssmall <- p$nstage - p$nstage * 2/3
  f <- calcEncounter(biomassend, p)$f
  
  bom <- t(t(p$theta[5:biomassstage, ]) * colMeans(biomass[Bin:nrow(biomass),])) 
  fbom <- f[5:biomassstage] / rowSums(bom)
  output <- bom * fbom
  
  if (length(p$ix)==5){
    
    p$SpId <- c('smallPel','mesoPel','largePel', 'bathyPel', 'Demersals')
    SpId <- c("smallZoo", "largeZoo", "smallBenthos", "largeBenthos", 
              rep(p$SpId[1], length(p$ix[[1]])),
              rep(p$SpId[2], length(p$ix[[2]])),
              rep(p$SpId[3], length(p$ix[[3]])),
              rep(p$SpId[4], length(p$ix[[4]])),
              rep(p$SpId[5], length(p$ix[[5]])))
    
    p$my_palette <- c("smallZoo" = "#BBDEFB",
                      "largeZoo" = "#2196F3",
                      "smallBenthos" = "#0D47A1",
                      "largeBenthos" = "#795548",
                      "smallPel" = "#F57C0D",
                      "mesoPel" = "#9E9E9E",
                      "largePel" = "#FFEE58",
                      "bathyPel" = "#000000", 
                      "Demersals" =  "#F9A825")
    
    p$RSpName <- c("Small zooplankton", "Big zooplankton", "Benthos", "Small pelagics",
                   "Mesopelagics", "Large pelagics", "Bathypelagics", "Demersals")
    
  } else {
    
    p$SpId <- c('smallPel','largePel', 'Demersals')
    SpId <- c("smallZoo", "largeZoo", "smallBenthos", "largeBenthos", 
              rep(p$SpId[1], length(p$ix[[1]])),
              rep(p$SpId[2], length(p$ix[[2]])),
              rep(p$SpId[3], length(p$ix[[3]])))
    
    p$my_palette <- c("smallZoo" = "#BBDEFB",
                      "largeZoo" = "#2196F3",
                      "smallBenthos" = "#0D47A1",
                      "largeBenthos" = "#795548",
                      "smallPel" = "#F57C0D",
                      "largePel" = "#FFEE58",
                      "Demersals" =  "#F9A825")
    
    p$RSpName <- c("Small zooplankton", "Big zooplankton", "Benthos", "Small pelagics",
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
    plots <- p1 + p2 + p3 + p5 + p6 +p7 & 
      theme(plot.background = element_blank())
    
  } else {
    
    plots <- p1 + p2 + p3 +p7 & 
      theme(plot.background = element_blank())
  }
  
  
  
  return(plots)
}


