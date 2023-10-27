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
# Makes a basic plot of the adult biomass (SSB) of all functional groups
# as a function time.
#-------------------------------------------------------------------------------

plotSSBtime = function(sim, bPlot=TRUE) {
  if (bPlot)
    defaultplot()
  p = sim$p
  
  semilogypanel(xlim=sim$t, ylim=sim$SSB+1e-10,
                xlab="Time (yr)", ylab = "SSB (gww  m-2)")
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
              xlab = "Mass (gww)", ylab="Biomass (gww m-2)")
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
