#
# Background tools for plotting
#
require(latex2exp)

singlewidth <- 8/2.54 # panel width in cm
doublewidth <- 13/2.54
height <- 8/1.6/2.54 + 2.5*0.138889

if(.Platform$OS.type=="windows") {
  quartz<-function(width, height) 
    windows(width=width, height=height)
}

# Base R plotting tools =====================================================

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
