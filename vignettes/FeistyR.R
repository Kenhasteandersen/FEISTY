## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
if (!require(FEISTY)) devtools::install_github("Kenhasteandersen/FEISTY@R-Package")
library(FEISTY)
palette("ggplot2")  # ggplot2-style palette

## -----------------------------------------------------------------------------
p <- setupBasic()
knitr::kable(p$resources, digits=2)
knitr::kable(p$fishes, digits=2)
knitr::kable(p$groups, digits=2)

## -----------------------------------------------------------------------------
p2 <- paramAddPhysiology(p, ac = 40, am = 8, ae=140)

p3 <- paramAddPhysiology(p, ac = 10, am = 2, ae=35)

## ---- fig.width=8, fig.height=8-----------------------------------------------
out1 <- simulateFEISTY(p=p,  times=seq(0, 200, length.out=1000),bCust=T)
out2 <- simulateFEISTY(p=p2, times=seq(0, 200, length.out=1000),bCust=T)
out3 <- simulateFEISTY(p=p3, times=seq(0, 200, length.out=1000),bCust=T)
# plot(out1, out2, out3, which=5:12, lty=1, lwd=2, subset=time>180)
# plot(out1, out2, out3, which=c("smallZoo", "largeZoo", "smallBenthos", 
#   "totBiomass.smallPel", "totBiomass.largePel",  "totBiomass.Demersals"), 
#   lty=1, lwd=2, subset=time>180)

