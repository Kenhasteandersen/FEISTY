## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE) # for printing the examples
library(FEISTY)
library(ggplot2)
palette("ggplot2")  # ggplot2-style palette

## -----------------------------------------------------------------------------
p <- setupBasic(szprod = 50,   # small mesozooplankton production
                lzprod = 50,   # large mesozooplankton production
                bprodin  = 8,  # benthos production
                depth  = 150,  # water column depth [m]
                Tp     = 17,   # pelagic layer averaged temperature [Celsius]
                Tb     = 12)   # sea floor temperature [Celsius]

## -----------------------------------------------------------------------------
sim <- simulateFEISTY(p=p,  times=seq(0, 500, length.out=500), USEdll = T)

## ----fig.width=8, fig.height=8------------------------------------------------
plotSimulation(sim)

## ----fig.width=8, fig.height=7------------------------------------------------
plotBiomasstime(sim)

## ----fig.width=8, fig.height=7------------------------------------------------
plotRates(sim)

## -----------------------------------------------------------------------------
p=setupVertical2(szprod= 120,
                 lzprod = 120,
                 dfpho=200, 
                 depth = 700, 
                 nStages = 9, 
                 F=0) # no fishing for all size classes
names(p$mortF)=p$stagenames
df=data.frame(mortF_original=c(p$mortF[p$ix[[1]]],p$mortF[p$ix[[5]]]))
# assign 0.2/year as the maximum fishing mortality to small pelagic fish
p=setFishing(p=p,F=0.2,etaF=0.05,groupidx=c(1))
# assign 0.3/year as the maximum fishing mortality to demersal fish
p=setFishing(p=p,F=0.3,etaF=0.05,groupidx=c(5))
df=cbind(df,data.frame(mortF_new=c(p$mortF[p$ix[[1]]],p$mortF[p$ix[[5]]])))
knitr::kable(df,caption="Fishing mortality before and after assignment")
df=data.frame("Stage"=1:length(p$ix[[1]]), "mortF"=p$mortF[p$ix[[1]]],"Groups"="smallPel")
df=rbind(df,data.frame("Stage"=1:length(p$ix[[5]]), 
                       "mortF"=p$mortF[p$ix[[5]]],"Groups"="demersals"))
df$Groups=factor(df$Groups,levels=c("smallPel","demersals"))
# plot of fishing mortality of small pelagics and demersals
fig=ggplot(df, aes(x = Stage, y = mortF, color = Groups))+
    geom_line(linewidth = 0.7,alpha=0.9)+
    geom_point(size=1.5,alpha=0.9)+
  labs(x = expression("Stage"), y = expression("Fishing mortality"~(yr^{-1}))) +
  scale_color_manual(values = c("red", "black"),labels=c("Small pelagics","Demersals")) +
  scale_x_continuous(breaks = unique(df$Stage))+
      theme(panel.background = element_rect(fill = "white"),
          panel.border = element_rect(color = "black", fill = NA),
          axis.line = element_line(color = "black"),
          #legend.title = element_blank(),
          legend.key = element_rect(fill = "transparent", color = "transparent"),
          legend.position = "bottom")
  
fig


## -----------------------------------------------------------------------------
sim=simulateFEISTY(p = p, tEnd = 200, USEdll = T, bCust = T)

## -----------------------------------------------------------------------------
plotYieldtime(sim)

## -----------------------------------------------------------------------------
plotSSBtime(sim)

## -----------------------------------------------------------------------------
p1=setupVertical2(depth=1000,szprod=5, lzprod=5,dfpho = 130) # oligotrophic 1000 meter
p2=setupVertical2(depth=1000,szprod=100, lzprod=100,dfpho =380) # eutrophic 1000 meter

## ----fig.width=8, fig.height=7------------------------------------------------
sim1=simulateFEISTY(p=p1,tEnd=500)
plotSimulation(sim1)

## ----fig.width=8, fig.height=7------------------------------------------------
sim2=simulateFEISTY(p=p2,tEnd=500)
plotSimulation(sim2)

