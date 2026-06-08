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

## ----fig.width = 8, fig.height = 6, cache = TRUE------------------------------
# load example data from package 
data(tsinput_example_1850_2014)
# setup time series based on setupBasic2
# 10-year time-series simulation after 2-year spinning up (loop 4 times) 
p = setupTimeseries(p = setupBasic2(depth = depth, nStages = 9),
                    tStep_ts = 1/12,
                    tSpin = 2,
                    nSpinloop = 4, 
                    Tp_ts = Tp[1:120],
                    Tm_ts = Tm[1:120],
                    Tb_ts = Tb[1:120],
                    szbio = Zbio[1:120]/2,
                    lzbio = Zbio[1:120]/2,
                    szprod_ts = Zhploss[1:120]/2,
                    lzprod_ts = Zhploss[1:120]/2,
                    dfbot_ts  = dfbot[1:120],
                    Fsmp_ts   = Fspel[1:120],
                    Flgp_ts   = Flpel[1:120],
                    Fdem_ts   = Fdem[1:120])
# Run simulation by Fortran
simF = simulateFEISTY(p = p, USEdll = T)

# Set dateframe for plotting
diagts <- data.frame(group = factor(rep(c("large zooplankton production",
                                          "large zooplankton consumption",
                                          "down-regulated large zooplankton consumption"),
                                           each = length(simF$t)),
                              levels = c("large zooplankton production",
                                         "large zooplankton consumption",
                                         "down-regulated large zooplankton consumption")), 
                     val    = c(as.numeric(p$lzprod_ts),
                                as.numeric(simF$lgzcsp),
                                as.numeric(simF$lgzcsp_dr)),
                     t      = rep(simF$t,3))
# Plot
ggplot(data = diagts, aes(x = t, y = val, color = group, group = group)) +
    geom_line(linewidth = 1) +  # Plot lines
    xlab("Time (yr)") +
    ylab(expression("production and consumption (g m"^"-2"*" yr"^"-1"*")")) +
    annotation_logticks(sides = "l", linewidth = 0.4, colour = "darkgrey") +
    coord_cartesian(ylim = c(1E-2, max(1E-2 * 100, max(diagts$val) * 5))) +
    #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
    #              labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(
    breaks = scales::breaks_log(n = 5, base = 10),
    labels = scales::label_log(base = 10)
  ) +
    theme_classic()+
    theme(legend.key = element_blank())+
    theme(legend.position = "inside",
          legend.position.inside = c(0.8, 0.2))+
    labs(color = NULL)


## ----fig.width = 8, fig.height = 6, cache = TRUE------------------------------
# load example data from package 
data(tsinput_example_1850_2014)
# setup time series based on setupVertical2 (9 stages)
# 165-year time-series simulation after 16.5-year spinning up (loop 4 times)
p = setupTimeseries(p = setupVertical2(photic = photic, depth = depth, nStages = 9),
                  tStep_ts = 1/12,
                  tSpin = 16.5, 
                  nSpinloop = 4,
                  Tp_ts = Tp,
                  Tm_ts = Tm,
                  Tb_ts = Tb,
                  szbio = Zbio/2,
                  lzbio = Zbio/2,
                  szprod_ts = Zhploss/2,
                  lzprod_ts = Zhploss/2,
                  dfbot_ts  = dfbot,
                  Fsmp_ts   = Fspel,
                  Fmesop_ts = Fmeso,
                  Flgp_ts   = Flpel,
                  Fmidwp_ts = FmidP,
                  Fdem_ts   = Fdem)
# Run simulation by Fortran
simF = simulateFEISTY(p = p, USEdll = T)
simR = simulateFEISTY(p = p, USEdll = F)

# Plot
plot(simR$t, rowSums(simR$totBiomass,2) ,type='l', xlab="year",ylab="tot biomass",
      log="", ylim=c(0.1,40), xlim=c(0,165), col='red')
lines(simF$t, rowSums(simF$totBiomass,2), type='l', col='black')

# setup time series based on setupVertical2 (3 stages)
# 165-year time-series simulation after 16.5-year spinning up (loop 4 times)
p = setupTimeseries(p = setupVertical2(photic = photic, depth = depth, nStages = 3),
                  tStep_ts = 1/12,
                  tSpin = 16.5, 
                  nSpinloop = 4,
                  Tp_ts = Tp,
                  Tm_ts = Tm,
                  Tb_ts = Tb,
                  szbio = Zbio/2,
                  lzbio = Zbio/2,
                  szprod_ts = Zhploss/2,
                  lzprod_ts = Zhploss/2,
                  dfbot_ts  = dfbot,
                  Fsmp_ts   = Fspel,
                  Fmesop_ts = Fmeso,
                  Flgp_ts   = Flpel,
                  Fmidwp_ts = FmidP,
                  Fdem_ts   = Fdem)
# Run simulation by Fortran
simF = simulateFEISTY(p = p, USEdll = T)
simR = simulateFEISTY(p = p, USEdll = F)

# Plot
plot(simR$t, rowSums(simR$totBiomass,2) ,type='l', xlab="year",ylab="tot biomass",
      log="", ylim=c(0.1,40), xlim=c(0,165), col='red')
lines(simF$t, rowSums(simF$totBiomass,2), type='l', col='black')



## ----fig.width = 8, fig.height = 6, cache = TRUE------------------------------
# load example data from package 
data(tsinput_example_1850_2014)
# add a perturbation
Zbio[100]
Zbio[100]=Zbio[100]+1E-5
# setup time series based on setupVertical2 (3 stages)
# 165-year time-series simulation after 16.5-year spinning up (loop 4 times)
p = setupTimeseries(p = setupVertical2(photic = photic, depth = depth, nStages = 3),
                  tStep_ts = 1/12,
                  tSpin = 16.5, 
                  nSpinloop = 4,
                  Tp_ts = Tp,
                  Tm_ts = Tm,
                  Tb_ts = Tb,
                  szbio = Zbio/2,
                  lzbio = Zbio/2,
                  szprod_ts = Zhploss/2,
                  lzprod_ts = Zhploss/2,
                  dfbot_ts  = dfbot,
                  Fsmp_ts   = Fspel,
                  Fmesop_ts = Fmeso,
                  Flgp_ts   = Flpel,
                  Fmidwp_ts = FmidP,
                  Fdem_ts   = Fdem)
# Run simulation by Fortran
simF2 = simulateFEISTY(p = p, USEdll = T)
# Plot
plot(simF$t, rowSums(simF$totBiomass,2) ,type='l', xlab="year",ylab="tot biomass",
      log="", ylim=c(0.1,40), xlim=c(0,165), col='green')
lines(simF2$t, rowSums(simF2$totBiomass,2), type='l', col='black')


