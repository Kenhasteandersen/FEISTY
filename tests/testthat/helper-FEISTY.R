#
# Shared helpers for the FEISTY tests. Sourced automatically by testthat.
#

setups <- list(setupBasic     = setupBasic,
               setupBasic2    = setupBasic2,
               setupVertical  = setupVertical,
               setupVertical2 = setupVertical2)

# Largest difference relative to the largest reference value
relDiff <- function(x, ref) {
  x <- as.vector(x); ref <- as.vector(ref)
  if (length(x) != length(ref)) return(Inf)
  if (!identical(is.finite(x), is.finite(ref))) return(Inf)
  ok <- is.finite(ref)
  max(abs(x[ok] - ref[ok])) / max(abs(ref[ok]), 1e-300)
}

expect_close <- function(x, ref, tol, label) {
  d <- relDiff(x, ref)
  expect_true(d <= tol, label = sprintf("%s: relative difference %.2e <= %.0e", label, d, tol))
}

# State vectors in which to compare derivatives. The spare resource (no. 4, with K=0) is kept
# at 0, as in all simulations with the prescribed setups.
testStates <- function(p) {
  set.seed(1)
  n <- length(p$u0)
  fish <- p$ixFish
  states <- list(
    initial  = p$u0,
    random   = p$u0*runif(n, 0.1, 3) + runif(n, 0, 0.5),
    small    = rep(1e-8, n),
    extinct  = replace(p$u0 + 0.1, fish[seq(1, length(fish), by = 2)], 0),
    noFish   = replace(p$u0, fish, 0),
    negative = replace(p$u0 + 0.1, fish[1], -1e-6))
  lapply(states, function(u) replace(u, 4, 0))
}

# Derivatives and all other outputs from R or Fortran, as one list of plain vectors
derivsR <- function(p, u) {
  lapply(simulateFEISTY(p = p, times = NA, yini = u, USEdll = FALSE), as.vector)
}

derivsFortran <- function(p, u, bCust = TRUE) {
  out <- simulateFEISTY(p = p, times = NA, yini = u, USEdll = TRUE, bCust = bCust)
  nF <- length(p$ixFish); nS <- p$nStages[1]; nG <- p$nGroups[1]
  len <- c(f=nF, mortpred=nS, g=nF, Repro=nF, Fin=nF, Fout=nF, totMort=nG, totGrazing=nG,
           totLoss=nG, totRepro=nG, totRecruit=nG, totBiomass=nG)
  stopifnot(sum(len) == length(out[[2]]))
  c(list(deriv = as.vector(out[[1]])),
    lapply(split(as.vector(out[[2]]), factor(rep(names(len), len), levels = names(len))), as.vector))
}

# The Fortran and R derivatives, and all other outputs, must be identical to round-off
expect_same_derivs <- function(p, label, bCust = TRUE, tol = 1e-12) {
  states <- testStates(p)
  for (s in names(states)) {
    r <- derivsR(p, states[[s]])
    f <- derivsFortran(p, states[[s]], bCust = bCust)
    for (x in names(f))
      expect_close(f[[x]], r[[x]], tol, sprintf("%s, state '%s', %s (bCust=%s)", label, s, x, bCust))
  }
}

# Known issues: "ok" tells whether the code behaves correctly. Skips (visibly) as long as it does
# not, and fails when the issue has been fixed, so that it can be made into a normal test.
known_issue <- function(ok, description) {
  if (isTRUE(ok))
    fail(paste("Known issue seems fixed; make this a normal test:", description))
  else
    skip(paste("KNOWN ISSUE:", description))
}

# ---- Time series -------------------------------------------------------------------------------

tsdata <- local({
  e <- new.env()
  utils::data("tsinput_example_1850_2014", package = "FEISTY", envir = e)
  as.list(e)
})

# Time-series parameters from the example data.
#  fishing: FALSE or TRUE; ix: indices of months to use.
makeTS <- function(p, fishing = FALSE, ix = 1501:1524, tSpin = NA, nSpinloop = 2,
                   bConstantT = FALSE, randomStart = TRUE) {
  d <- tsdata
  args <- list(p, tStep_ts = 1/12, tSpin = tSpin, nSpinloop = nSpinloop,
               szbio_ts  = d$Zbio[ix]/2,   lzbio_ts  = d$Zbio[ix]/2,
               szprod_ts = d$Zhploss[ix]/2, lzprod_ts = d$Zhploss[ix]/2,
               dfbot_ts  = d$dfbot[ix])
  if (!bConstantT) {
    args$Tp_ts <- d$Tp[ix]
    args$Tb_ts <- d$Tb[ix]
    if (p$setup == "setupVertical2") args$Tm_ts <- d$Tm[ix]
  }
  if (fishing) {
    args <- c(args, list(Fsmp_ts = d$Fspel[ix], Flgp_ts = d$Flpel[ix], Fdem_ts = d$Fdem[ix]))
    if (p$setup == "setupVertical2")
      args <- c(args, list(Fmesop_ts = d$Fmeso[ix], Fmidwp_ts = d$FmidP[ix]))
  }
  utils::capture.output(p <- do.call(setupTimeseries, args))
  if (randomStart) {
    set.seed(2)
    p$u0[p$ixFish] <- runif(length(p$ixFish), 0.01, 1)
  }
  p
}

runQuiet <- function(...) {
  utils::capture.output(sim <- simulateFEISTY(...))
  sim
}

# ---- Reference values --------------------------------------------------------------------------

# Simulations whose biomass and yield are compared with fixtures/reference.csv.
# After an intended change of the model, update with tests/update_reference.R
readCells <- function() utils::read.csv(test_path("fixtures", "global_cells.csv"))

referenceCases <- function() {
  cases <- list()
  for (n in names(setups)) local({
    s <- setups[[n]]
    cases[[paste0(n, "/free")]]    <<- function() simulateFEISTY(p = s(), tEnd = 200)
    cases[[paste0(n, "/fishing")]] <<- function() simulateFEISTY(p = setFishing(s(), 0.3), tEnd = 200)
  })
  for (n in c("setupBasic2", "setupVertical2")) local({
    s <- setups[[n]]
    cases[[paste0(n, "/27stages")]] <<- function() simulateFEISTY(p = s(nStages = 27), tEnd = 100)
  })
  cases[["timeseries/setupBasic2"]] <- function()
    runQuiet(p = makeTS(setupBasic2(depth = 500), fishing = TRUE, ix = 1501:1620))
  cases[["timeseries/setupVertical2"]] <- function()
    runQuiet(p = makeTS(setupVertical2(), fishing = TRUE, ix = 1501:1620))
  cases[["timeseries/spinup"]] <- function()
    runQuiet(p = makeTS(setupVertical2(), fishing = TRUE, ix = 1501:1620, tSpin = 5))
  cases
}

# Mean over the last 40% of the simulation
simSummary <- function(sim) {
  ix <- sim$t >= 0.6*max(sim$t)
  data.frame(group   = seq_len(sim$p$nGroups[1]),
             biomass = unname(colMeans(sim$totBiomass[ix, , drop = FALSE])),
             yield   = unname(colMeans(sim$yield[ix, , drop = FALSE])))
}

# ---- Global runs -------------------------------------------------------------------------------

# Fished groups in the global tests (1: small pel., 3: large pel., 5: demersals)
simulateCell <- function(cell, Fmax = 0, ixGroups = NULL, tEnd = 100) {
  simulatePosition(setupVertical2, cell$lat, cell$lon, Fmax = Fmax, ixGroups = ixGroups,
                   tEnd = tEnd, glob = cell)
}
