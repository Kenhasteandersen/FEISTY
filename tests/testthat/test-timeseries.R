#
# Time-series runs
#

outputs <- c("f", "mortpred", "g", "Repro", "Fin", "Fout", "totMort", "totGrazing", "totLoss",
             "totRepro", "totRecruit", "totBiomass", "smzcsp", "smzcsp_dr", "lgzcsp", "lgzcsp_dr", "yield")

# Compare with R at the first time, where the states are the same (identical to round-off), and 
# for the whole simulation (accuracy of the solvers)
compareTS <- function(p, label, tolSim = 2e-3) {
  fortran <- runQuiet(p = p, yini = p$u0)
  r       <- runQuiet(p = p, yini = p$u0, USEdll = FALSE)
  first <- function(sim, x) as.matrix(sim[[x]])[1, ]
  for (x in outputs)
    expect_close(first(fortran, x), first(r, x), 1e-10, paste(label, "first time,", x))
  for (x in c("u", "totBiomass", "yield", "smzcsp_dr", "mortpred"))
    expect_close(fortran[[x]], r[[x]], tolSim, paste(label, "simulation,", x))
  expect_close(fortran$p$mortF, r$p$mortF, 1e-10, paste(label, "mortF"))
}

tsSetups <- list("setupBasic2 no bET" = function() setupBasic2(bET = FALSE),
                 "setupBasic2 deep"   = function() setupBasic2(depth = 500),
                 "setupVertical2"     = function() setupVertical2(),
                 "setupVertical2 shelf" = function() setupVertical2(depth = 100, photic = 80),
                 "setupVertical2 12 stages" = function() setupVertical2(nStages = 12))

for (n in names(tsSetups)) local({
  setup <- tsSetups[[n]]
  test_that(paste("time series: R and Fortran agree for", n), {
    compareTS(makeTS(setup()), n)
  })
  test_that(paste("time series: R and Fortran agree with fishing for", n), {
    compareTS(makeTS(setup(), fishing = TRUE), paste(n, "fishing"))
  })
})

test_that("time series: R and Fortran agree for setupBasic", {
  compareTS(makeTS(setupBasic(depth = 500)), "setupBasic deep")
})

test_that("time series: R and Fortran agree with constant temperature and bET", {
  compareTS(makeTS(setupBasic(),  bConstantT = TRUE), "setupBasic constant T")
  compareTS(makeTS(setupBasic2(), bConstantT = TRUE), "setupBasic2 constant T")
})

test_that("time series: R and Fortran agree with temperature series and bET", {
  p <- makeTS(setupBasic2())
  g <- function(...) runQuiet(p = p, yini = p$u0, ...)$g[1, ]
  known_issue(relDiff(g(), g(USEdll = FALSE)) < 1e-10,
              "Fortran uses constant temperatures in updateET() (bET, depth < 200 m) in time-series runs")
})

test_that("time series: setupBasic with fishing", {
  ok <- !inherits(try(runQuiet(p = makeTS(setupBasic(), fishing = TRUE)), silent = TRUE), "try-error")
  known_issue(ok, "setupBasic does not define etaF, so time series of fishing fail")
})

test_that("time series: fishing varies in time and gives the yield", {
  ix <- 1501:1524
  p <- makeTS(setupVertical2(), fishing = TRUE, ix = ix)
  sim <- runQuiet(p = p)
  Fts <- cbind(tsdata$Fspel, tsdata$Fmeso, tsdata$Flpel, tsdata$FmidP, tsdata$Fdem)[c(ix, max(ix)), ]
  for (g in 1:5) {
    psi <- setFishing(p, 1, etaF = 0.05, groupidx = g)$mortF[p$ix[[g]]]
    expect_equal(unname(sim$p$mortF[, p$ix[[g]]]), unname(outer(Fts[, g], psi)), tolerance = 1e-10)
    expect_equal(unname(sim$yield[, g]), 
                 unname(rowSums(sim$p$mortF[, p$ix[[g]]]*pmax(sim$u[, p$ix[[g]]], 0))), tolerance = 1e-10)
  }
  # Fishing reduces the total biomass of the main fished groups:
  free <- runQuiet(p = makeTS(setupVertical2(), ix = ix))
  expect_lt(sim$totBiomass[25, 3], free$totBiomass[25, 3])
})

test_that("time series: zooplankton follows the input and consumption is limited by production", {
  ix <- 1501:1524
  p <- makeTS(setupVertical2(), ix = ix)
  p$u0[p$ixFish] <- 20   # many fish, so that the limitation is needed
  sim <- runQuiet(p = p, yini = p$u0)
  i <- c(ix, max(ix))
  expect_equal(unname(sim$R[, 1]), tsdata$Zbio[i]/2)
  expect_equal(unname(sim$R[, 2]), tsdata$Zbio[i]/2)
  prod <- tsdata$Zhploss[i]/2
  expect_true(all(sim$smzcsp_dr <= prod*(1 + 1e-10)))
  expect_true(all(sim$lgzcsp_dr <= prod*(1 + 1e-10)))
  expect_true(any(sim$smzcsp > prod | sim$lgzcsp > prod))
})

test_that("time series: no zooplankton production", {
  p <- makeTS(setupVertical2())
  p$szprod_ts <- p$szprod_ts*0
  p$lzprod_ts <- p$lzprod_ts*0
  sim <- runQuiet(p = p)
  expect_true(all(is.finite(sim$u)))
  expect_true(all(sim$smzcsp_dr == 0))
})

test_that("time series: spin-up is used as the initial condition", {
  ix <- 1501:1620
  p <- makeTS(setupVertical2(), fishing = TRUE, ix = ix, tSpin = 5, nSpinloop = 2)
  sim <- runQuiet(p = p)
  # One spin-up of 5 years, repeated by hand:
  y <- p$u0
  for (i in 1:2) {
    q <- makeTS(setupVertical2(), fishing = TRUE, ix = ix[1:60])
    u <- runQuiet(p = q, yini = y)$u
    y <- u[nrow(u), ]
  }
  expect_close(sim$u[1, p$ixFish], y[p$ixFish], 1e-3, "spin-up")
  expect_error(runQuiet(p = makeTS(setupVertical2(), tSpin = 1)), "longer than 10 years")
})

test_that("time series: repeated and interleaved runs are identical", {
  a <- makeTS(setupVertical2(), fishing = TRUE)
  b <- makeTS(setupBasic2(depth = 500), ix = 1201:1260)
  first <- runQuiet(p = a)$u
  runQuiet(p = b)
  simulateFEISTY(p = setupBasic(), tEnd = 5)
  expect_identical(runQuiet(p = a)$u, first)
})

test_that("time series: the input is checked", {
  d <- tsdata; i <- 1:24
  ts <- function(...) {
    args <- list(p = setupVertical2(), tSpin = NA, szbio_ts = d$Zbio[i], lzbio_ts = d$Zbio[i],
                 szprod_ts = d$Zhploss[i], lzprod_ts = d$Zhploss[i])
    new <- list(...)
    args[names(new)] <- new
    utils::capture.output(p <- do.call(setupTimeseries, args))
    p
  }
  expect_true(ts()$bTS)
  expect_equal(ts()$tEnd_ts, 2)
  expect_error(ts(szbio_ts = d$Zbio[1:12]), "same length")
  expect_error(ts(szbio_ts = replace(d$Zbio[i], 3, NA)), "NA values")
  expect_error(ts(szbio_ts = NA), "must be provided")
  expect_error(ts(bprodin_ts = d$dfbot[i], dfbot_ts = d$dfbot[i]), "Only one")
  expect_error(ts(p = setupBasic2(), Fmesop_ts = d$Fmeso[i]), "mesopelagic")
})
