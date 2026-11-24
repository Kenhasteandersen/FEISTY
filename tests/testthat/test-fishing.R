#
# Fishing: setFishing(), calcYield() and the effect on the simulations
#

for (n in names(setups)) local({
  setup <- setups[[n]]

  test_that(paste(n, "setFishing makes a trawl selectivity"), {
    p <- setFishing(setup(), 0.4, etaF = 0.1)
    expect_equal(unname(p$mortF[p$ixR]), rep(0, 4))
    for (g in seq_len(p$nGroups)) {
      ix <- p$ix[[g]]
      psi <- 1/(1 + (p$mc[ix]/(0.1*max(p$mUpper[ix])))^(-3))
      expect_equal(unname(p$mortF[ix]), unname(0.4*psi))
      expect_false(is.unsorted(p$mortF[ix]))
      expect_true(all(p$mortF[ix] <= 0.4))
    }
  })

  test_that(paste(n, "setFishing only changes the selected groups"), {
    p0 <- setFishing(setup(), 0.1)
    for (g in seq_len(p0$nGroups)) {
      p <- setFishing(p0, 0.7, groupidx = g)
      others <- setdiff(p0$ixFish, p0$ix[[g]])
      expect_identical(p$mortF[others], p0$mortF[others])
      expect_true(all(p$mortF[p$ix[[g]]] > p0$mortF[p$ix[[g]]]))
    }
    expect_equal(setFishing(setFishing(setup(), 0.5), 0)$mortF[p0$ixFish], 0*p0$mortF[p0$ixFish])
  })

  test_that(paste(n, "yield is fishing mortality times biomass"), {
    p <- setFishing(setup(), seq(0.1, 0.5, length.out = setup()$nGroups))
    sim <- simulateFEISTY(p = p, tEnd = 50)
    for (g in seq_len(p$nGroups)) {
      ix <- p$ix[[g]]
      expect_equal(unname(sim$yield[, g]), as.vector(pmax(sim$u[, ix], 0) %*% p$mortF[ix]),
                   tolerance = 1e-10)
    }
    last <- sim$t >= 0.6*max(sim$t)
    expect_equal(unname(sim$yieldAMean), unname(colMeans(sim$yield[last, ])))
    expect_true(all(sim$yieldMin <= sim$yieldAMean & sim$yieldAMean <= sim$yieldMax))
  })

  test_that(paste(n, "no fishing gives no yield"), {
    expect_true(all(simulateFEISTY(p = setup(), tEnd = 10)$yield == 0))
  })

  test_that(paste(n, "fishing a group reduces its spawning stock biomass"), {
    p <- setup()
    last <- function(sim, x) sim[[x]][sim$nTime, ]
    free <- last(simulateFEISTY(p = p, tEnd = 200), "SSB")
    for (g in seq_len(p$nGroups)) {
      if (free[g] < 1e-3) next   # group not present
      sim <- simulateFEISTY(p = setFishing(p, 0.5, groupidx = g), tEnd = 200)
      expect_lt(last(sim, "SSB")[[g]], free[[g]])
      expect_gt(last(sim, "yield")[[g]], 0)
    }
  })
})

test_that("Fmax can be given to the setup or to setFishing", {
  for (n in c("setupBasic2", "setupVertical2")) {
    a <- setups[[n]](Fmax = 0.4, etaF = 0.1)
    b <- setFishing(setups[[n]](), 0.4, etaF = 0.1)
    expect_equal(a$mortF, b$mortF)
  }
})

test_that("setFishing checks its input", {
  p <- setupVertical2()
  expect_error(setFishing(p, 0.1, groupidx = 6), "groupidx")
  expect_error(setFishing(p, 0.1, groupidx = 0), "groupidx")
  expect_error(setFishing(p, c(0.1, 0.2)), "Fmax")
  expect_error(setFishing(p, 0.1, etaF = c(0.1, 0.2)), "etaF")
  expect_error(setFishing(p, NA), "NA")
  expect_error(setFishing(p, 0.1, etaF = NA), "NA")
  expect_identical(setFishing(p, 0.1, groupidx = integer(0)), p)
})
