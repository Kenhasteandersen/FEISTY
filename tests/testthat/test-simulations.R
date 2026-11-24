#
# Free runs of the four setups
#

for (n in names(setups)) local({
  setup <- setups[[n]]

  test_that(paste(n, "R and Fortran simulations agree, without and with fishing"), {
    for (Fmax in c(0, 0.3)) {
      p <- setFishing(setup(), Fmax)
      fortran <- simulateFEISTY(p = p, tEnd = 50)
      r       <- simulateFEISTY(p = p, tEnd = 50, USEdll = FALSE)
      for (x in c("u", "totBiomass", "SSB", "yield", "mortpred", "f"))
        expect_close(fortran[[x]], r[[x]], 1e-5, sprintf("%s F=%g %s", n, Fmax, x))
    }
  })

  test_that(paste(n, "output is well-formed"), {
    p <- setup()
    sim <- simulateFEISTY(p = p, tEnd = 50)
    expect_s3_class(sim, "FEISTY")
    expect_equal(sim$nTime, 51)
    expect_equal(dim(sim$B), c(51, length(p$ixFish)))
    expect_equal(dim(sim$R), c(51, 4))
    expect_equal(dim(sim$totBiomass), c(51, p$nGroups))
    expect_true(all(is.finite(sim$u)))
    expect_true(min(sim$u) > -1e-6)
    expect_true(all(sim$f >= 0 & sim$f < 1))
    for (g in seq_len(p$nGroups))
      expect_equal(unname(sim$totBiomass[, g]), unname(rowSums(sim$u[, p$ix[[g]]])), tolerance = 1e-8)
    expect_equal(unname(sim$totRecruit),
                 unname(sweep(sim$totRepro, 2, rep(p$epsRepro, length.out = p$nGroups), "*")),
                 tolerance = 1e-10)
  })

  test_that(paste(n, "does not depend on previous simulations"), {
    p <- setup()
    first <- simulateFEISTY(p = p, tEnd = 20)$u
    for (m in names(setups)) simulateFEISTY(p = setFishing(setups[[m]](), 1), tEnd = 5)
    runQuiet(p = makeTS(setupVertical2(), fishing = TRUE))
    expect_identical(simulateFEISTY(p = p, tEnd = 20)$u, first)
  })
})

test_that("the initial conditions can be changed", {
  p <- setupVertical2()
  yini <- p$u0*2
  expect_equal(unname(simulateFEISTY(p = p, tEnd = 1, yini = yini)$u[1, ]), unname(yini))
  expect_error(simulateFEISTY(p = p, yini = 1:3), "yini")
})

test_that("too many stages are refused", {
  expect_error(simulateFEISTY(p = setupVertical2(nStages = 30), tEnd = 1), "27")
})

test_that("biomass and yield agree with the reference", {
  ref <- utils::read.csv(test_path("fixtures", "reference.csv"))
  cases <- referenceCases()
  expect_setequal(unique(ref$case), names(cases))
  for (k in names(cases)) {
    s <- simSummary(cases[[k]]())
    r <- ref[ref$case == k, ]
    expect_equal(s$biomass, r$biomass, tolerance = 1e-3, info = paste(k, "biomass"))
    expect_equal(s$yield,   r$yield,   tolerance = 1e-3, info = paste(k, "yield"))
  }
})
