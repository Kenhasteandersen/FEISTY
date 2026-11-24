#
# The derivatives, and all other outputs, must be the same from R and Fortran.
#

for (n in names(setups)) local({
  setup <- setups[[n]]

  test_that(paste(n, "R and Fortran derivatives are identical"), {
    expect_same_derivs(setup(), n)
  })

  test_that(paste(n, "R and Fortran derivatives are identical with fishing"), {
    p <- setup()
    expect_same_derivs(setFishing(p, 0.5), paste(n, "all groups"))
    expect_same_derivs(setFishing(p, 0.5, etaF = 0.2), paste(n, "etaF=0.2"))
    for (g in seq_len(p$nGroups))
      expect_same_derivs(setFishing(p, 0.5, groupidx = g), paste(n, "group", g))
    expect_same_derivs(setFishing(p, seq(0.1, 1, length.out = p$nGroups)), paste(n, "F per group"))
  })

  test_that(paste(n, "R and Fortran derivatives are identical at different depths"), {
    for (depth in c(60, 199, 201, 1000, 4000))
      expect_same_derivs(setFishing(setup(depth = depth), 0.3), paste(n, "depth", depth))
  })

  test_that(paste(n, "R and Fortran derivatives are identical for different productions"), {
    for (prod in c(1, 500))
      expect_same_derivs(setup(szprod = prod, lzprod = prod/2, dfbot = prod), paste(n, "prod", prod))
  })

  test_that(paste(n, "the built-in Fortran setup is identical to R"), {
    expect_same_derivs(setup(), n, bCust = FALSE)
  })
})

for (n in c("setupBasic2", "setupVertical2")) local({
  setup <- setups[[n]]
  test_that(paste(n, "R and Fortran derivatives are identical for different numbers of stages"), {
    for (nStages in c(6, 12, 21, 27)) {
      p <- setup(nStages = nStages, Fmax = 0.3)
      expect_same_derivs(p, paste(n, nStages, "stages"))
      expect_same_derivs(p, paste(n, nStages, "stages"), bCust = FALSE)
    }
  })

  test_that(paste(n, "the built-in Fortran setup is identical to R with fishing"), {
    expect_same_derivs(setup(Fmax = 0.5), n, bCust = FALSE)
    expect_same_derivs(setup(Fmax = 0.5, etaF = 0.2), n, bCust = FALSE)
    expect_same_derivs(setFishing(setup(), 0.5), n, bCust = FALSE)
  })
})

for (n in c("setupBasic", "setupVertical")) local({
  setup <- setups[[n]]
  test_that(paste(n, "the built-in Fortran setup warns that fishing is ignored"), {
    p <- setFishing(setup(), 0.5)
    expect_warning(simulateFEISTY(p = p, tEnd = 2, bCust = FALSE), "Fishing is ignored")
    expect_warning(simulateFEISTY(p = p, times = NA, bCust = FALSE), "Fishing is ignored")
    # No warning if fishing is included, or if there is no fishing:
    expect_no_warning(simulateFEISTY(p = p, tEnd = 2, bCust = TRUE))
    expect_no_warning(simulateFEISTY(p = p, tEnd = 2, USEdll = FALSE))
    expect_no_warning(simulateFEISTY(p = setup(), tEnd = 2, bCust = FALSE))
    expect_no_warning(simulateFEISTY(p = setFishing(p, 0), tEnd = 2, bCust = FALSE))
  })
})

test_that("temperature and bET give identical derivatives", {
  for (Tp in c(-1, 15, 28)) {
    expect_same_derivs(setupBasic(Tp = Tp, Tb = Tp/2), paste("setupBasic Tp", Tp))
    expect_same_derivs(setupBasic2(Tp = Tp, Tb = Tp/2, bET = TRUE), paste("setupBasic2 bET Tp", Tp))
    expect_same_derivs(setupBasic2(Tp = Tp, Tb = Tp/2, bET = FALSE), paste("setupBasic2 no bET Tp", Tp))
    expect_same_derivs(setupVertical2(Tp = Tp, Tm = Tp*0.6, Tb = Tp/2), paste("setupVertical2 Tp", Tp))
  }
  for (region in 1:4)
    expect_same_derivs(setupVertical(region = region), paste("setupVertical region", region))
})

test_that("group-specific fishing is refused by the built-in Fortran setups", {
  for (n in names(setups))
    expect_error(simulateFEISTY(p = setFishing(setups[[n]](), 0.5, groupidx = 1), tEnd = 2, bCust = FALSE),
                 "not supported")
})
